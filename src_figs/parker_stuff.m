function varargout = parker_stuff(varargin)
% Helper window to do Parker inversions and Reduce To the Pole 

%	Copyright (c) 2004-2016 by J. Luis
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

	hObject = figure('Vis','off');
	parker_stuff_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hObject,'center')

	handles.hCallingFig = [];		% Handles to the calling figure
	handles.geog = 0;				% Set this as default
	handles.is_redPole = 0;			% Flag to signal that options are for use in reduction to the pole
	handles.Z_bat = [];
	handles.Z_src = [];
	handles.zobs = 0;				% Default observation level
	handles.data = [];				% Date has no default value
	no_igrf = false;				% In reduction to the pole we don't need to compute the IGRF

	if (~isempty(varargin))
		if (strcmp(varargin{1},'parker_direct'))
			h_txt = findobj(hObject,'Tag','text_FieldMag');
			set(h_txt,'String','Mag')
			set(handles.edit_SourceGrid,'Tooltip','Enter Magnetization grid name')
			set(handles.edit_wlong,'Enable','off','BackgroundColor',get(0,'DefaultUicontrolBackgroundColor'))
			set(handles.edit_wshort,'Enable','off','BackgroundColor',get(0,'DefaultUicontrolBackgroundColor'))
			set(handles.edit_sdec,'String','0')			% Default is geocentric dipole
			set(handles.edit_sdip,'String','0')
			set(handles.check_CenterDipole,'Value',1)
			set(handles.check_fieldIsRTP,'Visible','off')
			set(hObject,'Name','Parker Direct')
			set([handles.check_showResidue handles.check_showDirect], 'Vis', 'off')
			handles.what_parker = 'direct';
			handles.hCallingFig = varargin{2};
		elseif (strcmp(varargin{1},'parker_inverse'))
			set(handles.edit_sdec,'String','0')			% Default is geocentric dipole
			set(handles.edit_sdip,'String','0')
			set(handles.check_CenterDipole,'Value',1)
			set(handles.check_fieldIsRTP,'Value',0)
			set(hObject,'Name','Parker Inverse')
			handles.what_parker = 'inverse';        
			handles.hCallingFig = varargin{2};
		elseif (strcmp(varargin{1},'redPole'))
			no_igrf = true;
			set(handles.edit_BatGrid,'Enable','off','BackgroundColor',get(0,'DefaultUicontrolBackgroundColor'))
			set(handles.push_BatGrid,'Enable','inactive')
			set(handles.edit_wlong,'Enable','off','BackgroundColor',get(0,'DefaultUicontrolBackgroundColor'))
			set(handles.edit_wshort,'Enable','off','BackgroundColor',get(0,'DefaultUicontrolBackgroundColor'))
			set(handles.edit_date,'Enable','off','BackgroundColor',get(0,'DefaultUicontrolBackgroundColor'))
			set([handles.check_CenterDipole handles.check_fieldIsRTP],'Visible','off')
			set([handles.text_wshort handles.text_wlong],'Enable','off')
			set([handles.check_showResidue handles.check_showDirect], 'Vis', 'off')
			% Reuse some edit boxes, but we have to change their text
			set(handles.text_Level,'String','Field dip')
			set(handles.edit_zobs,'String','','Tooltip','Inclination of the magnetic field.')
			set(handles.text_Thickness,'String','Field dec')
			set(handles.edit_thickness,'String','','Tooltip','Declination of the magnetic field.')
			set(hObject,'Name','Reduce to the Pole')
			handles.what_parker = 'redPole';
			handles.hCallingFig = varargin{2};
			
			handMir = guidata(handles.hCallingFig);
			if (handMir.validGrid)
				[X, Y, Z, head] = load_grd(handMir);
				if (~isempty(Z) && head(5) < 0 && head(6) > 0)	% Run some simple tests
					handles.X = X;	handles.Y = Y;	handles.Z_src = Z;	handles.head_src = head;
					handles.have_nans = handMir.have_nans;
					push_SourceGrid_CB(handles.push_SourceGrid, handles, '__inMEM')
					set(handles.edit_SourceGrid,'String',' In memory array','Enable','off')
					handles = guidata(hObject);		% Get the updated version
				end
			end
			
		else				% defaults to "direct"
			handles.what_parker = 'direct';
		end
	else
		% Else defaults to "direct" but without knowing hCallingFig
		handles.what_parker = 'direct';
	end

	if (~no_igrf)
		handles.start_stop_epoch = [1900 2015];     % ISTO TEM DE SER AUTOMATIZADO
	end

	% Set upt some useful tooltips
	str = sprintf(['The default value is the number of rows in the grid\n',...
		'However, for reducing border effects you may want to apply\n',...
		'a skirt to the grid. For that, select a value from the side\n',...
		'listbox. Extra points will be padded by mirroiong the west side.']);
	set(handles.edit_Nrows,'Tooltip',str)
	str = sprintf(['The default value is the number of cols in the grid\n',...
		'However, for reducing border effects you may want to apply\n',...
		'a skirt to the grid. For that, select a value from the side\n',...
		'listbox. Extra points will be padded by mirroiong the south side.']);
	set(handles.edit_Ncols,'Tooltip',str)

	str = sprintf('Good FFT numbers for padding the grid');
	set(handles.listbox_nnx,'Tooltip',str)
	set(handles.listbox_nny,'Tooltip',str)

	if (~isempty(handles.hCallingFig))                    % If we know the handle to the calling fig
		handMir = guidata(handles.hCallingFig);      % get handles of the calling fig
		handles.last_dir = handMir.last_dir;
		handles.home_dir = handMir.home_dir;
		handles.work_dir = handMir.work_dir;
	else
		handles.home_dir = cd;
		handles.last_dir = handles.home_dir;
		handles.work_dir = handles.home_dir;
	end
	handles.path_data = [handles.home_dir filesep 'data' filesep];

	% Import icons
	load([handles.path_data 'mirone_icons.mat'],'Mfopen_ico');
	set(handles.push_BatGrid,'CData',Mfopen_ico)
	set(handles.push_SourceGrid,'CData',Mfopen_ico)

	%------------ Give a Pro look (3D) to the frame boxes  --------
	new_frame3D(hObject, NaN)
	%------------- END Pro look (3D) ------------------------------

	guidata(hObject, handles);
	set(hObject,'Vis','on');
	if (nargout),   varargout{1} = hObject;     end

% -------------------------------------------------------------------------------------------------
function edit_BatGrid_CB(hObject, handles)
	fname = get(hObject,'String');
	if isempty(fname)
		handles.Z_bat = [];		guidata(handles.figure1,handles)
		return
	end
	% Let the push_BatGrid_CB do all the work
	push_BatGrid_CB(handles.push_BatGrid, handles, fname)

% -------------------------------------------------------------------------------------------------
function push_BatGrid_CB(hObject, handles, opt)
	if (nargin == 3),	fname = opt;
	else				opt = [];
	end

	if (isempty(opt))       % Otherwise 'opt' already transmited the file name.
		[FileName,PathName] = put_or_get_file(handles,{ ...
				'*.grd;*.nc', 'Grid files (*.grd,*.nc)';'*.*', 'All Files (*.*)'},'Select GMT grid','get');
		if isequal(FileName,0),		return,		end
		fname = [PathName FileName];
	end

	[handles, X, Y, handles.Z_bat, handles.head_bat] = read_gmt_type_grids(handles,fname);
	if (handles.have_nans)
		errordlg('Bathymetry grid error: Mr. Fourier (FFT) really does not stand the presence of NaNs. Did you know that?', 'Error')
		return
	end

	if (numel(handles.head_bat) > 9)		% May happen when grids are read by "grdread_m"
		handles.head_src = handles.head_bat(1:9);
	end

	% See if Source/Mag grid is already loaded and, if yes, if they are compatible
	if (~isempty(get(handles.edit_SourceGrid,'String')))
		difa_hdrs = abs( diff([handles.head_bat; handles.head_src]) );
		if ( any(difa_hdrs(1:4) > 1e-4) )
			errordlg('Error: Bathymetry & Source grids do not cover the same region','Error'),	return
		elseif ( any(difa_hdrs(8:9) > 1e-6) )
			errordlg('Error: Bathymetry & Source grids do not have the same size.','Error'),	return
		end
	end
	% Try to guess if bat is in meters. If yes convert to km
	if (abs(handles.head_bat(6) - handles.head_bat(5)) > 15)
		grdutils(handles.Z_bat, sprintf('-M%.3f',0.001));
		handles.head_bat(5) = handles.head_bat(5) / 1000;
		handles.head_bat(6) = handles.head_bat(6) / 1000;
	end
	set(handles.edit_BatGrid,'String',fname)
	guidata(handles.figure1,handles)

% -------------------------------------------------------------------------------------------------
function edit_SourceGrid_CB(hObject, handles)
	fname = get(hObject,'String');
	if isempty(fname)
		handles.Z_src = [];		guidata(handles.figure1,handles)
		return
	end
	% Let the push_SourceGrid_CB do all the work
	push_SourceGrid_CB(handles.push_SourceGrid, handles, fname)

% -------------------------------------------------------------------------------------------------
function push_SourceGrid_CB(hObject, handles, opt)
	if (nargin == 3),	fname = opt;
	else				opt = [];
	end

	if (isempty(opt))       % Otherwise 'opt' already transmited the file name.
		[FileName,PathName] = put_or_get_file(handles,{ ...
			'*.grd;*.nc', 'Grid files (*.grd,*.nc)';'*.*', 'All Files (*.*)'},'Select GMT grid','get');
		if isequal(FileName,0),		return,		end
		fname = [PathName FileName];	
	end

	if (~strcmp(fname,'__inMEM'))	% Otherwise we already know those (fished in Mirone handles)
		[handles, handles.X, handles.Y, handles.Z_src, handles.head_src] = read_gmt_type_grids(handles,fname);
		if (numel(handles.head_src) > 9)		% May happen when grids are read by "grdread_m"
			handles.head_src = handles.head_src(1:9);
		end
	end

	if (handles.have_nans)
		errordlg('Magnetic grid error: Mr. Fourier (FFT) really does not stand the presence of NaNs. Did you know that?', 'Error')
		return
	end

	% See if Bat grid is already loaded and, if yes, if both grids are compatible
	if (~isempty(get(handles.edit_BatGrid,'String')))
		difa_hdrs = abs( diff([handles.head_bat; handles.head_src]) );
		if ( any(difa_hdrs(1:4) > 1e-4) )
			errordlg('Error: Bathymetry & Source grids do not cover the same region','Error'),	return
		elseif ( any(difa_hdrs(8:9) > 1e-6) )
			errordlg('Error: Bathymetry & Source grids do not have the same size.','Error'),	return
		end
	end

	[handles.orig_nrows,handles.orig_ncols] = size(handles.Z_src);
	handles.nrows = handles.orig_nrows;             % Make them equal by default
	handles.ncols = handles.orig_ncols;
	set(handles.edit_SourceGrid,'String',fname)
	set(handles.edit_Nrows,'string',sprintf('%d',handles.nrows))
	set(handles.edit_Ncols,'string',sprintf('%d',handles.ncols))

	% Try to guess if grid is in geogs (restricting to 90 degrees spaning is more than enough as a test)
	if (abs(handles.head_src(2)-handles.head_src(1)) < 90 || abs(handles.head_src(4)-handles.head_src(3)) < 90)
		handles.geog = 1;		% We probably have a geog grid
		set(handles.check_geog,'Value',1)
	end

	% The easeast way of not leting the user screw things by selecting a nnx and/or nny lower
	% than nx or ny is to delete the forbiden numbers from the listboxes	
	[w,nlist] = mboard([], handles.ncols, handles.nrows, 0, 0);
	ind = cat(1,nlist{:}) > handles.ncols;
	nlist_t = [{handles.ncols}; nlist(ind)];
	ind = find(cat(1,nlist_t{:}) == w(1));        % Find index of new_nx
	set(handles.listbox_nnx,'String',nlist_t,'Value',ind)
	set(handles.edit_Ncols,'string',sprintf('%d',w(1)))

	ind = cat(1,nlist{:}) > handles.nrows;
	nlist_t = [{handles.orig_nrows}; nlist(ind)];
	ind = find(cat(1,nlist_t{:}) == w(2));        % Find index of new_ny
	set(handles.listbox_nny,'String',nlist_t,'Value',ind)
	set(handles.edit_Nrows,'string',sprintf('%d',w(2)))
	handles.ncols = w(1);
	handles.nrows = w(2);
	
	handles.rlon = (handles.head_src(2) + handles.head_src(1)) / 2;
	handles.rlat = (handles.head_src(4) + handles.head_src(3)) / 2;

	if (handles.geog)
		[sclat,sclon] = scltln(handles.rlat);
		dx = handles.head_src(8) * sclon;
		dy = handles.head_src(9) * sclat;
		handles.scaled_dx = dx;     handles.scaled_dy = dy;
	else
		handles.scaled_dx = handles.head_src(8);		handles.scaled_dy = handles.head_src(9);
	end

	% Compute the wshort & wlong default values (only for the 'inverse' case)
	if (strcmp(handles.what_parker,'inverse'))
		if (handles.geog)
			wshort = max(dx*2, dy*2);
			wlong = max(dx*handles.edit_Ncols, dy*handles.edit_Nrows);
		else
			wshort = max(2*handles.head_src(8),2*handles.head_src(9));
			wlong = max(handles.edit_Ncols*handles.head_src(8),handles.edit_Nrows*handles.head_src(9));
		end
		wlong = max(wlong,150);     % Beter use this as the wlong default 
		set(handles.edit_wshort,'string',sprintf('%.1f',wshort))
		set(handles.edit_wlong, 'string',sprintf('%.0f',wlong))
	end
	guidata(handles.figure1,handles)

% -------------------------------------------------------------------------------------------------
function check_mirror_CB(hObject, handles)
	if (get(hObject,'Value'))
		set([handles.edit_Ncols handles.edit_Nrows handles.listbox_nny handles.listbox_nnx],'Enble','off')
	else
		set([handles.edit_Ncols handles.edit_Nrows handles.listbox_nny handles.listbox_nnx],'Enble','on')
	end

% -------------------------------------------------------------------------------------------------
function edit_date_CB(hObject, handles)
	if (strcmp(handles.what_parker,'redPole'))
		return			% In this mode the box serves only to store a value
	end
	if (get(handles.check_fieldIsRTP,'Value'))
		set(hObject,'String','')
		return
	end
	xx = str2double(get(hObject,'String'));
	if (xx < handles.start_stop_epoch(1) || xx > handles.start_stop_epoch(2))
		errordlg('Date outside the current IGRF model limits','Error');
		handles.date = [];
		set(hObject,'String','');       return
	else
		handles.date = xx;
		elev = str2double(get(handles.edit_zobs,'String'));
		try			% Use a try here because some vars may not have been yet defined
			if (~get(handles.check_CenterDipole,'Value'))
				out = igrf_m(handles.rlon, handles.rlat, elev, xx);
				set(handles.edit_sdec,'String',sprintf('%.1f',out(6)))
				set(handles.edit_sdip,'String',sprintf('%.1f',out(7)))
			end
		end
	end
	guidata(handles.figure1,handles)

% -------------------------------------------------------------------------------------------------
function edit_zobs_CB(hObject, handles)
	zobs = str2double(get(hObject,'String'));
	if (isnan(zobs))
		set(hObject,'String','0'),		return
	else
		handles.zobs = zobs;
	end
	guidata(handles.figure1,handles)

% -------------------------------------------------------------------------------------------------
function edit_thickness_CB(hObject, handles)
	if (strcmp(handles.what_parker,'redPole'))
		return			% In this mode the box serves only to store a value
	end
	xx = str2double(get(hObject,'String'));
	if (xx <= 0 || isnan(xx))
		errordlg('You must be dreaming. What is a layer with zero or negative thickness?','Chico Clever');
		set(hObject,'String',''),	return
	else
		handles.thick = xx;
	end
	guidata(handles.figure1,handles)

% -------------------------------------------------------------------------------------------------
function edit_sdip_CB(hObject, handles)
	if (get(handles.check_fieldIsRTP,'Value'))
		set(hObject,'String','90');    return		% Don't let it be changed if Field is RTP
	end
	if (~strcmp(handles.what_parker,'redPole'))
		if (get(handles.check_CenterDipole,'Value'))     % sdip = 0 for centered dipole
			set(hObject,'String','0');    return
		end
	end
	xx = str2double(get(hObject,'String'));
	if (xx < -90 || xx > 90)
		errordlg('Inlinations are restricted to the [-90;+90] interval.','Error');
		set(hObject,'String','');       return
	else
		handles.sdip = xx;
	end
	guidata(handles.figure1,handles)

% -------------------------------------------------------------------------------------------------
function edit_sdec_CB(hObject, handles)
	if (get(handles.check_fieldIsRTP,'Value'))
		set(hObject,'String','0'),		return        % Don't let it be changed if Field is RTP
	end
	if (~strcmp(handles.what_parker,'redPole'))
		if (get(handles.check_CenterDipole,'Value'))     % sdec = 0 for centered dipole
			set(hObject,'String','0'),	return
		end
	end
	xx = str2double(get(hObject,'String'));
	if (xx < -90 || xx > 90)
		errordlg('Declinations are restricted to the [-90;+90] interval.','Error');
		set(hObject,'String',''),	return
	else
		handles.sdec = xx;
	end
	guidata(handles.figure1,handles)

% -------------------------------------------------------------------------------------------------
function listbox_nnx_CB(hObject, handles)
	contents = get(hObject,'String');
	nnx = str2double(contents{get(hObject,'Value')});
	set(handles.edit_Ncols,'String',sprintf('%d',nnx))
	handles.ncols = nnx;
	guidata(handles.figure1,handles)

% -------------------------------------------------------------------------------------------------
function listbox_nny_CB(hObject, handles)
	contents = get(hObject,'String');
	nny = str2double(contents{get(hObject,'Value')});
	set(handles.edit_Nrows,'String',sprintf('%d',nny))
	handles.nrows = nny;
	guidata(handles.figure1,handles)

% -------------------------------------------------------------------------------------------------
function edit_Ncols_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	if (isempty(get(hObject,'String')))
		try		set(hObject,'String',sprintf('%d',handles.ncols)),	return,		end
	end
	if (xx < handles.cols)
		set(hObject,'String',sprintf('%d',handles.ncols)),	return
	end
	handles.ncols = xx;     guidata(handles.figure1,handles)

% -------------------------------------------------------------------------------------------------
function edit_Nrows_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx))
		try		set(hObject,'String',sprintf('%d',handles.nrows)),	return,		end
	end
	if (xx < handles.nrows)
		set(hObject,'String',sprintf('%d',handles.nrows)),	return
	end
	handles.nrows = xx;     guidata(handles.figure1,handles)

% -------------------------------------------------------------------------------------------------
function check_fieldIsRTP_CB(hObject, handles)
	if (get(hObject,'Value'))
		set(handles.edit_sdip,'String','90');   set(handles.edit_sdec,'String','0')
		set(handles.edit_date,'String','');     set(handles.check_CenterDipole,'Value',0)
	else
		set(handles.edit_sdip,'String','0');   set(handles.edit_sdec,'String','0')
		set(handles.check_CenterDipole,'Value',1)
	end

% -------------------------------------------------------------------------------------------------
function check_CenterDipole_CB(hObject, handles)
	if (get(handles.check_fieldIsRTP,'Value'))
		set(hObject,'Value',0)
	end

% -------------------------------------------------------------------------------------------------
function [G, Hdr] = push_compute_CB(hObject, handles)
% When used with output args, no figure is created.
%
% Before asking the apropriate function to do the work we have to ... TEST
% Source grid, Nrows & Ncols are common to all options. So test them first
	if (isempty(handles.Z_src))
		errordlg('You didn''t give me a Source grid (Field or Magnetization). What do you want me to do?','Chico Clever')
		return
	end
	if (isempty(get(handles.edit_Nrows,'String')) || isempty(get(handles.edit_Ncols,'String')))
		errordlg('One or both of grid size dimensions are empty. What have you done?','Error')
		return
	end

	% Now those that are common to direct/inverse cases
	if (strcmp(handles.what_parker,'direct') || strcmp(handles.what_parker,'inverse'))
		if (isempty(handles.Z_bat))
			errordlg('Must give me a grid with the bathymetry','Error');    return
		end
		date = str2double(get(handles.edit_date,'String'));
		if (isnan(date) && ~get(handles.check_fieldIsRTP,'Value'))
			errordlg('I need to know the year of the survey (see Date box)','Error');   return
		end
		thick = str2double(get(handles.edit_thickness,'String'));
		if (isnan(thick))
			errordlg('I need to know the thickness of the magnetic layer (see Thickness box)','Error');   return
		end
		zobs = str2double(get(handles.edit_zobs,'String'));
		if (isnan(zobs))
			errordlg('I need to know the level of the observation of the survey (see Level box)','Error');   return
		end
		dx = handles.scaled_dx;     dy = handles.scaled_dy;
	end
	% Get Mag/Field Dec & Dip
	sdec = str2double(get(handles.edit_sdec,'String'));
	sdip = str2double(get(handles.edit_sdip,'String'));
	new_nx = str2double(get(handles.edit_Ncols,'String'));
	new_ny = str2double(get(handles.edit_Nrows,'String'));

	switch (handles.what_parker)
		case 'direct'
			if (get(handles.check_fieldIsRTP,'Value'))			% Case of a mag computed from an RTP Field
				handles.rlat = 90;      handles.rlon = 0;		% Tricky values to inform syn3d about the RTP fact
				sdip = 90;              sdec = 0;       date = 2000;
			end
			if (handles.orig_ncols < handles.ncols || handles.orig_nrows < handles.nrows)      % Padding was asked for
				if (get(handles.check_mirror,'Value'))   % Do mirror
					h = mboard(handles.Z_bat,handles.orig_ncols,handles.orig_nrows);
					f = mboard(handles.Z_src,handles.orig_ncols,handles.orig_nrows);
					f3d = syn3d(double(f),double(h),handles.rlat,handles.rlon,date,zobs,thick,0,dx,dy,sdip,sdec);
					f3d = f3d(1:handles.orig_nrows,1:handles.orig_ncols);   % Remove the mirror
				else
					[h,band] = mboard(handles.Z_bat,handles.orig_ncols,handles.orig_nrows,new_nx,new_ny);
					f = mboard(handles.Z_src,handles.orig_ncols,handles.orig_nrows,new_nx,new_ny);
					f3d = syn3d(double(f),double(h),handles.rlat,handles.rlon,date,zobs,thick,0,dx,dy,sdip,sdec);
					m1 = band(1)+1;     m2 = m1 + handles.orig_nrows - 1;
					n1 = band(3)+1;     n2 = n1 + handles.orig_ncols - 1;
					f3d = f3d(m1:m2,n1:n2);         % Remove the padding skirt
				end
			else
				f3d = syn3d(double(handles.Z_src),double(handles.Z_bat),handles.rlat,handles.rlon, ...
					date,zobs,thick,0,dx,dy,sdip,sdec);
			end
			z_min = min(f3d(:));    z_max = max(f3d(:));
			tmp.head = [handles.head_src(1) handles.head_src(2) handles.head_src(3) handles.head_src(4) ...
					z_min z_max 0 handles.head_src(8) handles.head_src(9)];
			tmp.X = handles.X;      tmp.Y = handles.Y;      tmp.name = 'Magnetic Anomaly (nT)';
			if (nargout == 0)				% Show result right away
				f3d = single(f3d);	mirone(f3d,tmp);
			else							% Send back the result to the caller
				G = f3d;	Hdr = tmp;
			end
		case 'inverse'
			if (get(handles.check_fieldIsRTP,'Value'))			% Case of a already RTP Field
				handles.rlat = 90;      handles.rlon = 0;		% Tricky values to inform inv3d about the RTP fact
				sdip = 90;              sdec = 0;       date = 2000;
			end
			ws = str2double(get(handles.edit_wshort,'String'));
			wl = str2double(get(handles.edit_wlong,'String'));
			if (handles.orig_ncols < handles.ncols || handles.orig_nrows < handles.nrows)      % Padding was asked for
				if (get(handles.check_mirror,'Value'))   % Do mirror
					h = mboard(handles.Z_bat,handles.orig_ncols,handles.orig_nrows);
					f = mboard(handles.Z_src,handles.orig_ncols,handles.orig_nrows);
					m3d = inv3d(double(f),double(h),wl,ws,handles.rlat,handles.rlon,date,zobs,thick,0,dx,dy,sdec,sdip);
					m3d = m3d(1:handles.orig_nrows,1:handles.orig_ncols);   % Remove the mirror
				else
					[h,band] = mboard(handles.Z_bat,handles.orig_ncols,handles.orig_nrows,new_nx,new_ny);
					f = mboard(handles.Z_src,handles.orig_ncols,handles.orig_nrows,new_nx,new_ny);
					m3d = inv3d(double(f),double(h),wl,ws,handles.rlat,handles.rlon,date,zobs,thick,0,dx,dy,sdec,sdip);
					m1 = band(1)+1;     m2 = m1 + handles.orig_nrows - 1;
					n1 = band(3)+1;     n2 = n1 + handles.orig_ncols - 1;
					m3d = m3d(m1:m2,n1:n2);         % Remove the padding skirt
				end
			else
				m3d = inv3d(double(handles.Z_src),double(handles.Z_bat),wl,ws,handles.rlat,handles.rlon, ...
					date,zobs,thick,0,dx,dy,sdec,sdip);
			end
			z_min = min(m3d(:));        z_max = max(m3d(:));
			tmp.head = [handles.head_src(1) handles.head_src(2) handles.head_src(3) handles.head_src(4) ...
					z_min z_max 0 handles.head_src(8) handles.head_src(9)];
			tmp.X = handles.X;      tmp.Y = handles.Y;      tmp.name = 'Magnetization (A/m^2)';

			if (get(handles.check_showDirect, 'Val') || get(handles.check_showResidue, 'Val'))
				handles.what_parker = 'direct';
				handles.Z_src = m3d;		% Since we want the forward solution we must now feed it with the mag
				[G, Hdr] = push_compute_CB([], handles);
				if (get(handles.check_showDirect, 'Val'))
					Hdr.name = 'Forward Anomaly (nT)';
					G = single(G);		mirone(G, Hdr)
				end
				if (get(handles.check_showResidue, 'Val'))
					G = handles.Z_src - G;
					Hdr.head(5:6) = [min(G(:)) max(G(:))];
					Hdr.name = 'Residual Anomaly (nT)';
					G = single(G);		mirone(G, Hdr)
				end
			end
			if (nargout == 0)				% Show result right away
				m3d = single(m3d);		mirone(m3d,tmp);
			else							% Send back the result to the caller
				G = m3d;	Hdr = tmp;
			end
		case 'redPole'
			incl_fld = str2double(get(handles.edit_zobs,'String'));
			decl_fld = str2double(get(handles.edit_thickness,'String'));
			incl_mag = str2double(get(handles.edit_sdip,'String'));
			decl_mag = str2double(get(handles.edit_sdec,'String'));
			if (isnan(incl_fld) || isnan(decl_fld))
				errordlg('You need to give me valid magnetic field Inclination and Declination.','Error');  return;
			end
			if (isnan(incl_mag) || isnan(decl_mag))
				errordlg('You need to give me valid magnetization Inclination and Declination.','Error');  return;
			end
			if (handles.orig_ncols < handles.ncols || handles.orig_nrows < handles.nrows)      % Padding was asked for
				if (get(handles.check_mirror,'Value'))   % Do mirror
					f = mboard(handles.Z_src,handles.orig_ncols,handles.orig_nrows);
					f = rtp3d(double(f),incl_fld,decl_fld,incl_mag,decl_mag);
					f = f(1:handles.orig_nrows,1:handles.orig_ncols);   % Remove the mirror
				else
					[f,band] = mboard(handles.Z_src,handles.orig_ncols,handles.orig_nrows,new_nx,new_ny);
					f = rtp3d(double(f),incl_fld,decl_fld,incl_mag,decl_mag);
					m1 = band(1)+1;     m2 = m1 + handles.orig_nrows - 1;
					n1 = band(3)+1;     n2 = n1 + handles.orig_ncols - 1;
					f = f(m1:m2,n1:n2);         % Remove the padding skirt
				end
			else
				f = rtp3d(double(handles.Z_src),incl_fld,decl_fld,incl_mag,decl_mag);
			end
			z_min = min(f(:));    z_max = max(f(:));
			tmp.head = [handles.head_src(1:4) z_min z_max 0 handles.head_src(8:9)];
			tmp.X = handles.X;      tmp.Y = handles.Y;      tmp.name = 'Reduce to the Pole anomaly (nT)';
			if (nargout == 0)				% Show result right away
				f = single(f);		mirone(f,tmp);
			else
				G = f;		Hdr = tmp;
			end
	end

% -------------------------------------------------------------------------------------------------
function [sclat,sclon] = scltln(orlat)
% Routine to determine lat-lon scales, km/deg, for ellipsoids
% of revolution,  using equations of:
%       Snyder, J.P., 1987, Map Projections -- A Working Manual,
%       USGS Professional Paper 1395, Washington DC, 383p. cf. pp 24-25.
%
% Currently, this is hard-wired for the WGS-84 ellipsoid.
%
% The basic equations are:
% 	sclat = a * (1-e*e)    /  (1 - e*e * sin(orlat)*sin(orlat))**(1.5)
%	sclon = a * cos(orlat) /  (1 - e*e * sin(orlat)*sin(orlat))**(0.5)
%
% where:    a  is the equatorial radius
%           b  is the polar radius
%           e  is the eccentricity
%           f  is the flattening
% Also:
%	e*e = 1. - b*b/a*a
%	f   = 1. - b/a
%
% Dan Scheirer, 21 May 1991

	% These constants belong to the: WGS, 1984 ellipsoid (gmt_defaults.h)
	a = 6378.137;   b = 6356.7521;

	% Now, do the calculations...
	e2 = 1 - (b*b)/(a*a);
	sinlat = sin(orlat*pi/180);
	denom  = sqrt(1 - e2 * sinlat * sinlat);
	sclat = (pi/180) * a * (1 - e2)  / denom / denom / denom;
	sclon = (pi/180) * a * cos(orlat*pi/180) / denom;


% --- Creates and returns a handle to the GUI figure. 
function parker_stuff_LayoutFcn(h1)

set(h1, 'PaperUnits',get(0,'defaultfigurePaperUnits'),...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'MenuBar','none',...
'Name','parker_stuff',...
'NumberTitle','off',...
'Position',[520 498 400 302],...
'Resize','off',...
'Tag','figure1');

uicontrol('Parent',h1,'Position',[10 11 301 71],'String',{''},'Style','frame');
uicontrol('Parent',h1,'Position',[10 101 381 81],'String',{''},'Style','frame');
uicontrol('Parent',h1,'Position',[10 201 381 91],'String',{''},'Style','frame');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@parker_stuff_uiCB,...
'HorizontalAlignment','left',...
'Position',[50 232 311 21],...
'Style','edit',...
'Tooltip','Enter bathymetry grid name (km +ve up)',...
'Tag','edit_BatGrid');

uicontrol('Parent',h1,...
'Call',@parker_stuff_uiCB,...
'Position',[361 231 23 23],...
'Tag','push_BatGrid');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@parker_stuff_uiCB,...
'HorizontalAlignment','left',...
'Position',[50 262 311 21],...
'Style','edit',...
'Tooltip','Enter magnetic field grid name',...
'Tag','edit_SourceGrid');

uicontrol('Parent',h1,'HorizontalAlignment','left','Position',[17 236 31 15],...
'String','Bat','Style','text','Tag','text1');

uicontrol('Parent',h1,'HorizontalAlignment','left','Position',[18 265 31 15],...
'String','Field','Style','text','Tag','text_FieldMag');

uicontrol('Parent',h1,...
'Call',@parker_stuff_uiCB,...
'Position',[361 261 23 23],...
'Tag','push_SourceGrid');

uicontrol('Parent',h1,...
'Position',[50 210 140 15],...
'String','Geographic coords?',...
'Style','checkbox',...
'Tooltip','Are the grids in geographical coordinates?',...
'Tag','check_geog');

uicontrol('Parent',h1,...
'Call',@parker_stuff_uiCB,...
'Position',[225 210 140 15],...
'String','Field is already RTP',...
'Style','checkbox',...
'Tooltip','Check this box if the anomalous field is already Reduced To the Pole',...
'Tag','check_fieldIsRTP');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@parker_stuff_uiCB,...
'Position',[263 53 41 21],...
'Style','edit',...
'Tooltip','Decimal year for IGRF field calculation',...
'Tag','edit_date');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@parker_stuff_uiCB,...
'Position',[162 53 41 21],...
'String','0.5',...
'Style','edit',...
'Tooltip','Thickness of magnetic source layer (km)',...
'Tag','edit_thickness');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@parker_stuff_uiCB,...
'Position',[55 53 41 21],...
'String','0',...
'Style','edit',...
'Tooltip','Observation level above sealevel (km +ve up)',...
'Tag','edit_zobs');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@parker_stuff_uiCB,...
'Position',[55 18 41 21],...
'Style','edit',...
'Tooltip','inclination of magnetization',...
'Tag','edit_sdip');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@parker_stuff_uiCB,...
'Position',[161 18 41 21],...
'Style','edit',...
'Tooltip','declination of magnetization',...
'Tag','edit_sdec');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Position',[335 150 41 21],...
'Style','edit',...
'Tooltip','filter short wavelength cutoff (km)',...
'Tag','edit_wshort');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Position',[335 110 41 21],...
'Style','edit',...
'Tooltip','filter long wavelength cutoff (km)',...
'Tag','edit_wlong');

uicontrol('Parent',h1,'HorizontalAlignment','right','Position',[229 56 31 15],...
'String','Date','Style','text');

uicontrol('Parent',h1,'Position',[110 56 49 15],'String','Thickness','Style','text','Tag','text_Thickness');

uicontrol('Parent',h1,'HorizontalAlignment','right','Position',[12 56 40 15],...
'String','Level','Style','text','Tag','text_Level');

uicontrol('Parent',h1,'HorizontalAlignment','right','Position',[107 21 51 15],...
'String','Mag dec','Style','text');

uicontrol('Parent',h1,'HorizontalAlignment','right','Position',[12 22 40 15],...
'String','Mag dip','Style','text');

uicontrol('Parent',h1,'Position',[291 153 41 15],'String','Wshort','Style','text','Tag','text_wshort');

uicontrol('Parent',h1,'HorizontalAlignment','right','Position',[295 114 35 15],...
'String','Wlong','Style','text','Tag','text_wlong');

uicontrol('Parent',h1,...
'Call',@parker_stuff_uiCB,...
'Position',[19 139 60 15],...
'String','Mirror',...
'Style','checkbox',...
'Tooltip','Check this to Mirror the grid before FFT',...
'Tag','check_mirror');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@parker_stuff_uiCB,...
'Position',[126 102 51 78],'Style','listbox','Value',1,'Tag','listbox_nny');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@parker_stuff_uiCB,...
'Position',[232 102 51 78],...
'Style','listbox',...
'Value',1,...
'Tag','listbox_nnx');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@parker_stuff_uiCB,...
'Position',[86 135 41 21],...
'Style','edit',...
'Tooltip','Number of grid rows',...
'Tag','edit_Nrows');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@parker_stuff_uiCB,...
'Position',[192 135 41 21],...
'Style','edit',...
'Tooltip','Number of grid columns',...
'Tag','edit_Ncols');

uicontrol('Parent',h1,'Position',[86 159 39 15],'String','# Rows','Style','text','Tag','text11');
uicontrol('Parent',h1,'Position',[192 159 39 15],'String','# Cols','Style','text','Tag','text12');

uicontrol('Parent',h1,...
'Call',@parker_stuff_uiCB,...
'Position',[227 19 90 15],...
'String','Geocentric?',...
'Style','checkbox',...
'Tooltip','Check this to assume geocentric dipole hypothesis',...
'Tag','check_CenterDipole');

uicontrol('Parent',h1,...
'Position',[320 70 90 15],...
'String','Show res',...
'Style','checkbox',...
'Tooltip','Show also the residue grid (original anomaly - forward solution)',...
'Tag','check_showResidue');

uicontrol('Parent',h1,...
'Position',[320 50 90 15],...
'String','Show direct',...
'Style','checkbox',...
'Tooltip','Show also the forward solution grid',...
'Tag','check_showDirect');

uicontrol('Parent',h1,...
'Call',@parker_stuff_uiCB,...
'FontWeight','bold',...
'Position',[326 11 65 21],...
'String','Compute',...
'Tag','push_compute');

function parker_stuff_uiCB(hObject, eventdata)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));

function varargout = fft_stuff(varargin)
% Helper window to do FFT operations 

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

	hObject = figure('Tag','figure1','Visible','off');
	fft_stuff_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hObject,'right')

	handles.hMirFig = [];			% Handle to the calling figure
	handles.geog = 0;				% Set this as default
	handles.Z1 = [];
	handles.Z2 = [];
	grid_in_continue = 0;

	if (~isempty(varargin))         % When called from a Mirone window
		if (ishandle(varargin{1}))
			handles.hMirFig = varargin{1};
			handles.Z1 = varargin{2};
			handles.head_Z1 = varargin{3};
			handles.geog = varargin{4};
			mode = varargin{5};
		elseif (isa(varargin{1}, 'struct'))			% Currently only used by ECRAN -> Spector & Grant
			handles.Z1 = varargin{1}.Z;
			handles.head_Z1 = varargin{1}.head;
			handles.geog = varargin{1}.geog;
			mode = varargin{1}.mode;
			handles.bandpass = varargin{1}.bandpass;	% FREQUENCY in 1/m
			if (isfield(varargin{1},'hMirFig'))		handles.hMirFig = varargin{1}.hMirFig;		end % MANDATORY
		end
		if (strcmp(mode,'Allopts')),	grid_in_continue = 1;	end     % "Slow mode" show all options in this figure
		[handles.orig_nrows,handles.orig_ncols] = size(handles.Z1);
		[w,nlist] = mboard([],handles.orig_ncols,handles.orig_nrows,0,0);
		handles.new_nx = w(1);
		handles.new_ny = w(2);
		rlat = (handles.head_Z1(4) + handles.head_Z1(3)) / 2;
		if (handles.geog)
			[sclat,sclon] = scltln(rlat);
			dx = handles.head_Z1(8) * sclon;
			dy = handles.head_Z1(9) * sclat;
			handles.scaled_dx = dx;     handles.scaled_dy = dy;
			handles.is_meters = 0;      handles.is_km = 0;
		else					% Guess if grid units are meters or km
			dx = handles.head_Z1(2) - handles.head_Z1(1);
			dy = handles.head_Z1(4) - handles.head_Z1(3);
			len = sqrt(dx.*dx + dy.*dy);         % Distance in user unites
			if (len > 1e4 || handles.head_Z1(8) > 10)      % If grid's diagonal > 1e4 consider we have meters
				handles.is_meters = 1;     handles.is_km = 0;
			else				% km
				handles.is_meters = 0;     handles.is_km = 1;
			end
			% It's actually very hard to guess if m or km so I'll set it to meters always
			% but leave the above lines to remind the case which can only be adressed well
			% if the grid has header info on X,Y units.
			handles.is_meters = 1;     handles.is_km = 0;
			handles.scaled_dx = handles.head_Z1(8);
			handles.scaled_dy = handles.head_Z1(9);
			if (handles.is_km)
				handles.scaled_dx = handles.scaled_dx * 1000;
				handles.scaled_dy = handles.scaled_dy * 1000;
			end
		end
		if (~grid_in_continue)       % Called in the "quick mode" 
			if (any( strcmp( mode,{'Power' 'Autocorr' 'Amplitude' 'lpass' 'hpass' 'bpass'} ) ))
				sectrumFun(handles, handles.Z1, handles.head_Z1, mode)
			end
			delete(hObject),		return
		end
	end

	if (~isempty(handles.hMirFig))
		handMir = guidata(handles.hMirFig);
		handles.home_dir = handMir.home_dir;
		handles.work_dir = handMir.work_dir;
		handles.last_dir = handMir.last_dir;
		handles.path_data = handMir.path_data;
	else
		handles.home_dir = cd;
		handles.work_dir = cd;		handles.last_dir = cd;	% To not compromize put_or_get_file
		handles.path_data = [cd filesep 'data' filesep];	% Wil fail if called from outside Mir home
	end

	if (grid_in_continue)       % Grid recieved in argument. Fill the listboxes
		% The easeast way of not leting the user screw things by selecting a nnx and/or nny
		% lower than nx or ny is to delete the forbiden numbers from the listboxes
		ind = cat(1,nlist{:}) > handles.orig_ncols;
		nlist_t = [{handles.orig_ncols}; nlist(ind)];
		ind = find(cat(1,nlist_t{:}) == handles.new_nx);        % Find index of new_nx
		set(handles.listbox_nnx,'String',nlist_t,'Value',ind)
		set(handles.edit_Ncols,'string',sprintf('%d',handles.new_nx))
		
		ind = cat(1,nlist{:}) > handles.orig_nrows;
		nlist_t = [{handles.orig_nrows}; nlist(ind)];
		ind = find(cat(1,nlist_t{:}) == handles.new_ny);        % Find index of new_ny
		set(handles.listbox_nny,'String',nlist_t,'Value',ind)
		set(handles.edit_Nrows,'string',sprintf('%d',handles.new_ny))
		
		handles.X = linspace(handles.head_Z1(1),handles.head_Z1(2),handles.orig_ncols);
		handles.Y = linspace(handles.head_Z1(3),handles.head_Z1(4),handles.orig_nrows);
		set(handles.edit_Grid1,'String',' In Memory array')
	end

	% Import icons
	load([handles.path_data 'mirone_icons.mat'],'Mfopen_ico');
	set(handles.push_Grid2,'CData',Mfopen_ico)
	set(handles.push_Grid1,'CData',Mfopen_ico)
	clear Mfopen_ico;

	% Set upt some useful tooltips
	str = sprintf(['The default value is the number of rows in the grid\n',...
		'However, for reducing border effects you may want to apply\n',...
		'a skirt to the grid. For that, select a value from the side\n',...
		'listbox. Extra points are obtained by data extention and tapered to zero.']);
	set(handles.edit_Nrows,'TooltipString',str)
	str = sprintf(['The default value is the number of cols in the grid\n',...
		'However, for reducing border effects you may want to apply\n',...
		'a skirt to the grid. For that, select a value from the side\n',...
		'listbox. Extra points are obtained by data extention and tapered to zero.']);
	set(handles.edit_Ncols,'TooltipString',str)

	str = sprintf('Good FFT numbers for padding the grid');
	set(handles.listbox_nnx,'TooltipString',str)
	set(handles.listbox_nny,'TooltipString',str)

	%------------ Give a Pro look (3D) to the frame boxes  -------------------------------
	new_frame3D(hObject, NaN)
	%------------- END Pro look (3D) -----------------------------------------------------

	% Add this figure handle to the carraças list
	if (~isempty(handles.hMirFig))
		plugedWin = getappdata(handles.hMirFig,'dependentFigs');
		plugedWin = [plugedWin hObject];
		setappdata(handles.hMirFig,'dependentFigs',plugedWin);
	end

	guidata(hObject, handles);
	set(hObject,'Visible','on');
	if (nargout),	varargout{1} = hObject;		end

% -------------------------------------------------------------------------------------------------
function edit_Grid1_CB(hObject, handles)
	fname = get(hObject,'String');
	if isempty(fname),   handles.Z1 = [];    return;     end
	% Let the push_Grid1_CB do all the work
	push_Grid1_CB(handles.push_Grid1, handles, fname)

% -------------------------------------------------------------------------------------------------
function push_Grid1_CB(hObject, handles, opt)
	if (nargin == 3)	fname = opt;
	else				opt = [];
	end

	if (isempty(opt))       % Otherwise 'opt' already transmited the file name.
		[FileName,PathName] = put_or_get_file(handles, ...
			{'*.grd;*.nc', 'Grid files (*.grd,*.nc)';'*.*', 'All Files (*.*)'},'Select GMT grid','get');
		if isequal(FileName,0),		return,		end
		fname = [PathName FileName];
	end

	[handles, handles.X, handles.Y, handles.Z1, handles.head_Z1] = read_gmt_type_grids(handles,fname);

	if (grdutils(handles.Z1,'-N'))
		errordlg('This grid has NaNs. That is not allowed in FFTs','Error');    return;
	end

	% See if Grid2 is already loaded and, if yes, if both grids are compatible
	if (~isempty(get(handles.edit_Grid2,'String')))
		difa_hdrs = abs( diff([handles.head_Z1; handles.head_Z2]) );
		if ( any(difa_hdrs(1:4) > 1e-4) )
			errordlg('Error: Grid1 & Grid2 do not cover the same region','Error'),	return
		elseif ( any(difa_hdrs(8:9) > 1e-6) )
			errordlg('Error: Grid1 & Grid2 do not have the same size.','Error'),	return
		end
	end
	[handles.orig_nrows,handles.orig_ncols] = size(handles.Z1);
	set(handles.edit_Grid1,'String',fname)

	[ns,nlist] = mboard([],handles.orig_ncols,handles.orig_nrows,0,0);
	handles.new_nx = ns(1);
	handles.new_ny = ns(2);

	% The easeast way of not leting the user screw things by selecting a nnx and/or nny lower
	% than nx or ny is to delete the forbiden numbers from the listboxes
	ind = cat(1,nlist{:}) > handles.orig_ncols;
	nlist_t = [{handles.orig_ncols}; nlist(ind)];
	ind = find(cat(1,nlist_t{:}) == handles.new_nx);        % Find index of new_nx
	set(handles.listbox_nnx,'String',nlist_t,'Value',ind)
	set(handles.edit_Ncols,'string',sprintf('%d',handles.new_nx))

	ind = cat(1,nlist{:}) > handles.orig_nrows;
	nlist_t = [{handles.orig_nrows}; nlist(ind)];
	ind = find(cat(1,nlist_t{:}) == handles.new_ny);        % Find index of new_ny
	set(handles.listbox_nny,'String',nlist_t,'Value',ind)
	set(handles.edit_Nrows,'string',sprintf('%d',handles.new_ny))

	% Try to guess if grid is in geogs
	if (abs(handles.head_Z1(2)-handles.head_Z1(1)) < 180 || abs(handles.head_Z1(4)-handles.head_Z1(3)) < 170)
		handles.geog = 1;	% We probably have a geog grid
		handles.is_meters = 0;     handles.is_km = 0;
		set(handles.popup_GridCoords,'Value',1)
	else
		dx = handles.head_Z1(2) - handles.head_Z1(1);
		dy = handles.head_Z1(4) - handles.head_Z1(3);
		len = sqrt(dx.*dx + dy.*dy);         % Distance in user unites
		if (len > 1e5)		% If grid's diagonal > 1e5 consider we have meters
			handles.is_meters = 1;     handles.is_km = 0;   handles.geog = 0;
			set(handles.popup_GridCoords,'Value',2)
		else				% km
			handles.is_meters = 0;     handles.is_km = 1;   handles.geog = 0;
			set(handles.popup_GridCoords,'Value',3)
		end
	end

	rlat = (handles.head_Z1(4) + handles.head_Z1(3)) / 2;

	if (handles.geog)
		[sclat,sclon] = scltln(rlat);
		dx = handles.head_Z1(8) * sclon;
		dy = handles.head_Z1(9) * sclat;
		handles.scaled_dx = dx;     handles.scaled_dy = dy;
	else
		handles.scaled_dx = handles.head_Z1(8);
		handles.scaled_dy = handles.head_Z1(9);
		if (handles.is_km)
			handles.scaled_dx = handles.scaled_dx * 1000;
			handles.scaled_dy = handles.scaled_dy * 1000;
		end
	end
	guidata(hObject,handles)

% -------------------------------------------------------------------------------------------------
function edit_Grid2_CB(hObject, handles)
	fname = get(hObject,'String');
	if isempty(fname),    handles.Z2 = [];    return;     end
	% Let the push_Grid2_CB do all the work
	push_Grid2_CB(handles.push_Grid2, handles, fname)

% -------------------------------------------------------------------------------------------------
function push_Grid2_CB(hObject, handles, opt)
	if (nargin == 3)	fname = opt;
	else				opt = [];
	end

	if (isempty(opt))       % Otherwise 'opt' already transmited the file name.
		[FileName,PathName] = put_or_get_file(handles,{...
			'*.grd;*.nc', 'Grid files (*.grd,*.nc)';'*.*', 'All Files (*.*)'},'Select GMT grid','get');
		if isequal(FileName,0),		return,		end
		fname = [PathName FileName];
	end

	[handles, X, Y, handles.Z2, handles.head_Z2] = read_gmt_type_grids(handles,fname);

	if (grdutils(handles.Z2,'-N'))
		errordlg('This grid has NaNs. That is not allowed in FFTs','Error'),	return;
	end
	% See if Grid1 grid is already loaded and, if yes, if they are compatible
	if (~isempty(get(handles.edit_Grid1,'String')))
		difa_hdrs = abs( diff([handles.head_Z1; handles.head_Z2]) );
		if ( any(difa_hdrs(1:4) > 1e-4) )
			errordlg('Error: Grid1 & Grid2 do not cover the same region','Error'),	return
		elseif ( any(difa_hdrs(8:9) > 1e-6) )
			errordlg('Error: Grid1 & Grid2 do not have the same size.','Error'),	return
		end
	end
	set(handles.edit_Grid2,'String',fname)
	guidata(handles.figure1, handles)

% -------------------------------------------------------------------------------------------------
function listbox_nnx_CB(hObject, handles)
	contents = get(hObject,'String');
	nnx = str2double(contents{get(hObject,'Value')});
	set(handles.edit_Ncols,'String',sprintf('%d',nnx))
	handles.new_nx = nnx;     guidata(hObject,handles)

% -------------------------------------------------------------------------------------------------
function listbox_nny_CB(hObject, handles)
	contents = get(hObject,'String');
	nny = str2double(contents{get(hObject,'Value')});
	set(handles.edit_Nrows,'String',sprintf('%d',nny))
	handles.new_ny = nny;     guidata(hObject,handles)

% -------------------------------------------------------------------------------------------------
function edit_Ncols_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	if (isempty(get(hObject,'String')))
		try    set(hObject,'String',sprintf('%d',handles.ncols)),	return,		end
	end
	if (xx < handles.cols)   
		set(hObject,'String',sprintf('%d',handles.ncols));    return;
	end
	handles.ncols = xx;     guidata(hObject,handles)

% -------------------------------------------------------------------------------------------------
function edit_Nrows_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx))
		try    set(hObject,'String',sprintf('%d',handles.nrows)),	return,		end
	end
	if (xx < handles.nrows)
		set(hObject,'String',sprintf('%d',handles.nrows)),	return
	end
	handles.nrows = xx;     guidata(hObject,handles)

% -------------------------------------------------------------------------------------------------
function edit_dirDerivative_CB(hObject, handles)
	azim = str2double(get(hObject,'String'));
	if (isnan(azim)),   set(hObject,'String','0');  end

% -------------------------------------------------------------------------------------------------
function edit_derivative_CB(hObject, handles)
	n_der = str2double(get(hObject,'String'));
	if (isnan(n_der)),		set(hObject,'String','1'),		end

% -------------------------------------------------------------------------------------------------
function push_powerSpectrum_CB(hObject, handles)
	if (isempty(handles.Z1))   
		errordlg('No grid loaded. Rebooting ...','Error');  return
	end
	sectrumFun(handles, handles.Z1, handles.head_Z1, 'Power')

% --------------------------------------------------------------------
function push_crossSpectra_CB(hObject, handles)
	if (isempty(handles.Z1) || isempty(handles.Z2))
		errordlg('"Cross" means one grid against the other. You got the idea?','Error');
		return
	end
	sectrumFun(handles, handles.Z1, handles.head_Z1, 'CrossPower', handles.Z2)

% -------------------------------------------------------------------------------------------------
function push_autoCorr_CB(hObject, handles)
	sectrumFun(handles, handles.Z1, handles.head_Z1, 'Autocorr')

% -------------------------------------------------------------------------------------------------
function push_crossCorr_CB(hObject, handles)
	if (isempty(handles.Z1) || isempty(handles.Z2))
		errordlg('"Cross" means one grid against the other. You got the idea?','Error');
		return
	end
	sectrumFun(handles, handles.Z1, handles.head_Z1, 'CrossCorrel', handles.Z2)

% -------------------------------------------------------------------------------------------------
function push_radialPowerAverage_CB(hObject, handles)
% This is modeled on the 1-D case, using the following ideas:
% In 1-D, we ensemble average over samples of length L = n * dt.
% This gives n/2 power spectral estimates, at frequencies i/L,
% where i = 1, n/2.  If we have a total data set of ndata,
% we can make nest=ndata/n independent estimates of the power.
% Our standard error is then 1/sqrt(nest). In making 1-D estimates
% from 2-D data, we may choose n and L from nx2 or ny2 and delta_kx, delta_ky 
% as appropriate.  In this routine, we are giving the sum over all frequencies
% in the other dimension; that is, an approximation of the integral.
%   ORIGINAL TEXT FROM WALTER SMITH IN GRDFFT

	if (isempty(handles.Z1))	errordlg('No grid loaded yet.','Error'),	return,		end

	nx2 = handles.new_nx;		ny2 = handles.new_ny;
	delta_kx = 2*pi / (nx2 * handles.scaled_dx);
	delta_ky = 2*pi / (ny2 * handles.scaled_dy);
	if (delta_kx < delta_ky)
		delta_k = delta_kx;		nk = fix(nx2/2);
	else
		delta_k = delta_ky;		nk = fix(ny2/2);
	end
	r_delta_k = 1 / delta_k;
	set(handles.figure1,'pointer','watch')
	if (~isa(handles.Z1, 'double'))		handles.Z1 = double(handles.Z1);	end
	if (get(handles.checkbox_leaveTrend,'Value'))       % Remove trend
		handles.Z1 = grdtrend_m(handles.Z1, handles.head_Z1, '-D', '-N3');
	end
	[Z,band,modk] = wavenumber_and_mboard(handles);
	ifreq = round(modk * r_delta_k) + 1;     clear modk;
	Z = fft2(Z);    Z(1,1) = 0;
	Z = Z.* conj(Z);

	ifreq(ifreq < 1) = 1;           % Might happen when doing r spectrum
	nk1 = nk + 1;
	power = zeros(nk1,1);
	n_used = 0;
	for (m=1:ny2)
		for(n=1:nx2)
			if (ifreq(m,n) > nk1)		continue,		end		% Might happen when doing r spectrum
			power(ifreq(m,n)) = power(ifreq(m,n)) + Z(m,n);
			n_used = n_used + 1;
		end
	end
	power(1) = [];						% Remove DC component

	delta_k = delta_k / (2*pi);			% Write out frequency, not wavenumber
	powfactor = 1 / (nx2*ny2)^2;
	power = power * powfactor;
	freq  = (1:nk) * delta_k;

	if (handles.geog || handles.is_km)
		freq = freq * 1000;
		x_label = 'Frequency (1/km)';	 % Report frequency in 1/km
	else
		x_label = 'Frequency (1/m)';
	end

	if (~isempty(handles.hMirFig))			arg1 = handles.hMirFig;
	else									arg1 = 'reuse';
	end

	set(handles.figure1,'pointer','arrow')
	ecran(arg1, freq, power,'','Radial average power spectrum',x_label,'Log(Power)','','semilogy')

% -------------------------------------------------------------------------------
function calSave(obj,eventdata,h_fig)
% Save data in file
	ud = get(h_fig,'UserData');          % Get userdata
	handles = guidata(ud.h_mir_fig);
	[FileName,PathName] = put_or_get_file(handles, ...
		{'*.dat;*.DAT', 'radial spectrum file (*.dat,*.DAT)'; '*.*', 'All Files (*.*)'}, 'Select File name','put','.dat');
	if isequal(FileName,0),		return,		end
	fname = [PathName FileName];

	fid = fopen(fname, 'w');
	if (fid < 0),    errordlg(['Can''t open file:  ' fname],'Error');    return;     end
	fprintf(fid,'%f\t%f\t%f\n',[ud.data(1,:); ud.data(2,:); ud.data(3,:)]);
	fclose(fid);

% -------------------------------------------------------------------------------------------------
function push_integrate_CB(hObject, handles, opt)
	if (isempty(handles.Z1)),    errordlg('No grid loaded yet.','Error'),	return,		end
	if (isempty(opt))
		scale = 1;
	else
		scale = 980619.9203;    % Moritz's 1980 IGF value for gravity in mGal at 45 degrees latit
	end
	set(handles.figure1,'pointer','watch')
	if (~isa(handles.Z1, 'double'))		handles.Z1 = double(handles.Z1);	end
	if (get(handles.checkbox_leaveTrend,'Value'))       % Remove trend
		handles.Z1 = grdtrend_m(handles.Z1,handles.head_Z1,'-D','-N3');
	end
	[Z,band,k] = wavenumber_and_mboard(handles);
	k(1,1) = eps;       % Avoid a devided by zero warning
	Z = fft2(Z) ./ (k*scale);    Z(1,1) = 0;    clear k;
	Z = real(ifft2(Z));
	[Z,tmp] = unband(handles,Z,band);
	if (scale == 1),	tmp.name = 'Integrated grid';
	else				tmp.name = 'Geoid height';
	end
	set(handles.figure1,'pointer','arrow')
	mirone(single(Z),tmp);

% --------------------------------------------------------------------
function sectrumFun(handles, Z, head, opt1, Z2)
% OPT1 = 'Amplitude'   -> compute amplitude spectrum
% OPT1 = 'Power'       -> compute power spectrum
% OPT1 = 'Autocorr'    -> compute autocorrelation
% OPT1 = 'CrossPower'  -> compute cross power spectra between Z & Z2
% OPT1 = 'CrossCorrel' -> compute cross correlation between Z & Z2
% OPT1 = 'lpass|hpass' -> compute low pass|high pass filtering of Z
% Z2 -> needed when OPT1 = 'CrossPower' OR OPT1 = 'CrossCorrel'

	two_grids = 0;			image_type = 1;				% Used by the inverse transform case (default case)
	if (nargin == 5),		two_grids = 1;		end
	set(handles.figure1,'pointer','watch')
	if (~isa(Z, 'double'))		Z = double(Z);		end
	if (get(handles.checkbox_leaveTrend,'Value'))       % Remove trend
		Z = grdtrend_m(Z,handles.head_Z1,'-D','-N3');
		if (two_grids)
			if (~isa(Z2, 'double'))		Z2 = double(Z2);		end
			Z2 = grdtrend_m(Z2,handles.head_Z1,'-D','-N3');
		end
	end
	nx = handles.orig_ncols;        ny = handles.orig_nrows;
	[Z,band] = mboard(Z,nx,ny,handles.new_nx,handles.new_ny);
	if (two_grids)
		if (~isa(Z2, 'double'))		Z2 = double(Z2);		end
		Z2 = mboard(Z2,nx,ny,handles.new_nx,handles.new_ny);
	end
	m1 = band(1)+1;     m2 = m1 + ny - 1;
	n1 = band(3)+1;     n2 = n1 + nx - 1;

	if (any(strcmp(opt1,{'Amplitude' 'Power' 'CrossPower'})))   % For these this is more efficient
		Z = fftshift(fft2(Z));
		Z = Z(m1:m2,n1:n2);                 % Remove the padding skirt
		if (two_grids)
			Z2 = fftshift(fft2(Z2));        Z2 = Z2(m1:m2,n1:n2);
		end
	end

	if (strcmp(opt1,'Power'))
		Z = log10( (Z .* conj(Z) / (nx*ny)) + 1);       % An offset of 1 is added to avoid log(0)
		tmp.name = 'Power spectrum';
	elseif (strcmp(opt1,'CrossPower'))                  % Cross spectra
		Z = log10( (real(Z).*real(Z2) + imag(Z).*imag(Z2)) / (nx*ny) + 1);  % An offset of 1 is added to avoid log(0)
		tmp.name = 'Cross spectra';    clear Z2;
	elseif (strcmp(opt1,'Amplitude'))                   % Amplitude spectrum
		Z = log10(abs(Z) / (nx*ny) + 1);
		tmp.name = 'Amplitude spectrum';
	elseif (strcmp(opt1,'Autocorr'))                    % Autocorrelation
		Z = fftshift( single(real(ifft2(abs(fft2(Z)).^2))) );
		Z = Z(m1:m2,n1:n2);                             % Remove the padding skirt
		tmp.name = 'Autocorrelation';
	elseif (strcmp(opt1,'CrossCorrel'))                 % Cross Correlation
		Z = ifft2(abs(fft2(Z) .* fft2(Z2)));
		Z = fftshift(single(real(Z)));
		Z = Z(m1:m2,n1:n2);                             % Remove the padding skirt
		tmp.name = 'Cross Correlation';
	elseif ( any(strcmp(opt1,{'lpass' 'hpass'})) )		% Filtering
		handMir = guidata(handles.hMirFig);				% Handles of the space domain figure
		image_type = handMir.image_type;				% Used by the inverse transform case
		x = get(handMir.XXXhLine, 'XData');			y = get(handMir.XXXhLine, 'YData');
		hand = guidata(handMir.XXXhLine);				% Handles of the Spectrum figure
		xy_lims = getappdata(hand.axes1,'ThisImageLims');
		mask = ifftshift(img_fun('roipoly_j',xy_lims(1:2), xy_lims(3:4), Z, x, y));
		if (strcmp(opt1,'lpass')),	mask = ~mask;	end
		Z = fft2(Z);
		Z(mask) = 0;		% Do the filtering
		Z = single(real(ifft2(Z)));
		Z = Z(m1:m2,n1:n2);                             % Remove the padding skirt
		tmp.name = [upper(opt1(1))  'pass filtering'];
	end

	if ( (strcmp(opt1,'Autocorr') || strcmp(opt1,'CrossCorrel')) )
		delta_kx = 1;   delta_ky = 1;
	else				% make wave number array
		delta_kx = 2*pi / (handles.new_nx * handles.scaled_dx);
		delta_ky = 2*pi / (handles.new_ny * handles.scaled_dy);
	end
	if (handles.is_km)      % In km case scaled_dx|dy were in meters
		delta_kx = delta_kx * 1000;
		delta_ky = delta_ky * 1000;
	end

	nx2 = fix(nx/2);		ny2 = fix(ny/2);	sft_x = 0;		sft_y = 0;
	if (rem(nx,2) == 0),	sft_x = 1;			end
	if (rem(ny,2) == 0),	sft_y = 1;			end

	% This is currently only activated by the experimental Spector & Grant Filtering in R.A.S.
	if ( strcmp(opt1, 'bpass') )	% The bandpass (Do it here because only now that we know the freq vectors)
		x_lim = [-nx2 (nx2-sft_x)] * delta_kx / (2*pi);		y_lim = [-ny2 (ny2-sft_y)] * delta_ky / (2*pi);
		bp = handles.bandpass;			% Use a shorter name
		sqLow = [-bp(1) -bp(1); -bp(1) bp(1); bp(1) bp(1); bp(1) -bp(1); -bp(1) -bp(1)];
		maskLow = ifftshift(img_fun('roipoly_j',x_lim, y_lim, Z, sqLow(:,1), sqLow(:,2)));
		sqHigh = [-bp(4) -bp(4); -bp(4) bp(4); bp(4) bp(4); bp(4) -bp(4); -bp(4) -bp(4)];
		maskHigh = ~ifftshift(img_fun('roipoly_j',x_lim, y_lim, Z, sqHigh(:,1), sqHigh(:,2)));
		mask = maskLow | maskHigh;		% Combine Low and High cuts
		clear maskLow maskHigh
		% Need to add the cos tapper code to LP & HP
		Z = fft2(Z);
		Z(mask) = 0;
		Z = single(real(ifft2(Z)));
		Z = Z(m1:m2,n1:n2);                             % Remove the padding skirt
		handMir = guidata(handles.hMirFig);				% Handles of the space domain figure
	end

	if (isa(Z, 'double'))	Z = single(Z);		end
	if ( ~any(strcmp(opt1,{'lpass' 'hpass' 'bpass'})) )
		tmp.X = (-nx2:nx2-sft_x).*delta_kx;		tmp.Y = (-ny2:ny2-sft_y).*delta_ky;
		tmp.geog = 0;
	else
		tmp.X = getappdata(handMir.figure1,'dem_x');	tmp.Y = getappdata(handMir.figure1,'dem_y');
		if (isempty(tmp.X))			% We are operating on an image
			tmp.X = linspace(handMir.head(1), handMir.head(2), size(Z,2));
			tmp.Y = linspace(handMir.head(3), handMir.head(4), size(Z,1));
			Z = scaleto8(Z,-8);
		end
		tmp.name = 'FFT filtered';
		delta_kx = handMir.head(8);						delta_ky = handMir.head(9);
	end

	set(handles.figure1,'pointer','arrow')
	if (isa(Z,'uint8')),	tmp.cmap = gray(256);	end
	if (image_type ~= 2)
		tmp.head = [tmp.X(1) tmp.X(end) tmp.Y(1) tmp.Y(end) min(Z(:)) max(Z(:)) 0 delta_kx delta_ky];
		h = mirone(Z,tmp);
	else
		setappdata(0,'CropedColormap', gray(256));		% Force it to use a gray colormap
		h = mirone(Z);
		set(h,'Name','FFT filtered image')
	end

	if ( ~any(strcmp(opt1,{'lpass' 'hpass' 'bpass'})) )	% FFTs are never geogs
		handMir = guidata(h);		handMir.geog = 0;
		guidata(handMir.figure1, handMir)
	end
	if ( (strcmp(opt1,'Power') || strcmp(opt1,'Amplitude')) )
		setappdata(h, 'ParentFig', handles.hMirFig);	% Save the space domain Fig handle to eventual future use
	end

% --------------------------------------------------------------------
function push_goUDcont_CB(hObject, handles)
	if (isempty(handles.Z1)),    errordlg('No grid loaded yet.','Error');    return; end
	zup = str2double(get(handles.edit_UDcont,'String'));
	if (isnan(zup)),     return;     end
	set(handles.figure1,'pointer','watch')
	if (~isa(handles.Z1, 'double'))		handles.Z1 = double(handles.Z1);	end
	if (get(handles.checkbox_leaveTrend,'Value'))       % Remove trend
		handles.Z1 = grdtrend_m(handles.Z1,handles.head_Z1,'-D','-N3');
	end
	[Z,band,k] = wavenumber_and_mboard(handles);
	Z = fft2(Z) .* exp(-k.*zup);    clear k;
	Z = real(ifft2((Z)));
	[Z,tmp] = unband(handles,Z,band);
	tmp.name = 'U/D Continuation';
	set(handles.figure1,'pointer','arrow')
	mirone(single(Z),tmp);

% --------------------------------------------------------------------
function push_goDerivative_CB(hObject, handles, opt)
% Compute N-vertical derivative
	if (isempty(handles.Z1)),   errordlg('No grid loaded yet.','Error');    return; end
	if (isempty(opt)),	scale = 1;
	else				scale = 980619.9203;
	end					% Moritz's 1980 IGF value for gravity in mGal at lat = 45 degrees
	n_der = fix(str2double(get(handles.edit_derivative,'String')));
	set(handles.figure1,'pointer','watch')
	if (~isa(handles.Z1, 'double'))		handles.Z1 = double(handles.Z1);	end
	if (get(handles.checkbox_leaveTrend,'Value'))       % Remove trend
		handles.Z1 = grdtrend_m(handles.Z1,handles.head_Z1,'-D','-N3');
	end
	[Z,band,k] = wavenumber_and_mboard(handles,false,'taper');
	if (scale == 980619.9203)	n_der = 1;		end
	if (n_der > 1),  k = k.^n_der;   end
	Z = fft2(Z) .* k * scale;	Z(1,1) = 0;		clear k;
	Z = real(ifft2((Z)));
	[Z,tmp] = unband(handles,Z,band);
	if (scale == 1),	tmp.name = sprintf('%dth Vertical Derivative',n_der);
	else				tmp.name = 'Gravity anomaly (mGal)';
	end
	set(handles.figure1,'pointer','arrow')
	mirone(single(Z),tmp);

% --------------------------------------------------------------------
function push_AnalyticSig_CB(hObject, handles)
% Compute the 3D Analytic Signal
	if (isempty(handles.Z1)),    errordlg('No grid loaded yet.','Error'),	return,	end
	[Z,band,k] = wavenumber_and_mboard(handles,1);
	Z1 = fft2(Z);	Z1(1,1) = 0;				clear Z;
	dTdx = complex(-imag(Z1) .* k.x, real(Z1) .* k.x);
	dTdx = real(ifft2((dTdx))) .^2;
	dTdy = complex(-imag(Z1) .* k.y, real(Z1) .* k.y);
	dTdy = real(ifft2((dTdy))) .^2 + dTdx;		clear dTdx;
	k = sqrt(k.x .^ 2 + k.y .^ 2);
	dTdz = real(ifft2((Z1 .* k))) .^2 + dTdy;	clear Z1 dTdy k;
	dTdz = sqrt(dTdz);
	[Z,tmp] = unband(handles,dTdz,band);		clear dTdz;
	tmp.name = '3D Analytic Signal';
	set(handles.figure1,'pointer','arrow')
	mirone(single(Z),tmp);

% --------------------------------------------------------------------
function push_goDirDerivative_CB(hObject, handles)
% Compute a directional derivative
	if (isempty(handles.Z1)),    errordlg('No grid loaded yet.','Error'),	return,	end
	azim = str2double(get(handles.edit_dirDerivative,'String'));
	set(handles.figure1,'pointer','watch')
	if (~isa(handles.Z1, 'double'))		handles.Z1 = double(handles.Z1);	end
	if (get(handles.checkbox_leaveTrend,'Value'))       % Remove trend
		handles.Z1 = grdtrend_m(handles.Z1,handles.head_Z1,'-D','-N3');
	end
	[Z,band,k] = wavenumber_and_mboard(handles,1);
	fact = (sin(azim*pi/180) * k.x + cos(azim*pi/180) * k.y);   clear k;
	Z = fft2(Z);	Z(1,1) = 0;
	Z = complex(-imag(Z) .* fact, real(Z) .* fact);
	Z = real(ifft2((Z)));
	[Z,tmp] = unband(handles,Z,band);
	tmp.name = [num2str(azim) ' Azimuthal Derivative'];
	set(handles.figure1,'pointer','arrow')
	mirone(single(Z),tmp);

% --------------------------------------------------------------------
function [Z,band,k] = wavenumber_and_mboard(handles,opt,mode)
% Compute the wavenumber array and call mboard
% Z is the expanded array (to handles.new_nx & handles.new_ny)
% BAND is the second output of mboard
% K is the wavenumber array
% If opt = 1 K will be instead a structure with KX & KY fields
% MODE is either 'taper' (the default) or 'mirror'

	if (nargin == 1),	opt = false;	mode = 'taper';		end
	if (nargin == 2),	mode = 'taper';		end
	if (strcmp(mode, 'taper'))
		new_nx = handles.new_nx;			new_ny = handles.new_ny;
	else			% Mirror
		new_nx = handles.orig_ncols * 2;	new_ny = handles.orig_nrows * 2;
	end
	% calculate wavenumber array
	nx2 = fix(new_nx/2);		ny2 = fix(new_ny/2);

	if (rem(new_nx,2) == 0),	sft_x = 1;
	else						sft_x = 0;
	end
	if (rem(new_ny,2) == 0),	sft_y = 1;
	else						sft_y = 0;
	end

	%nx2 = handles.new_nx/2;     ny2 = handles.new_ny/2;    % Tivey way -> but and if n is odd?
	dkx = 2*pi / (new_nx * handles.scaled_dx);
	dky = 2*pi / (new_ny * handles.scaled_dy);
	%kx = (-nx2:nx2-1).*dkx;		ky = (-ny2:ny2-1).*dky;    % Tivey way
	kx = (-nx2:nx2-sft_x).*dkx;		ky = (-ny2:ny2-sft_y).*dky;
	%X = ones(size(ky))'*kx;		Y = ky'*ones(size(kx));
	X = repmat(kx,length(ky),1);	Y = repmat(ky',1,length(kx));
	if (opt)
		k.x = ifftshift(X);     clear X;
		k.y = ifftshift(Y);     clear Y;
	else
		k = ifftshift(sqrt(X.^2+Y.^2));      % wavenumber array   (Tivey used fftshift)
		clear X Y;
	end
	if (~isa(handles.Z1, 'double'))		handles.Z1 = double(handles.Z1);	end
	[Z,band] = mboard(handles.Z1,handles.orig_ncols,handles.orig_nrows,new_nx,new_ny, mode);

% --------------------------------------------------------------------
function [Z,hdr] = unband(handles,Z,band)
% Take the BAND vector and remove the padding skirt previously applyed
% (by mboard) to the Z matrix
% Return also the HDR structure used the call a new Mirone window
	if (~isempty(band))
		m1 = band(1)+1;     m2 = m1 + handles.orig_nrows - 1;
		n1 = band(3)+1;     n2 = n1 + handles.orig_ncols - 1;
	else					% The padding was made by replication
		m1 = 1;				m2 = handles.orig_nrows;
		n1 = 1;				n2 = handles.orig_ncols;
	end
	Z = Z(m1:m2,n1:n2);                 % Remove the padding skirt
	hdr.X = handles.X;		hdr.Y = handles.Y;
	hdr.head = handles.head_Z1;
	hdr.head(5) = min(min(Z));		hdr.head(6) = max(max(Z));

% --------------------------------------------------------------------
function popup_GridCoords_CB(hObject, handles)
	xx = get(hObject,'Value');
	if (xx == 1),		handles.geog = 1;		handles.is_meters = 0;  handles.is_km = 0;
	elseif (xx == 2)	handles.is_meters = 1;	handles.is_geog = 0;	handles.is_km = 0;
	elseif (xx == 3)
		handles.is_km = 1;      handles.is_geog = 0;    handles.is_meters = 0;
		handles.scaled_dx = handles.scaled_dx * 1000;
		handles.scaled_dy = handles.scaled_dy * 1000;
	end
	guidata(hObject,handles)

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
%a = 6378.137;   b = 6356.7521;
a = 6378137;    b = 6356752.1;      % Use meters

% Now, do the calculations...
e2 = 1 - (b*b)/(a*a);
sinlat = sin(orlat*pi/180);
denom  = sqrt(1 - e2 * sinlat * sinlat);
sclat = (pi/180) * a * (1 - e2)  / denom / denom / denom;
sclon = (pi/180) * a * cos(orlat*pi/180) / denom;

% --- Creates and returns a handle to the GUI figure. 
function fft_stuff_LayoutFcn(h1)

set(h1,'PaperUnits',get(0,'defaultfigurePaperUnits'),...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'MenuBar','none',...
'Name','FFT stuff',...
'NumberTitle','off',...
'Position',[520 579 570 221],...
'Resize','off',...
'Tag','figure1');

uicontrol('Parent',h1,'Position',[330 139 231 79],'String',{''},'Style','frame');
uicontrol('Parent',h1,'Position',[10 139 301 79],'String',{''},'Style','frame');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',{@fft_stuff_uiCB,h1,'edit_Grid1_CB'},...
'HorizontalAlignment','left',...
'Position',[50 187 231 22],...
'Style','edit',...
'TooltipString','Enter grid 1',...
'Tag','edit_Grid1');

uicontrol('Parent',h1, 'Position',[278 185 23 23],...
'Call',{@fft_stuff_uiCB,h1,'push_Grid1_CB'},...
'Tag','push_Grid1');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',{@fft_stuff_uiCB,h1,'edit_Grid2_CB'},...
'HorizontalAlignment','left',...
'Position',[50 157 231 22],...
'Style','edit',...
'TooltipString','Enter grid 2 (ONLY USED IN "CROSS" OPERATIONS)',...
'Tag','edit_Grid2');

uicontrol('Parent',h1, 'Position',[278 155 23 23],...
'Call',{@fft_stuff_uiCB,h1,'push_Grid2_CB'},...
'Tag','push_Grid2');

uicontrol('Parent',h1,'HorizontalAlignment','left',...
'Position',[17 160 31 15],'String','Grid2','Style','text');

uicontrol('Parent',h1,'HorizontalAlignment','left',...
'Position',[18 188 26 16],'String','Grid1','Style','text','Tag','text_FieldMag');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',{@fft_stuff_uiCB,h1,'listbox_nny_CB'},...
'Position',[390 147 51 60],'Style','listbox','Value',1,'Tag','listbox_nny');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',{@fft_stuff_uiCB,h1,'listbox_nnx_CB'},...
'Position',[498 150 51 60],'Style','listbox','Value',1,'Tag','listbox_nnx');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@fft_stuff_uiCB,h1,'edit_Nrows_CB'},...
'Position',[350 168 41 21],...
'Style','edit',...
'TooltipString','Number of grid rows',...
'Tag','edit_Nrows');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',{@fft_stuff_uiCB,h1,'edit_Ncols_CB'},...
'Position',[456 168 41 21],...
'Style','edit',...
'TooltipString','Number of grid columns',...
'Tag','edit_Ncols');

uicontrol('Parent',h1,'Position',[351 192 39 15],'String','# Rows', 'Style','text');
uicontrol('Parent',h1,'Position',[456 192 39 15],'String','# Cols', 'Style','text');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Position',[127 73 41 21],...
'Style','edit',...
'TooltipString','Height of continuation (in meters)',...
'Tag','edit_UDcont');

uicontrol('Parent',h1,'HorizontalAlignment','right','Position',[-5 75 130 17],...
'String','Up/Down Continuation','Style','text','HorizontalAlignment','right');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',{@fft_stuff_uiCB,h1,'edit_derivative_CB'},...
'Position',[127 43 41 21],...
'String','1',...
'Style','edit',...
'TooltipString','Order of vertical differentiation',...
'Tag','edit_derivative');

uicontrol('Parent',h1,'HorizontalAlignment','right','Position',[10 43 115 19],...
'String','Derivative (N-order)','Style','text', 'HorizontalAlignment','right');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',{@fft_stuff_uiCB,h1,'edit_dirDerivative_CB'},...
'Position',[127 13 41 21],...
'String','0',...
'Style','edit',...
'TooltipString','Angle of directional derivative',...
'Tag','edit_dirDerivative');

uicontrol('Parent',h1,'HorizontalAlignment','right','Position',[0 13 125 19],...
'String','Directional derivative','Style','text','HorizontalAlignment','right');

uicontrol('Parent',h1, 'Position',[168 73 31 21],...
'Call',{@fft_stuff_uiCB,h1,'push_goUDcont_CB'},...
'String','Go',...
'Tag','push_goUDcont');

uicontrol('Parent',h1, 'Position',[168 43 31 21],...
'Call',{@fft_stuff_uiCB3,h1,[],'push_goDerivative_CB'},...
'String','Go',...
'Tag','push_goDerivative');

uicontrol('Parent',h1, 'Position',[168 13 31 21],...
'Call',{@fft_stuff_uiCB,h1,'push_goDirDerivative_CB'},...
'String','Go',...
'Tooltip','Azimuth of directional derivative CCW from north',...
'Tag','push_goDirDerivative');

uicontrol('Parent',h1, 'Position',[220 109 101 21],...
'Call',{@fft_stuff_uiCB,h1,'push_powerSpectrum_CB'},...
'String','Power Spectrum',...
'Tag','push_powerSpectrum');

uicontrol('Parent',h1, 'Position',[220 76 101 21],...
'Call',{@fft_stuff_uiCB,h1,'push_autoCorr_CB'},...
'String','Auto Correlation',...
'Tag','push_autoCorr');

uicontrol('Parent',h1, 'Position',[220 43 101 21],...
'Call',{@fft_stuff_uiCB3,h1,[],'push_integrate_CB'},...
'String','Integrate',...
'Tag','push_integrate');

uicontrol('Parent',h1, 'Position',[220 10 101 21],...
'Call',{@fft_stuff_uiCB3,h1,'FAA2Geoid','push_integrate_CB'},...
'String','FAA2Geoid',...
'Tag','push_faa2geoid');

uicontrol('Parent',h1, 'Position',[340 109 121 21],...
'Call',{@fft_stuff_uiCB,h1,'push_crossSpectra_CB'},...
'String','Cross Spectra',...
'Tag','push_crossSpectra');

uicontrol('Parent',h1, 'Position',[340 76 121 21],...
'Call',{@fft_stuff_uiCB,h1,'push_crossCorr_CB'},...
'String','Cross Correlation',...
'Tag','push_crossCorr');

uicontrol('Parent',h1, 'Position',[340 43 121 21],...
'Call',{@fft_stuff_uiCB,h1,'push_radialPowerAverage_CB'},...
'String','Radial power average',...
'Tag','push_radialPowerAverage');

uicontrol('Parent',h1, 'Position',[340 10 121 21],...
'Call',{@fft_stuff_uiCB3,h1,'Geoid2FAA','push_goDerivative_CB'},...
'String','Geoid2FAA',...
'Tag','push_geoid2faa');

uicontrol('Parent',h1, 'Position',[470 109 100 21],...
'Call',{@fft_stuff_uiCB,h1,'push_AnalyticSig_CB'},...
'String','Analytic Signal',...
'Tag','push_AnalyticSig');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',{@fft_stuff_uiCB,h1,'popup_GridCoords_CB'},...
'Position',[86 106 82 22],...
'String',{'Geogs'; 'Meters'; 'Kilometers'},...
'Style','popupmenu',...
'Tooltip','GRID COORDINATES: IT IS YOUR RESPONSABILITY THAT THIS IS CORRECT',...
'Value',1,...
'Tag','popup_GridCoords');

uicontrol('Parent',h1,'ForegroundColor',[1 0 0],'Position',[27 109 58 15],...
'String','CONFIRM','Style','text');

uicontrol('Parent',h1, 'Position',[50 140 105 15],...
'String','Remove trend',...
'Style','checkbox',...
'Tooltip','If checked remove a plane before transformations',...
'Tag','checkbox_leaveTrend');

function fft_stuff_uiCB(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
	feval(callback_name,hObject,guidata(h1));

function fft_stuff_uiCB3(hObject, eventdata, h1, opt, callback_name)
% This function is executed by the callback and than the handles is allways updated.
	feval(callback_name,hObject,guidata(h1),opt);

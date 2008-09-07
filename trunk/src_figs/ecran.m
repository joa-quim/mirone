function varargout = ecran(varargin)

%	Copyright (c) 2004-2008 by J. Luis
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

%#function select_cols

	hObject = figure('Tag','figure1','Visible','off');
	ecran_LayoutFcn(hObject);
	handles = guihandles(hObject);
	set(hObject,'RendererMode','auto')
	movegui(hObject,'east');

	global home_dir

	if (~isempty(home_dir)),		handles.home_dir = home_dir;
	else							handles.home_dir = cd;
	end
	handles.d_path = [handles.home_dir filesep 'data' filesep];
	load([handles.d_path 'mirone_pref.mat']);
	try,		handles.last_dir = directory_list{1};
	catch		handles.last_dir = cd;
	end

	if (~nargin),   varargin(1) = {[]};   end

	handles.handMir = [];
	if ( ~isempty(varargin{1}) && isa(varargin{1},'struct') )
		handles.handMir = varargin{1};
		varargin{1} = 'Image';			% For backward compatibility sake
		handles.work_dir = handles.handMir.work_dir;
	else
		handles.work_dir = handles.last_dir;
	end

	% Load some icons from mirone_icons.mat
	load([handles.d_path 'mirone_icons.mat'],'zoom_ico','zoomx_ico', 'clipcopy_ico');
	link_ico = imread([handles.d_path 'link.png']);

	hTB = uitoolbar('parent',hObject,'Clipping', 'on', 'BusyAction','queue','HandleVisibility','on',...
		'Interruptible','on','Tag','FigureToolBar','Visible','on');
	uitoggletool('parent',hTB,'Click',{@zoom_CB,''}, 'cdata',zoom_ico,'Tooltip','Zoom');
	uitoggletool('parent',hTB,'Click',{@zoom_CB,'x'}, 'cdata',zoomx_ico,'Tooltip','Zoom X');
	if (strncmp(computer,'PC',2))
		uipushtool('parent',hTB,'Click',@copyclipbd_CB, 'cdata',clipcopy_ico,'Tooltip', 'Copy to clipboard ','Sep','on');
	end
	if (~isempty(handles.handMir))
		uitoggletool('parent',hTB,'Click',@pick_CB, 'cdata',link_ico,'Tooltip', ...
			'Pick data point in curve and plot it the mirone figure','Sep','on');
	end
	uitoggletool('parent',hTB,'Click',@isocs_CB, 'Tooltip','Enter ages & plot a geomagnetic barcode','Sep','on');

	handles.n_plot = 0;         % Counter of the number of lines. Used for line color painting

	if (strcmp(varargin{1},'reuse') && nargin < 3)
		errordlg('Error calling ecran: Minimum arguments are "type",X,Y','Error');      return
	end
	if ( strcmp(varargin{1},'Image') && nargin < 5 )   
		errordlg('Error calling ecran: Minimum arguments are "type",X,Y,Z','Error');    return
	end

	handles.ageStart = 0;		handles.ageEnd = nan;
	handles.dist = [];			% It will contain cumulated distance if input is (x,y,z)
	handles.hLine = [];			% Handles to the ploted line
	handles.polyFig = [];		% Handles to the (eventual) figure for trend1d polyfit

	axes(handles.axes1)		% Make it the active one

	% Choose the default ploting mode
	if isempty(varargin{1})          % When the file will be read latter
		set([handles.checkbox_geog handles.popup_selectPlot handles.popup_selectSave], 'Visible','off')	% Hide those
		set(hObject,'Name','XY view')
	
	elseif strcmp(varargin{1},'Image')
		handles.data(:,1) = varargin{2};    handles.data(:,2) = varargin{3};    handles.data(:,3) = varargin{4};
		set(handles.popup_selectSave,'String',{'Save Profile on disk';'Distance,Z (data units -> ascii)';
			'Distance,Z (data units -> binary)';'X,Y,Z (data units -> ascii)';'X,Y,Z (data units -> binary)';
			'Distance,Z (data units -> mat file)'});
		rd = dist_along_profile(handles.data(:,1), handles.data(:,2));
		handles.dist = rd;				% This one is by default, so save it in case user wants to save it to file
		handles.hLine = plot(rd,handles.data(:,3));		axis(handles.axes1,'tight');
		set(hObject,'Name',varargin{5})
	
	elseif strcmp(varargin{1},'reuse')                      % Case of auto-referenced call
		varargin(nargin+1:9) = cell(1,9-nargin);			% So that varargin{1:9} allways exists.
		set([handles.checkbox_geog handles.popup_selectPlot handles.popup_selectSave], 'Visible','off')	% Hide those
		handles.data(:,1) = varargin{2};        handles.data(:,2) = varargin{3};
		if ~isempty(varargin{9}) && strcmp(varargin{9},'semilogy')
			handles.hLine = semilogy(handles.data(:,1),handles.data(:,2));
		elseif ~isempty(varargin{9}) && strcmp(varargin{9},'semilogx')
			handles.hLine = semilogx(handles.data(:,1),handles.data(:,2));
		else
			handles.hLine = plot(handles.data(:,1),handles.data(:,2));
		end
		axis(handles.axes1,'tight');
		
		if ~isempty(varargin{5}),    set(hObject,'Name',varargin{5});	end		% Figure Name
		if ~isempty(varargin{6}),    xlabel(varargin{6});				end		% XLabel
		if ~isempty(varargin{7}),    ylabel(varargin{7});				end		% YLabel
		if ~isempty(varargin{8}),    title(varargin{8});				end		% Title
	end

	guidata(hObject, handles);
	set(hObject,'Visible','on');
	if (nargout),	varargout{1} = hObject;		end

% --------------------------------------------------------------------------------------------------
function zoom_CB(obj,eventdata,opt)
	if (strcmp(get(obj,'State'),'on'))
		if (strcmp(opt,'x')),		zoom_j xon;
		else						zoom_j on;
		end
	else
		zoom_j off;
	end

% --------------------------------------------------------------------------------------------------
function copyclipbd_CB(obj,eventdata)
% Copy curve(s) to the Windows clipboard.
	handles = guidata(obj);
	if (isempty(handles.hLine)),	return,		end
	hLines = findobj(handles.axes1,'type','line');
	x = get(hLines(1), 'xdata');	y = get(hLines, 'ydata');
	if (~isa(y,'cell'))
		data = [x(:) y(:)];
	else
		n_pts = numel(x);
		nys = zeros(numel(y),1);
		for (k = 1:numel(y)),		nys(k) = numel(y{k});	end		% Count each curve number of points
		ind = find(nys == nys(1));			% See if they are all of the same size
		if (numel(ind) == numel(nys))		% All curves have the same number of points. Good.
			data = x(:);
			for (k = 1:numel(y))			% Loop over number of curves
				data = [data y{k}(:)];		% Slower due to non-preallocation but safe
			end
		else
			% Search for the first curve that has as many points as the X vector 
			for (k = 1:numel(y))			% Loop over number of curves
				if (numel(y{k}) == n_pts)
					data = [x(:) y{k}(:)];	% Found one (We must have at least one)
					break
				end
			end
			warndlg('Several curves but not all of the same size. Only one was copyied.','Warning')
		end
	end
	mat2clip(data)

% --------------------------------------------------------------------------------------------------
function pick_CB(obj,eventdata)
	handles = guidata(obj);
	if (strcmp(get(obj,'State'),'on'))
		set(handles.figure1,'WindowButtonDownFcn',@add_MarkColor)
	else
		set(handles.figure1,'WindowButtonDownFcn','')	
	end

% --------------------------------------------------------------------
function add_MarkColor(obj, eventdata)
% Add a red Marker over the closest (well, near closest) clicked point.

	handles = guidata(obj);     % get handles

	button = get(handles.figure1, 'SelectionType');
	if (~strcmp(button,'normal')),		return,		end     % Accept only left-clicks

	pt = get(handles.axes1, 'CurrentPoint');

    hM = findobj(handles.figure1,'Type','Line','tag','marker');
	x = handles.dist;			y = handles.data(:,3);

	x_lim = get(handles.axes1,'XLim');		y_lim = get(handles.axes1,'YLim');
	dx = diff(x_lim) / 20;					% Search only betweem +/- 1/10 of x_lim
	id = (x < (pt(1,1)-dx) | x > (pt(1,1)+dx));
	x(id) = [];             y(id) = [];     % Clear outside-2*dx points to speed up the search code
	x_off = find(~id);		x_off = x_off(1);	% Get the index of the first non killed x element
	XScale = diff(x_lim);	YScale = diff(y_lim);

	r = sqrt(((pt(1,1)-x) ./ XScale).^2 + ((pt(1,2)-y) ./ YScale).^2);
	[temp,i] = min(r);
	pt_x = x(i);                pt_y = y(i);
	
	xr = get(hM,'XData');       yr = get(hM,'YData');
	id = find(xr == pt_x);
	if (isempty(id))            % New Marker
		if (~isempty(handles.handMir))
			% Get the X,Y coordinates to plot this point in the Mirone figure
			% We need to add x_off since "i" counts only inside the +/- 1/10 of x_lim centered on current point
			mir_pt_x = handles.data(i+x_off-1,1);	mir_pt_y = handles.data(i+x_off-1,2);
			h = line(mir_pt_x, mir_pt_y,'Parent',handles.handMir.axes1,'Marker','o','MarkerFaceColor','k', ...
				'MarkerSize',6,'LineStyle','none', 'Tag','picked');
			draw_funs(h,'DrawSymbol')		% Set uicontexts
		end

		if (isempty(hM))        % First red Marker on this axes
			line(pt_x, pt_y,'Marker','s','MarkerFaceColor','r','MarkerSize',5,'LineStyle','none','Tag','marker','UserData',h);
		else
			xr = [xr pt_x];     yr = [yr pt_y];
			set(hM,'XData',xr, 'YData', yr)
			ud = get(hM,'UserData');
			set(hM,'UserData', [ud h]);		% Save the Mirone symbol handle here
		end
	else                        % Marker already exists. Kill it
		xr(id) = [];            yr(id) = [];
		set(hM,'XData',xr, 'YData', yr)
		ud = get(hM,'UserData');
		try,	delete(ud(id)),		end
		ud(id) = [];
		set(hM,'UserData', ud)
	end

% --------------------------------------------------------------------------------------------------
function isocs_CB(obj,eventdata)
	handles = guidata(obj);
	if (strcmp(get(obj,'State'),'on'))
		set([handles.checkbox_geog handles.popup_selectPlot handles.popup_selectSave], 'Visible','off')
		set([handles.edit_ageStart handles.edit_ageEnd handles.push_magBar handles.text_ageStart ...
				handles.text_ageEnd], 'Visible','on')
	else
		set([handles.checkbox_geog handles.popup_selectPlot handles.popup_selectSave], 'Visible','on')
		set([handles.edit_ageStart handles.edit_ageEnd handles.push_magBar handles.text_ageStart ...
				handles.text_ageEnd], 'Visible','off')
	end

% -------------------------------------------------------------------------------
function checkbox_geog_Callback(hObject, eventdata, handles)
	if get(hObject,'Value')
		set(handles.popup_selectPlot,'String',{'Distance along profile (data units)';
		'Distance along profile (km)';'Distance along profile (NM)'});
		set(handles.popup_selectSave,'String',{'Save Profile on disk';'Distance,Z (data units -> ascii)';
			'Distance,Z (data units -> binary)';'Distance,Z (km -> ascii)';'Distance,Z (km -> binary)';
			'Distance,Z (NM -> ascii)';'Distance,Z (NM -> binary)';
			'X,Y,Z (data units -> ascii)';'X,Y,Z (data units -> binary)';
			'Distance,Z (NM -> mat file)';'Distance,Z (km -> mat file)';'Distance,Z (data units -> mat file)'});
	else
		set(handles.popup_selectPlot,'String','Distance along profile (data units)','Value',1);
		set(handles.popup_selectSave,'String',{'Save Profile on disk';'Distance,Z (data units -> ascii)';
			'Distance,Z (data units -> binary)';'X,Y,Z (data units -> ascii)';'X,Y,Z (data units -> binary)';
			'Distance,Z (data units -> mat file)'},'Value',1);
	end

% ---------------------------------------------------------------------------------
function popup_selectPlot_Callback(hObject, eventdata, handles)
	val = get(hObject,'Value');     str = get(hObject, 'String');
	D2R = pi/180;
	h = findobj(handles.axes1,'Type','line');
	switch str{val};
	
		case 'Distance along profile (data units)'  % Compute the accumulated distance along profile in data units
			xd = handles.data(:,1);			yd = handles.data(:,2);
	
		case 'Distance along profile (km)'           % Compute the accumulated distance along profile in km
            deg2km = 111.1949;
			xd = (handles.data(:,1)*deg2km) .* cos(handles.data(:,2)*D2R);		yd = handles.data(:,2) * deg2km;
	
		case 'Distance along profile (NM)'            % Compute the accumulated distance along profile in Nmiles
            deg2nm = 60.04;
			xd = (handles.data(:,1)*deg2nm) .* cos(handles.data(:,2)*D2R);		yd = handles.data(:,2) * deg2nm;
	end
	rd = dist_along_profile(xd, yd);
	set(handles.hLine,'XData',rd);				axis tight;
	guidata(hObject, handles);

% ---------------------------------------------------------------------------------
function popup_selectSave_Callback(hObject, eventdata, handles)
val = get(hObject,'Value');     str = get(hObject, 'String');
D2R = pi/180;
deg2km = 111.1949;		deg2nm = 60.04;
switch str{val};
    case 'Save Profile on disk'                    %
    case 'Distance,Z (data units -> ascii)'                    % Save profile in ascii data units
		[FileName,PathName] = put_or_get_file(handles,{'*.dat', 'Dist Z (*.dat)'; '*.*', 'All Files (*.*)'},'Distance,Z (ascii)','put','.dat');
		if isequal(FileName,0),		set(hObject,'Value',1),		return,		end     % User gave up
		double2ascii([PathName FileName],[handles.dist handles.data(:,3)],'%f\t%f');
    case 'Distance,Z (data units -> binary)'                    % Save profile in binary data units
		[FileName,PathName] = put_or_get_file(handles,{'*.dat', 'Dist Z (*.dat)'; '*.*', 'All Files (*.*)'},'Distance,Z (binary float)','put');
		if isequal(FileName,0),		set(hObject,'Value',1),		return,		end     % User gave up
		fid = fopen([PathName FileName],'wb');
		fwrite(fid,[handles.dist handles.data(:,3)]','float');  fclose(fid);
    case 'Distance,Z (km -> ascii)'                    % Save profile in ascii (km Z) 
		[FileName,PathName] = put_or_get_file(handles,{'*.dat', 'Dist Z (*.dat)'; '*.*', 'All Files (*.*)'},'Distance (km),Z (ascii)','put','.dat');
		if isequal(FileName,0),		set(hObject,'Value',1),		return,		end     % User gave up
		xd = (handles.data(:,1)*deg2km) .* cos(handles.data(:,2)*D2R);		yd = handles.data(:,2) * deg2km;
		rd = dist_along_profile(xd, yd);
		double2ascii([PathName FileName],[rd handles.data(:,3)],'%f\t%f')
    case 'Distance,Z (km -> binary)'                    % Save profile in binary (km Z) 
		[FileName,PathName] = put_or_get_file(handles,{'*.dat', 'Dist Z (*.dat)'; '*.*', 'All Files (*.*)'},'Distance (km),Z (binary float)','put');
		if isequal(FileName,0),		set(hObject,'Value',1),		return,		end     % User gave up
		xd = (handles.data(:,1)*deg2km) .* cos(handles.data(:,2)*D2R);		yd = handles.data(:,2) * deg2km;
		rd = dist_along_profile(xd, yd);
		fid = fopen([PathName FileName],'wb');
		fwrite(fid,[rd handles.data(:,3)]','float');  fclose(fid);
    case 'Distance,Z (NM -> ascii)'                    % Save profile in ascii (NM Z) 
		[FileName,PathName] = put_or_get_file(handles,{'*.dat', 'Dist Z (*.dat)'; '*.*', 'All Files (*.*)'},'Distance (m),Z (ascii)','put','.dat');
		if isequal(FileName,0),		set(hObject,'Value',1),		return,		end     % User gave up
		xd = (handles.data(:,1)*deg2nm) .* cos(handles.data(:,2)*D2R);		yd = handles.data(:,2) * deg2nm;
		rd = dist_along_profile(xd, yd);
        double2ascii([PathName FileName],[rd handles.data(:,3)],'%f\t%f')
    case 'Distance,Z (NM -> binary)'                    % Save profile in binary (NM Z) 
		[FileName,PathName] = put_or_get_file(handles,{'*.dat', 'Dist Z (*.dat)'; '*.*', 'All Files (*.*)'},'Distance (m),Z (binary float)','put');
		if isequal(FileName,0),		set(hObject,'Value',1),		return,		end     % User gave up
		xd = (handles.data(:,1)*deg2nm) .* cos(handles.data(:,2)*D2R);		yd = handles.data(:,2) * deg2nm;
		rd = dist_along_profile(xd, yd);
		fid = fopen([PathName FileName],'wb');
		fwrite(fid,[rd handles.data(:,3)]','float');  fclose(fid);
    case 'X,Y,Z (data units -> ascii)'						% Save profile in ascii (km Z) 
		[FileName,PathName] = put_or_get_file(handles,{'*.dat', 'x,y,z (*.dat)';'*.*', 'All Files (*.*)'},'X,Y,Z (ascii)','put','.dat');
		if isequal(FileName,0),		set(hObject,'Value',1),		return,		end     % User gave up
		double2ascii([PathName FileName],[handles.data(:,1) handles.data(:,2) handles.data(:,3)],'%f\t%f\t%f')
    case 'X,Y,Z (data units -> binary)'						% Save profile in binary (km Z) 
		[FileName,PathName] = put_or_get_file(handles,{'*.dat', 'x,y,z (*.dat)';'*.*', 'All Files (*.*)'},'X,Y,Z (binary float)','put');
		if isequal(FileName,0),		set(hObject,'Value',1),		return,		end     % User gave up
		fid = fopen([PathName FileName],'wb');
		fwrite(fid,[handles.data(:,1) handles.data(:,2) handles.data(:,3)]','float');  fclose(fid);
    case 'Distance,Z (NM -> mat file)'						% Save profile in mat file (m Z) 
		[FileName,PathName] = put_or_get_file(handles,{'*.mat', 'Dist Z (*.mat)';'*.*', 'All Files (*.*)'},'Distance (m),Z (Matlab mat file)','put');
		if isequal(FileName,0),		set(hObject,'Value',1),		return,		end     % User gave up
		xd = (handles.data(:,1)*deg2nm) .* cos(handles.data(:,2)*D2R);		yd = handles.data(:,2) * deg2nm;
		rd = dist_along_profile(xd, yd);
		R = rd;		Z = handles.data(:,3);					% More one BUG, handles.data(:,3) canot be saved
		save([PathName FileName],'R','Z')
    case 'Distance,Z (km -> mat file)'						% Save profile in mat file (km Z)
		[FileName,PathName] = put_or_get_file(handles,{'*.mat', 'Dist Z (*.mat)';'*.*', 'All Files (*.*)'},'Distance (km),Z (Matlab mat file)','put');
		if isequal(FileName,0),		set(hObject,'Value',1),		return,		end     % User gave up
		xd = (handles.data(:,1)*deg2km) .* cos(handles.data(:,2)*D2R);		yd = handles.data(:,2) * deg2km;
		rd = dist_along_profile(xd, yd);
		R = rd;		Z = handles.data(:,3);
		save([PathName FileName],'R','Z')
    case 'Distance,Z (data units -> mat file)'				% Save profile in binary data units
		[FileName,PathName] = put_or_get_file(handles,{'*.mat', 'Dist Z (*.mat)';'*.*', 'All Files (*.*)'},'Distance,Z (Matlab mat file)','put');
		if isequal(FileName,0),		set(hObject,'Value',1),		return,		end     % User gave up
		R = handles.dist';   Z = handles.data(:,3)';      % More one BUG, handles.data(:,3) canot be saved
		save([PathName FileName],'R','Z')
end
set(hObject,'Value',1);

% --------------------------------------------------------------------
function FileExport_Callback(hObject, eventdata, handles)
	filemenufcn(handles.figure1,'FileExport')

% --------------------------------------------------------------------
function FilePrint_Callback(hObject, eventdata, handles)
	if (ispc),		print -v
	else			print
	end

% --------------------------------------------------------------------
function FileOpen_Callback(hObject, eventdata, handles)
% Read the file and select what columns to plot
	str1 = {'*.dat;*.DAT;', 'Data files (*.dat,*.DAT)';'*.*', 'All Files (*.*)'};
	[FileName,PathName] = put_or_get_file(handles,str1,'Select input data','get');
	if isequal(FileName,0),		return,		end	
	fname = [PathName FileName];

	data = text_read(fname,NaN);
	if (isempty(data) || size(data,2) == 1)
		errordlg('File doesn''t have any recognized nymeric data (Quiting) or one column only.','Error');    return
	end

	% If msgbox exist we have to move it from behind the main window. So get it's handle
	hMsgFig = gcf;
	if (handles.figure1 ~= hMsgFig)
		figure(hMsgFig);		% If error msgbox exists, bring it forward
		% Here we have a stupid problem. If don't kill the message window before the
		% select_cols is called this later wont work. FDS I have no more patiente for this.
		pause(0.5)
		try		delete(hMsgFig),		end
	end

% 	if (size(data,2) > 2)			% Not yet mature
% 		out = select_cols(data,'xyz',fname,1000);
% 	else
% 		out = select_cols(data,'xy',fname,1000);
% 	end
	out = select_cols(data,'xy',fname,1000);
	if (isempty(out)),		return,		end

	if (numel(out) == 4)			% Not yet in use
		rd = dist_along_profile( data(:,out(1)), data(:,out(2)) );
		handles.dist = rd;				% Save it in case user wants to save it to file
		handles.hLine = line('Parent',handles.axes1,'XData',rd, 'YData',data(:,out(3)));
	else
		handles.hLine = line('Parent',handles.axes1,'XData',data(:,out(1)),'YData',data(:,out(2)));
	end
	axis(handles.axes1,'tight');
	handles.n_plot = handles.n_plot + 1;
	if (handles.n_plot > 1)
		c_order = get(handles.axes1,'ColorOrder');
		if (handles.n_plot <= 7)
			nc = handles.n_plot;
		else
			nc = rem(handles.n_plot,7);     % recycle through the default colors
		end
		cor = c_order(nc,:);
		set(handles.hLine,'Color',cor)
	end
	handles.data = [data(:,out(1)) data(:,out(2))];     % NOTE, if handles.n_plot > 1 only last data is saved
	guidata(hObject, handles);

% --------------------------------------------------------------------
function FileSave_Callback(hObject, eventdata, handles)
	if (isempty(handles.hLine)),	return,		end
	[FileName,PathName] = put_or_get_file(handles,{'*.dat', 'X,Y (*.dat)';'*.*', 'All Files (*.*)'},'X,Y (ascii)','put', '.dat');
	if isequal(FileName,0),		return,		end     % User gave up	
	xx = get(handles.hLine,'XData');		yy = get(handles.hLine,'YData');
	double2ascii([PathName FileName],[xx(:) yy(:)],'%f\t%f');

% --------------------------------------------------------------------
function AnalysisFFT_AmpSpectrum_Callback(hObject, eventdata, handles)
	if (isempty(handles.hLine)),	return,		end
	x = get(handles.hLine,'XData');
	Fs = 1 / (x(2) - x(1));			% Sampling frequency
	Fn=Fs/2;						% Nyquist frequency
	x = get(handles.hLine,'YData');
	NFFT=2.^(ceil(log(length(x))/log(2)));		% Next highest power of 2 greater than or equal to length(x)
	FFTX=fft(x,NFFT);							% Take fft, padding with zeros, length(FFTX)==NFFT
	NumUniquePts = ceil((NFFT+1)/2);
	FFTX=FFTX(1:NumUniquePts);					% fft is symmetric, throw away second half
	MX=abs(FFTX);								% Take magnitude of X
	% Multiply by 2 to take into account the fact that we threw out second half of FFTX above
	MX=MX*2;				MX(1)=MX(1)/2;		% Account for endpoint uniqueness
	MX(length(MX))=MX(length(MX))/2;			% We know NFFT is even
	% Scale the FFT so that it is not a function of the length of x.
	MX=MX/length(x);        f=(0:NumUniquePts-1)*2*Fn/NFFT;
	ecran('reuse',f,MX,[],'Amplitude Spectrum','frequency',[],'Amplitude Spectrum')

% --------------------------------------------------------------------
function AnalysisFFT_PSD_Callback(hObject, eventdata, handles)
	if (isempty(handles.hLine)),	return,		end
	x = get(handles.hLine,'XData');			Fs = 1 / (x(2) - x(1));		% Sampling frequency
	x = get(handles.hLine,'YData');
	[Pxx,w] = psd(x,Fs);
	% We want to guarantee that the result is an integer if X is a negative power of 10.
	% To do so, we force some rounding of precision by adding 300-300.
	Pxx = (10.*log10(Pxx)+300)-300;    % Compute db
	ecran('reuse',w,Pxx,[],'Power Spectrum','Frequency (Hz)','Power Spectral Density (dB/Hz)','Periodogram PSD Estimate')

% --------------------------------------------------------------------
function AnalysisAutocorrelation_Callback(hObject, eventdata, handles)
	if (isempty(handles.hLine)),	return,		end
	xx = get(handles.hLine,'XData');		yy = get(handles.hLine,'YData');
	c = autocorr(yy);				n = length(yy);
	ecran('reuse',xx,c(n:end),[],'Normalized Autocorrelation','Lag in user X units')

% --------------------------------------------------------------------
function AnalysisRemoveMean_Callback(hObject, eventdata, handles)
	if (isempty(handles.hLine)),	return,		end
	xx = get(handles.hLine,'XData');		yy = get(handles.hLine,'YData');
	ecran('reuse',xx,yy-mean(yy),[],'Mean Removed')

% --------------------------------------------------------------------
function AnalysisRemoveTrend_Callback(hObject, eventdata, handles)
	if (isempty(handles.hLine)),	return,		end
	xx = get(handles.hLine,'XData');		yy = get(handles.hLine,'YData');
	p = polyfit(xx,yy,1);			y = polyval(p,xx);
	ecran('reuse',xx,yy-y,[],'Trend Removed')

% --------------------------------------------------------------------
function AnalysisFitPoly_Callback(hObject, eventdata, handles)
	if (isempty(handles.hLine)),	return,		end
	xx = get(handles.hLine,'XData');		yy = get(handles.hLine,'YData');
	handles.polyFig = ecran_trend1d(handles.axes1, [xx(:) yy(:)]);
	guidata(handles.figure1, handles)

% --------------------------------------------------------------------
function AnalysisSmoothSpline_Callback(hObject, eventdata, handles)
	if (isempty(handles.hLine)),	return,		end
	xx = get(handles.hLine,'XData');		yy = get(handles.hLine,'YData');
	[pp,p] = spl_fun('csaps',xx,yy);		% This is just to get csaps's p estimate
	y = spl_fun('csaps',xx,yy,p,xx);
	hold on;	h = plot(xx,y);		hold off;

	smoothing_param(p, [xx(1) xx(2)-xx(1) xx(end)], handles.figure1, handles.axes1, handles.hLine, h);
	guidata(hObject, handles);

% --------------------------------------------------------------------
function Analysis1derivative_Callback(hObject, eventdata, handles)
	if (isempty(handles.hLine)),	return,		end
	xx = get(handles.hLine,'XData');		yy = get(handles.hLine,'YData');
	pp = spl_fun('csaps',xx,yy,1);		% Use 1 for not smoothing, just interpolate
	v = spl_fun('ppual',pp,xx,'l','first');
	ecran('reuse',xx,v,[],'First derivative')

% --------------------------------------------------------------------
function Analysis2derivative_Callback(hObject, eventdata, handles)
	if (isempty(handles.hLine)),	return,		end
	xx = get(handles.hLine,'XData');		yy = get(handles.hLine,'YData');
	pp = spl_fun('csaps',xx,yy,1);		% Use 1 for not smoothing, just interpolate
	v = spl_fun('ppual',pp,xx,'l','second');
	ecran('reuse',xx,v,[],'Second derivative')

% --------------------------------------------------------------------
function edit_ageStart_Callback(hObject, eventdata, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx))
		set(hObject,'String','')
	elseif (xx < 0)
		xx = 0;
		set(hObject,'String',xx)
	end
	handles.ageStart = xx;
	guidata(handles.figure1, handles)

% --------------------------------------------------------------------
function edit_ageEnd_Callback(hObject, eventdata, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx))
		set(hObject,'String','')
	end
	handles.ageEnd = xx;
	guidata(handles.figure1, handles)

% --------------------------------------------------------------------
function push_magBar_Callback(hObject, eventdata, handles)
	if (isnan(handles.ageStart) || isnan(handles.ageEnd))
		errordlg('Take a second look to what you are asking for. Wrong ages','Error'),		return
	end

	set(handles.axes2, 'Vis', 'on','XTick',[], 'YTick',[])

	reverse_XDir = false;		% If first age > last age, we'll revert the sense of the X axis
	if ( handles.ageStart >= handles.ageEnd )
		reverse_XDir = true;
		tmp = handles.ageStart;		handles.ageStart = handles.ageEnd;		handles.ageEnd = tmp;
		set(handles.axes2,'XDir','reverse')
	end

	chron_file = [handles.d_path 'Cande_Kent_95.dat'];
	fid = fopen(chron_file,'r');
	todos = fread(fid,'*char');     [chron age_start age_end age_txt] = strread(todos,'%s %f %f %s');
	fclose(fid);    clear todos

	id_ini = (age_start >= handles.ageStart);		id_ini = find(id_ini);		id_ini = id_ini(1);
	id_fim = (age_start <= handles.ageEnd);			id_fim = find(~id_fim);		id_fim = id_fim(1) - 1;
	age_start = age_start(id_ini:id_fim);
	age_end = age_end(id_ini:id_fim);
	age_txt = age_txt(id_ini:id_fim);

	% Take care of end ages which certainly do not cuincide with what was asked
	if (age_start(1) > handles.ageStart)
		age_start = [handles.ageStart; age_start];	age_end = [handles.ageStart; age_end];
		age_txt = cat(1,{'a'}, age_txt(:));
	end
	if (age_end(end) < handles.ageEnd)
		age_start(end+1) = handles.ageEnd;		age_end(end+1) = handles.ageEnd;		age_txt{end+1} = 'a';
	end

	x = [age_start'; age_end'];
	x = x(:);
	set(handles.axes2, 'xlim', [age_start(1) age_end(end)])
	y = [zeros(numel(x),1); ones(numel(x),1)]; 
	x = [x; x(end:-1:1)];

	n_ages = numel(age_start);
	n2 = 2 * n_ages;
	c1 = (1:n2-1)';     c3 = n2*2 - c1;
	c2 = c3 + 1;        c4 = c1 + 1;
	faces = [c1 c2 c3 c4];

	cor = repmat([0 0 0; 1 1 1],n_ages-1,1);    cor = [cor; [0 0 0]];
	hp = patch('Parent',handles.axes2,'Faces',faces,'Vertices',[x y],'FaceVertexCData',cor,'FaceColor','flat');   
	set(handles.figure1,'renderer','Zbuffer')	% The patch command above set it to OpenGL, which is f... bugged

	% Get the index of anomalies that have names. We'll use them to plot those anomaly names
	ind = false(1,n_ages);
	for (i = 1:n_ages)
		if (~strcmp(age_txt(i),'a')),	ind(i) = true;		end
	end
	ages = age_start(ind);

	axes(handles.axes1)		% Put it back as current axes
	% Since the two axes are touching, axes1 would hide the XTickLabels of axes2.
	% So the trck is to plot what would have been axes2 XTickLabels as a text in axes1
	DX1 = diff(get(handles.axes1,'xlim'));
	y_lim = get(handles.axes1,'ylim');
	DX2 = age_end(end) - age_start(1);
	x_pos = (ages - x(1)) * DX1 / DX2;
	ha = 'Left';
	if (reverse_XDir)
		x_pos = (age_end(end) - x(1)) * DX1 / DX2 - x_pos;
		ha = 'Right';
	end
	text(x_pos,repmat(y_lim(2),numel(x_pos),1),age_txt(ind),'Parent',handles.axes1, ...
		'VerticalAlignment','top', 'HorizontalAlignment', ha)

% --------------------------------------------------------------------
function ac = autocorr(x)
%AUTOCORR Computes normalized auto-correlation of vector X.
	[x,nshift] = shiftdim(x);
	maxlag = size(x,1) - 1;
	x = x(:);   m = size(x,1);
	% Compute Autocorrelation via FFT
	X = fft(x,2^nextpow2(2*m-1));
	ac = ifft(abs(X).^2);

	ac = real(ac);				% We want only the real part
	% Move negative lags before positive
	ac = [ac(end-maxlag+1:end,:);ac(1:maxlag+1,:)];
	ac = ac./ac(maxlag+1);		% Normalize by ac[0]

	% If first vector is a row, return a row
	ac = shiftdim(ac,-nshift);

% --------------------------------------------------------------------
function [Pxx,w] = psd(xw,Fs)
%Power Spectral Density estimate via periodogram method.
	N = length(xw);     xw = xw(:);
	nfft  = max(256, 2^nextpow2(N));

	nx = size(xw,2);
	xw = [xw; zeros(nfft-N,1)];     % pad with zeros (I REALY don't like this)
	if (nx~=1),  xw = xw.';  end;    clear nx;

	% Compute the periodogram power spectrum [Power] estimate
	Sxx =(abs(fft(xw)).^2)./N; 

	% Generate the frequency vector in [rad/sample] at which Sxx was computed
	w = 2.*pi.*(0 : 1./nfft : 1-1./nfft);

	% Compute the Power/freq (PSD), the Power and the frequency at which it is computed
	w = w(:);

	% Generate the spectrum
	if rem(nfft,2),		select = 1:(nfft+1)/2;		% odd
	else				select = 1:nfft/2+1;		% even
	end
	Sxx_unscaled = Sxx(select);
	w = w(select);
	Sxx = [Sxx_unscaled(1); 2*Sxx_unscaled(2:end-1); Sxx_unscaled(end)];

	Pxx = Sxx./Fs;      % Scale by the sampling frequency to obtain the psd
	w = w.*Fs./(2.*pi); % Scale the frequency vector from rad/sample to Hz   

% -------------------------------------------------------------------------------------
function rd = dist_along_profile(x, y)
	xd = diff(x);		yd = diff(y);
	tmp = sqrt(xd.*xd + yd.*yd);
	rd = [0; cumsum(tmp(:))];

% -----------------------------------------------------------------------------
function figure1_KeyPressFcn(hObject, eventdata)
	if isequal(get(hObject,'CurrentKey'),'escape')
		handles = guidata(hObject);
		delete(handles.figure1)
	end

% -----------------------------------------------------------------------------
function figure1_CloseRequestFcn(hObject, eventdata)
	handles = guidata(hObject);
	if (~isempty(handles.polyFig))
		try,	delete(handles.polyFig),	end
	end
	delete(handles.figure1)

% --- Creates and returns a handle to the GUI figure. 
function ecran_LayoutFcn(h1)

set(h1,'Units','centimeters',...
'PaperUnits',get(0,'defaultfigurePaperUnits'),...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'KeyPressFcn',@figure1_KeyPressFcn,...
'CloseRequestFcn',@figure1_CloseRequestFcn,...
'MenuBar','none',...
'Name','ecran',...
'NumberTitle','off',...
'PaperPosition',[1 2.5 15.625 6.953],...
'PaperSize',[20.98404194812 29.67743169791],...
'PaperType',get(0,'defaultfigurePaperType'),...
'Position',[13.72 12.08 21.44 9.04],...
'RendererMode','manual',...
'Tag','figure1');

axes('Parent',h1,...
'CameraPosition',[0.5 0.5 9.16025403784439],...
'CameraPositionMode',get(0,'defaultaxesCameraPositionMode'),...
'NextPlot','Add',...
'Color',get(0,'defaultaxesColor'),...
'ColorOrder',get(0,'defaultaxesColorOrder'),...
'Position',[0.0493218249075216 0.14619883040935672 0.939580764488286 0.7923976608187134],...
'XColor',get(0,'defaultaxesXColor'),...
'YColor',get(0,'defaultaxesYColor'),...
'Tag','axes1');

axes('Parent',h1,...
'CameraPosition',[0.5 0.5 9.16025403784439],...
'CameraPositionMode',get(0,'defaultaxesCameraPositionMode'),...
'Color',get(0,'defaultaxesColor'),...
'ColorOrder',get(0,'defaultaxesColorOrder'),...
'Position',[0.0493218249075216 0.93859649122807 0.939580764488286 0.0584795321637427],...
'Tag','axes2',...
'Visible','off');

uicontrol('Parent',h1,...
'Units','normalized',...
'Callback',{@ecran_uicallback,h1,'checkbox_geog_Callback'},...
'Position',[0.0468557336621455 0.03216374269005848 0.2379778051787916 0.049707602339181284],...
'String','Geographical coordinates',...
'TooltipString',sprintf(['Check this if your data is in geographical coordinates.\n' ...
			'You will than be able to see and save the profile in km (or m) vs z.']),...
'Style','checkbox',...
'Tag','checkbox_geog');

uicontrol('Parent',h1,...
'Units','normalized',...
'BackgroundColor',[1 1 1],...
'Callback',{@ecran_uicallback,h1,'popup_selectPlot_Callback'},...
'Position',[0.3341553637484587 0.02046783625730994 0.3316892725030826 0.06432748538011696],...
'String','Distance along profile (data units)', ...
'Style','popupmenu',...
'Value',1,...
'TooltipString', 'Select different ways of seeing the profile', ...
'Tag','popup_selectPlot');

uicontrol('Parent',h1,...
'Units','normalized',...
'BackgroundColor',[1 1 1],...
'Callback',{@ecran_uicallback,h1,'popup_selectSave_Callback'},...
'Position',[0.7028360049321825 0.02046783625730994 0.2836004932182491 0.06432748538011696],...
'String',{  'Save Profile on disk'; 'distance Z (data units -> ascii)'; 'distance Z (data units -> binary)'; 'distance Z (km -> ascii)'; 'distance Z (km -> binary)'; 'distance Z (NM -> ascii)'; 'distance Z (NM -> binary)'; 'X Y Z (data units -> ascii)'; 'X Y Z (data units -> binary)' },...
'Style','popupmenu',...
'Value',1,...
'TooltipString', 'Choose how to save the profile', ...
'Tag','popup_selectSave');

h10 = uimenu('Parent',h1,'Label','File','Tag','menuFile');

uimenu('Parent',h10,...
'Callback',{@ecran_uicallback,h1,'FileOpen_Callback'},...
'Label','Open',...
'Tag','FileOpen');

uimenu('Parent',h10,...
'Callback',{@ecran_uicallback,h1,'FileSave_Callback'},...
'Label','Save',...
'Tag','FileSave');

uimenu('Parent',h10,'Callback','ecran','Label','New','Separator','on');

uimenu('Parent',h10,...
'Callback',{@ecran_uicallback,h1,'FileExport_Callback'},...
'Label','Export...',...
'Separator','on');

uimenu('Parent',h10, 'Callback','print -dsetup', 'Label','Print Setup', 'Separator','on');

uimenu('Parent',h10,...
'Callback',{@ecran_uicallback,h1,'FilePrint_Callback'},...
'Label','Print...',...
'Tag','FilePrint');

h17 = uimenu('Parent',h1, 'Label','Analysis');

uimenu('Parent',h17,...
'Callback',{@ecran_uicallback,h1,'AnalysisRemoveMean_Callback'},...
'Label','Remove Mean',...
'Tag','AnalysisRemoveMean');

uimenu('Parent',h17,...
'Callback',{@ecran_uicallback,h1,'AnalysisRemoveTrend_Callback'},...
'Label','Remove Trend',...
'Tag','AnalysisRemoveTrend');

uimenu('Parent',h17,...
'Callback',{@ecran_uicallback,h1,'AnalysisFitPoly_Callback'},...
'Label','Fit polynomial',...
'Tag','AnalysisFitPoly');

h20 = uimenu('Parent',h17, 'Label','FFT', 'Separator','on');

uimenu('Parent',h20,...
'Callback',{@ecran_uicallback,h1,'AnalysisFFT_AmpSpectrum_Callback'},...
'Label','Amplitude Spectrum',...
'Tag','AnalysisFFT_AmpSpectrum');

uimenu('Parent',h20,...
'Callback',{@ecran_uicallback,h1,'AnalysisFFT_PSD_Callback'},...
'Label','Power Spectrum Density',...
'Tag','AnalysisFFT_PSD');

uimenu('Parent',h17,...
'Callback',{@ecran_uicallback,h1,'AnalysisAutocorrelation_Callback'},...
'Label','Autocorrelation',...
'Tag','AnalysisAutocorrelation');

uimenu('Parent',h17,...
'Callback',{@ecran_uicallback,h1,'AnalysisSmoothSpline_Callback'},...
'Label','Smoothing Spline',...
'Separator','on',...
'Tag','AnalysisSmoothSpline');

uimenu('Parent',h17,...
'Callback',{@ecran_uicallback,h1,'Analysis1derivative_Callback'},...
'Label','1 st derivative',...
'Tag','Analysis1derivative');

uimenu('Parent',h17,...
'Callback',{@ecran_uicallback,h1,'Analysis2derivative_Callback'},...
'Label','2 nd derivative',...
'Tag','Analysis2derivative');

uicontrol('Parent',h1,...
'Units','normalized',...
'BackgroundColor',[1 1 1],...
'Callback',{@ecran_uicallback,h1,'edit_ageStart_Callback'},...
'Position',[0.319358816276202 0.0204678362573099 0.0579531442663379 0.0614035087719298],...
'String','0',...
'Style','edit',...
'Tooltip', sprintf(['Age at which we start ploting the bars. If older than "End age"\n' ...
	'the bar code is plotted reversed. That is, from older to younger ages']),...
'Tag','edit_ageStart',...
'Visible','off');

uicontrol('Parent',h1,...
'Units','normalized',...
'BackgroundColor',[1 1 1],...
'Callback',{@ecran_uicallback,h1,'edit_ageEnd_Callback'},...
'Position',[0.498150431565968 0.0233918128654971 0.0579531442663379 0.0614035087719298],...
'String','',...
'Style','edit',...
'Tooltip', 'Age at which we stop ploting the bars.',...
'Tag','edit_ageEnd',...
'Visible','off');

uicontrol('Parent',h1,...
'Units','normalized',...
'Callback',{@ecran_uicallback,h1,'push_magBar_Callback'},...
'FontName','Helvetica',...
'FontSize',9,...
'ListboxTop',0,...
'Position',[0.589395807644883 0.0204678362573099 0.129469790382244 0.0672514619883041],...
'String','Create Mag Bar',...
'TooltipString','Create a magnetic code bar on top of figure',...
'Tag','push_magBar',...
'Visible','off');

uicontrol('Parent',h1,...
'Units','normalized',...
'HorizontalAlignment','right',...
'ListboxTop',0,...
'Position',[0.250308261405672 0.0321637426900585 0.0678175092478422 0.043859649122807],...
'String','Start age',...
'Style','text',...
'Tag','text_ageStart',...
'Visible','off');

uicontrol('Parent',h1,...
'Units','normalized',...
'HorizontalAlignment','right',...
'ListboxTop',0,...
'Position',[0.429099876695438 0.0321637426900585 0.0678175092478422 0.043859649122807],...
'String','End age',...
'Style','text',...
'Tag','text_ageEnd',...
'Visible','off');

function ecran_uicallback(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
feval(callback_name,hObject,[],guidata(h1));

%============================================================================
function varargout = ecran_trend1d(varargin)
% Little help figure where to select the order of polynomial to fit
% The fit is done with the trend1d_m MEX because it allows robust fitting

	hObject = figure('Tag','figure1','Visible','off');
	ecran_trend1d_LayoutFcn(hObject);
	handles = guihandles(hObject);
	movegui(hObject,'center');

	handles.hCallingAx = varargin{1};
	handles.xy = varargin{2};
	if (size(handles.xy, 2) > 2)		% X,Y must be column vectors
		handles.xy = handles.xy';
	end

	handles.polyDeg = 1;

	guidata(hObject, handles);
	set(hObject,'Visible','on');
	if (nargout),	varargout{1} = hObject;		end

% --------------------------------------------------------------------
function edit_polDeg_Callback(hObject, eventdata, handles)
	xx = abs(fix(str2double(get(hObject,'String'))));
	if (isnan(xx))
		set(hObject,'String', handles.polyDeg)
		return
	end
	handles.polyDeg = xx;
	guidata(handles.figure1, handles)

% --------------------------------------------------------------------
function push_OK_Callback(hObject, eventdata, handles)
	opt_N = sprintf('-N%d', handles.polyDeg + 1);
	if (get(handles.check_robust, 'Val'))
		opt_N = [opt_N 'r'];
	end
	
	out = trend1d_m(handles.xy, '-Fxm', opt_N);
	h = line('XData', out(:,1), 'YData', out(:,2), 'Parent', handles.hCallingAx, 'Tag','fitted');
	
	% Compute the model parameters (trend1d_m only computes them in the linear case)
	p = polyfit(out(:,1), out(:,2), handles.polyDeg);
	
	% and put them on the line's uicontextmenu
	cmenuHand = uicontextmenu('Parent',get(handles.hCallingAx,'Parent'));
	set(h, 'UIContextMenu', cmenuHand);
 	uimenu(cmenuHand, 'Label', 'Poly Coefficients');
	uimenu(cmenuHand, 'Label', num2str(p));
	uimenu(cmenuHand, 'Label', 'Delete this line', 'Callback', 'delete(gco)', 'Sep', 'on');


% --- Creates and returns a handle to the GUI figure. 
function ecran_trend1d_LayoutFcn(h1)

set(h1,...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'MenuBar','none',...
'Name','Fit polynomial',...
'NumberTitle','off',...
'Position',[520 755 241 60],...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1');

uicontrol('Parent',h1,...
'Position',[161 38 75 15],...
'String','Robust Fit',...
'Style','checkbox',...
'TooltipString','Do a robust fit. See trend1d (GMT) manual to further details',...
'Tag','check_robust');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@ecran_trend1d_uicallback,h1,'edit_polDeg_Callback'},...
'Position',[111 34 30 21],...
'String','1',...
'Style','edit',...
'TooltipString','"1" means linear trend; "2" a quadratic model, and so on.',...
'Tag','edit_polDeg');

uicontrol('Parent',h1,...
'HorizontalAlignment','left',...
'Position',[10 37 100 15],...
'String','Polynomial degree',...
'Style','text');

uicontrol('Parent',h1,...
'Callback',{@ecran_trend1d_uicallback,h1,'push_OK_Callback'},...
'FontName','Helvetica',...
'FontSize',9,...
'Position',[164 6 66 21],...
'String','OK',...
'Tag','push_OK');

function ecran_trend1d_uicallback(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
feval(callback_name,hObject,[],guidata(h1));

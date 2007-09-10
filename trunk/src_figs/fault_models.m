function varargout = fault_models(varargin)
% M-File changed by desGUIDE 
 
	if (nargin)
		handMir = varargin{1};
	else
		return
	end

	[hObject, handles] = fault_models_LayoutFcn;
	movegui(hObject,'north')
	handles.handMir = handMir;
 
	load([handMir.path_data 'mirone_icons.mat'],'Mfopen_ico');
	set(handles.pushbutton_externalFile,'CData',Mfopen_ico)
	str = {'SUBFAULT FORMAT'; 'SRCMOD (.fsp)'};
	set(handles.listbox_readFilter,'String',str)

	% Add this figure handle to the carraças list
	plugedWin = getappdata(handMir.figure1,'dependentFigs');
	plugedWin = [plugedWin hObject];
	setappdata(handMir.figure1,'dependentFigs',plugedWin);

	handles.what_next = 'mansinha';

	guidata(hObject, handles);
	set(hObject,'Visible','on');
	if (nargout),	varargout{1} = hObject;		end

% -----------------------------------------------------------------------------
function push_externalFile_CB(hObject, handles)
	item = get(handles.listbox_readFilter,'Value');     % Get the reading filter number
	switch item
		case 1			% subfault
			str1 = {'*.dat;*.DAT', 'Data files (*.dat,*.DAT)';'*.*', 'All Files (*.*)'};
        case 2			% srcmod
			str1 = {'*.fsp;*.FSP', 'evTAG data files (*.fsp,*.FSP)';'*.*', 'All Files (*.*)'};
	end

	% Get file name
	[FileName,PathName] = put_or_get_file(handles.handMir,str1,'Select input file name','get');
	if isequal(FileName,0);     return,     end
	pause(0.05);

	cmap = hot(86);		cmap = cmap(84:-1:23,:);
	switch item
		case 1
			subfault(handles.handMir, handles, [PathName FileName], cmap)
		case 2
			evtag(handles.handMir, handles, [PathName FileName], cmap)
	end

% -----------------------------------------------------------------------------
function radio_Okada_CB(hObject, handles)
	if (get(hObject,'Val'))
		set([handles.radio_Mansinha handles.radio_nada],'Val',0)
		handles.what_next = 'okada';
	else
		set(hObject,'Val',1)
	end
	guidata(handles.figure1)

% -----------------------------------------------------------------------------
function radio_Mansinha_CB(hObject, handles)
	if (get(hObject,'Val'))
		set([handles.radio_Okada handles.radio_nada],'Val',0)
		handles.what_next = 'mansinha';
	else
		set(hObject,'Val',1)
	end
	guidata(handles.figure1)

% -----------------------------------------------------------------------------
function radio_nada_CB(hObject, handles)
	if (get(hObject,'Val'))
		set([handles.radio_Okada handles.radio_Mansinha],'Val',0)
		handles.what_next = 'nada';
	else
		set(hObject,'Val',1)
	end
	guidata(handles.figure1)

% ------------------------------------------------------------------------------------------------
function subfault(handles, localHandles, fname, cmap)
	% Read in a file in the subfault format
	% WARNING, this handles is the Mirone handles

	D2R = pi / 180;
	nCores = size(cmap,1) - 1;
	[numeric_data,multi_segs_str] = text_read(fname,NaN,5,'#');

	if (handles.no_file)		% I no_file, create a background
		w = min(numeric_data{1}(:,1));		e = max(numeric_data{1}(:,1));
		s = min(numeric_data{1}(:,2));		n = max(numeric_data{1}(:,2));
		mirone('FileNewBgFrame_CB',[],handles, [w e s n 1])   % Create a background
	end

	ind = strfind(multi_segs_str{1},'=');
	nSeg = str2double(multi_segs_str{1}(ind+1:end));

	%Fault_segment =   1 nx(Along-strike)=  30 Dx= 15.00km ny(downdip)=  15 Dy= 12.00km
	ind = strfind(multi_segs_str{2},'=');
	[t, r] = strtok(multi_segs_str{2}(ind(2)+1:end));		% t = 30;	r = [Dx= 15.00km ny(downdip)=  15 Dy= 12.00km]
	nx = str2double(t);
	ind = strfind(r,'=');
	[t, r] = strtok(r(ind(1)+1:end));		% t = 15.00km;	r = [ny(downdip)=  15 Dy= 12.00km]
	ind = strfind(t,'k');
	if (~isempty(ind)),		Dx = str2double(t(1:ind-1));
	else					Dx = str2double(t);
	end
	[t, r] = strtok(r);					% t = ny(downdip)=;		r = [15 Dy= 12.00km]
	[t, r] = strtok(r);					% t = 15;				r = [Dy= 12.00km]
	ny = str2double(t);
	ind = strfind(r,'=');
	t = strtok(r(ind(1)+1:end));		% t = 12.00km;
	ind = strfind(t,'k');
	if (~isempty(ind)),		Dy = str2double(t(1:ind-1));
	else					Dy = str2double(t);
	end
	
	rng = repmat((Dx / 2) / 6371 / D2R, nx, 1);
	numeric_data = numeric_data{2};			% First cell element contains the all patches rectangle
	numeric_data(:,4) = numeric_data(:,4) * 1e-2;	% they were in cm
	maxSlip = max(numeric_data(:,4));		minSlip = min(numeric_data(:,4));
	for (k = 1:ny)			% Loop over downdip patches (that is, rows)
        tmpx = numeric_data((k-1)*nx+1:k*nx, 2);		tmpy = numeric_data((k-1)*nx+1:k*nx, 1);
        [lat1,lon1] = circ_geo(tmpy, tmpx, rng, numeric_data((k-1)*nx+1:k*nx, 6)+180, 1);
        [lat2,lon2] = circ_geo(tmpy(end), tmpx(end), rng(end), numeric_data((k-1)*nx+1, 6), 1);
        lat = [lat1;lat2];      lon = [lon1;lon2];
        hLine(k) = line('XData',lon,'YData',lat,'Parent',handles.axes1,'Color','r','Tag','FaultTrace');
		
        %Lon. Lat. depth slip rake strike dip
        depth{k} = numeric_data((k-1)*nx+1:k*nx,3);		slip{k} = numeric_data((k-1)*nx+1:k*nx,4);
        rake{k}  = numeric_data((k-1)*nx+1:k*nx,5);		strike{k} = numeric_data((k-1)*nx+1:k*nx,6);
        dip{k}   = numeric_data((k-1)*nx+1:k*nx,7);
        for (i = 1:nx)		% Draw patches
            rng_p = Dy * cos(dip{k}(i)*D2R) / 6371 / D2R;
            [lat1,lon1] = circ_geo(lat(i),lon(i),rng_p,strike{k}(i)+90,1);
            [lat2,lon2] = circ_geo(lat(i+1),lon(i+1),rng_p,strike{k}(i)+90,1);
            x = [lon(i) lon(i+1) lon2 lon1];    y = [lat(i) lat(i+1) lat2 lat1];
			%z = -(depth{k}(i) * [1 1 1 1] + [0 0 [Dy Dy]*sin(dip{k}(1)*D2R)]);
			cor = cmap( round( (slip{k}(i) - minSlip) / (maxSlip - minSlip) * nCores + 1 ), :);
            hp(i) = patch('XData',x,'YData',y,'Parent',handles.axes1,'FaceColor',cor);
        end
        setappdata(hLine(k),'PatchHand',hp);
	end

	if (~handles.no_file)
		localHandles = guidata(localHandles.figure1);
		if (strcmp(localHandles.what_next,'okada'))
			deform_okada(handles,hLine,strike,depth,slip,dip,rake);
		elseif (strcmp(localHandles.what_next,'mansinha'))
			deform_mansinha(handles,hLine,strike,depth,slip,dip,rake);
		end
	end

% ------------------------------------------------------------------------------------------------
function evtag(handles, localHandles, fname, cmap)
	% Read in a file in the subfault format
	% WARNING, this handles is the Mirone handles

	D2R = pi / 180;
	nCores = size(cmap,1) - 1;
	[bin,n_column,multi_seg,n_headers] = guess_file(fname,5000,80);
	[numeric_data,multi_segs_str] = text_read(fname,NaN,n_headers,'%');

	% Search for a line like -> % Mech : STRK =  320         DIP =  11          RAKE =  92      Htop = 5.92 km
	k = 4;
	while (~strncmp(multi_segs_str{k},'% Mech :',8)),	k = k + 1;	end
	ind = strfind(multi_segs_str{k},'=');
	t = strtok(multi_segs_str{k}(ind(1)+1:end));		% t = 320;
	avgSTRK = str2double(t);
	t = strtok(multi_segs_str{k}(ind(2)+1:end));		% t = 11;
	avgDIP = str2double(t);
	t = strtok(multi_segs_str{k}(ind(3)+1:end));		% t = 92;
	avgRAKE = str2double(t);
	
	% Search for a line like -> % Invs :   Nx  =   30           Nz  =   15        Fmin = 0.01 Hz     Fmax = 0.5 Hz
	while (~strncmp(multi_segs_str{k},'% Invs :',8)),	k = k + 1;	end
	ind = strfind(multi_segs_str{k},'=');
	t = strtok(multi_segs_str{k}(ind(1)+1:end));		% t = 30;
	nx = str2double(t);
	t = strtok(multi_segs_str{k}(ind(2)+1:end));		% t = 15;
	nz = str2double(t);
	
	% Now, the two next lines are of this type
		% Invs :   Dx  =  15.00 km      Dz  = 12.00 km    
		% Invs :   Ntw =   1            Nsg =   1         (# of time-windows,# of fault segments)
	k = k + 1;
	ind = strfind(multi_segs_str{k},'=');
	t = strtok(multi_segs_str{k}(ind(1)+1:end));		% t = 15.00;
	Dx = str2double(t);
	t = strtok(multi_segs_str{k}(ind(2)+1:end));		% t = 12.00;
	Dy = str2double(t);
	k = k + 1;
	ind = strfind(multi_segs_str{k},'=');
	t = strtok(multi_segs_str{k}(ind(2)+1:end));		% t = 1;	We jump Ntw for now
	nSeg = str2double(t);

	% I no_file, create a background
	if (handles.no_file)
		w = 360;	e = -180;	s = 90;		n = -90;
		for (k = 1:nSeg)
			w = min(min(numeric_data{k}(:,2)), w);		e = max(max(numeric_data{k}(:,2)), e);
			s = min(min(numeric_data{k}(:,1)), s);		n = max(max(numeric_data{k}(:,1)), n);
		end
		offx = Dx / 6371 / D2R;		offy = Dy / 6371 / D2R;
		w = w - offx;		e = e + offx;
		s = s - offy;		n = n + offy;
		mirone('FileNewBgFrame_CB',[],handles, [w e s n 1])   % Create a background
	end

	% If we have a multi segment file, there is a lot more to fish in this squizophrenic format
	if (nSeg > 1)
		% Start a bottom and go up until we find the info for all nSeg segments
		Nsbfs = zeros(1,nSeg);		avgSTRK = zeros(1,nSeg);
		avgDIP = zeros(1,nSeg);		nx = zeros(1,nSeg);
		n_seg = nSeg;		k = numel(multi_segs_str);
		while (n_seg)
			while (~strncmp(multi_segs_str{k},'%   Nsbfs',9)),	k = k - 1;	end
			ind = strfind(multi_segs_str{k},'=');
			t = strtok(multi_segs_str{k}(ind(1)+1:end));
			Nsbfs(n_seg) = str2double(t);
			while (~strncmp(multi_segs_str{k},'% SEGMENT',9)),	k = k - 1;	end
			ind = strfind(multi_segs_str{k},'=');
			t = strtok(multi_segs_str{k}(ind(1)+1:end));
			avgSTRK(n_seg) = str2double(t);
			t = strtok(multi_segs_str{k}(ind(2)+1:end));
			avgDIP(n_seg) = str2double(t);
			nx(n_seg) = Nsbfs(n_seg) / nz;
			n_seg = n_seg - 1;
		end
	end

	% Let us see if we have the rake the 7th column
	if ( strfind(multi_segs_str{end-1}(55:end), 'RAKE') )
		rake_in_7 = true;
	else
		rake_in_7 = false;
	end
	
	if (nSeg > 1 && ~rake_in_7)
		warndlg('I don''t know where those guys have put the RAKE. Using the average (bad) one.','Warning')
	end
	
	for (s = 1:nSeg)
		data = numeric_data{s};
		strk = repmat(avgSTRK(s),nx(s),1);
		strike = repmat({strk},nz,1);
		dip  = repmat({repmat(avgDIP(s),nx(s),1)},nz,1);
		rng = repmat((Dx / 2) / 6371 / D2R, nx(s), 1);
		if (~rake_in_7),		rake  = repmat({repmat(avgRAKE,nx(s),1)},nz,1);		end
		maxSlip = max(data(:,6));	minSlip = min(data(:,6));
		for (k = 1:nz)			% Loop over downdip patches (that is, rows)
            tmpx = data((k-1)*nx(s)+1:k*nx(s), 2);		tmpy = data((k-1)*nx(s)+1:k*nx(s), 1);
            [lat1,lon1] = circ_geo(tmpy, tmpx, rng, strk+180, 1);
            [lat2,lon2] = circ_geo(tmpy(end), tmpx(end), rng(end), strk, 1);
            lat = [lat1;lat2];      lon = [lon1;lon2];
            hLine(k) = line('XData',lon,'YData',lat,'Parent',handles.axes1,'Color','r','Tag','FaultTrace');
			
            depth{k} = data((k-1)*nx(s)+1:k*nx(s),5);		slip{k} = data((k-1)*nx(s)+1:k*nx(s),6);
			if (rake_in_7),		rake{k}  = data((k-1)*nx(s)+1:k*nx(s),7);		end
			
            for (i = 1:nx(s))		% Draw patches
                rng_p = Dy * cos(dip{k}(i)*D2R) / 6371 / D2R;
                [lat1,lon1] = circ_geo(lat(i),lon(i),rng_p,strike{k}(i)+90,1);
                [lat2,lon2] = circ_geo(lat(i+1),lon(i+1),rng_p,strike{k}(i)+90,1);
                x = [lon(i) lon(i+1) lon2 lon1];    y = [lat(i) lat(i+1) lat2 lat1];
				%z = (depth{k}(i) * [1 1 1 1] + [0 0 [Dy Dy]*sin(dip{k}(1)*D2R)]);
				cor = cmap( round( (slip{k}(i) - minSlip) / (maxSlip - minSlip) * nCores + 1 ), :);
                hp(i) = patch('XData',x,'YData',y,'Parent',handles.axes1,'FaceColor',cor);
            end
            setappdata(hLine(k),'PatchHand',hp);
		end
	end

	if (~handles.no_file)
		localHandles = guidata(localHandles.figure1);
		if (strcmp(localHandles.what_next,'okada'))
			deform_okada(handles,hLine,strike,depth,slip,dip,rake);
		elseif (strcmp(localHandles.what_next,'mansinha'))
			deform_mansinha(handles,hLine,strike,depth,slip,dip,rake);
		end
	end

% -----------------------------------------------------------------------------
function figure1_KeyPressFcn(hObject, eventdata)
	if isequal(get(hObject,'CurrentKey'),'escape')
        delete(hObject);
	end

% --- Creates and returns a handle to the GUI figure. 
function [h1, handles] = fault_models_LayoutFcn()

h1 = figure('Tag','figure1','Visible','off', ...
'PaperUnits',get(0,'defaultfigurePaperUnits'),...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'KeyPressFcn',@figure1_KeyPressFcn,...
'MenuBar','none',...
'Name','Fault models',...
'NumberTitle','off',...
'PaperPosition',get(0,'defaultfigurePaperPosition'),...
'PaperSize',[20.98404194812 29.67743169791],...
'PaperType',get(0,'defaultfigurePaperType'),...
'Position',[520 600 240 140],...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Position',[10 73 201 61],...
'Style','listbox',...
'Value',1,...
'Tag','listbox_readFilter');

uicontrol('Parent',h1, 'Position',[209 93 23 23],...
'Callback',{@fault_models_uicallback,h1,'push_externalFile_CB'},...
'FontWeight','bold',...
'TooltipString','Browse for wanted file',...
'Tag','pushbutton_externalFile');

uicontrol('Parent',h1, 'Position',[10 49 160 15],...
'Callback',{@fault_models_uicallback,h1,'radio_Okada_CB'},...
'Style','radiobutton',...
'String', 'At the end run ''Okada'' tool', ...
'TooltipString','When finished reading call the Okada window',...
'Tag','radio_Okada');

uicontrol('Parent',h1, 'Position',[10 29 160 15],...
'Callback',{@fault_models_uicallback,h1,'radio_Mansinha_CB'},...
'Style','radiobutton',...
'Value',1,...
'String', 'At the end run ''Mansinha'' tool', ...
'TooltipString','When finished reading call the Mansinha window',...
'Tag','radio_Mansinha');

uicontrol('Parent',h1, 'Position',[10 9 160 15],...
'Callback',{@fault_models_uicallback,h1,'radio_nada_CB'},...
'Style','radiobutton',...
'String', 'Do nothing at the end', ...
'TooltipString','When finished reading ... do nothing else',...
'Tag','radio_nada');

% Create handles
handles = guihandles(h1);

function fault_models_uicallback(hObject, eventdata, h1, callback_name)
	% This function is executed by the callback and than the handles is allways updated.
	feval(callback_name,hObject,guidata(h1));

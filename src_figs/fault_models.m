function varargout = fault_models(varargin)
% Helper window to import fault models

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

% $Id$

	if (nargin)
		handMir = varargin{1};
	end

	[hObject, handles] = fault_models_LayoutFcn;
	handles.handMir = varargin{1};
	move2side(handMir.figure1, hObject)

	load([handMir.path_data 'mirone_icons.mat'],'Mfopen_ico');
	set(handles.push_externalFile,'CData',Mfopen_ico)
	str = {'SUBFAULT FORMAT'; 'SRCMOD (.fsp)'};
	set(handles.listbox_readFilter,'String',str)
	handles.last_dir = handles.handMir.last_dir;
	handles.work_dir = handles.handMir.work_dir;

	% Add this figure handle to the carra?as list
	plugedWin = getappdata(handMir.figure1,'dependentFigs');
	plugedWin = [plugedWin hObject];
	setappdata(handMir.figure1,'dependentFigs',plugedWin);

	handles.what_next = 'none';

	guidata(hObject, handles);
	set(hObject,'Vis','on');
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
	[FileName,PathName] = put_or_get_file(handles,str1,'Select input file name','get');
	if isequal(FileName,0);     return,     end
	pause(0.05);

	[pato, fname, EXT] = fileparts(FileName);
	cmap = hot(86);		cmap = cmap(84:-1:23,:);
	switch item
		case 1
			if (strcmpi(EXT,'.fsp'))
				errordlg('This file is in the .fsp format and not in the subfault''s one.','ERROR')
				return
			end
			subfault(handles.handMir, handles, [PathName FileName], cmap)
		case 2
			if (strcmpi(EXT,'.slp'))
				errordlg('This file is in the .slp format and not in the .fsp','ERROR')
				return
			end
			evtag(handles.handMir, handles, [PathName FileName], cmap)
	end

% -----------------------------------------------------------------------------
function radio_Okada_CB(hObject, handles)
	if (get(hObject,'Val'))
		set([handles.radio_Mansinha handles.radio_nada],'Val',0)
		handles.what_next = 'okada';
		guidata(handles.figure1, handles)
	else
		set(hObject,'Val',1)
	end

% -----------------------------------------------------------------------------
function radio_Mansinha_CB(hObject, handles)
	if (get(hObject,'Val'))
		set([handles.radio_Okada handles.radio_nada],'Val',0)
		handles.what_next = 'mansinha';
		guidata(handles.figure1, handles)
	else
		set(hObject,'Val',1)
	end

% -----------------------------------------------------------------------------
function radio_nada_CB(hObject, handles)
	if (get(hObject,'Val'))
		set([handles.radio_Okada handles.radio_Mansinha],'Val',0)
		handles.what_next = 'nada';
		guidata(handles.figure1, handles)
	else
		set(hObject,'Val',1)
	end

% ------------------------------------------------------------------------------------------------
function subfault(handles, localHandles, fname, cmap)
% Read in a file in the subfault format
% WARNING, the handles is the Mirone handles

	D2R = pi / 180;
	nCores = size(cmap,1) - 1;
	[numeric_data,multi_segs_str] = text_read(fname,NaN,0,'#');		% We force 0 header lines because the'r all multi-segs

	w = min(numeric_data{1}(:,1));		e = max(numeric_data{1}(:,1));
	s = min(numeric_data{1}(:,2));		n = max(numeric_data{1}(:,2));
	if (handles.no_file)		% If no_file, create a background
		mirone('FileNewBgFrame_CB', handles, [w e s n 1]);			% Create a background
		hFig = handles.figure1;
		figName = get(hFig,'Name');
		ind = strfind(figName,' @');
		if (~isempty(ind)),		figName = [fname figName(ind(1):end)];
		else					figName = fname;
		end
		set(hFig,'Name',figName)
	elseif (all(~insideRect(handles.head(1:4), [[w; e] [s; n]])))
		errordlg('Time to go to the oculist? The fault model is entirely out of the map region.','ERROR')
		return
	end

	%ind = strfind(multi_segs_str{1},'=');
	%nSeg = str2double(multi_segs_str{1}(ind+1:end));

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
	
	% Guess if the grid is in km or meters. If the guess fails so will the 3D flederization
	fact = 1;
	if (~handles.no_file && handles.validGrid)
		rng = diff(handles.head(5:6));
		if (rng > 12),		fact = 1e3;		end		% Grid is meters (but depth is in km)
	end

	rng = repmat((Dx / 2) / 6371 / D2R, nx, 1);
	numeric_data = numeric_data{2};			% First cell element contains the all patches rectangle
	numeric_data(:,4) = numeric_data(:,4) * 1e-2;	% they were in cm
	maxSlip = max(numeric_data(:,4));		minSlip = min(numeric_data(:,4));
	% Pre-allocations
	depth = cell(1,ny);		slip = cell(1,ny);		rake = cell(1,ny);
	strike = cell(1,ny);	dip = cell(1,ny);		hp = zeros(1,nx);		hLine = zeros(1,ny);
	
	for (k = 1:ny)			% Loop over downdip patches (that is, rows)
        tmpx = numeric_data((k-1)*nx+1:k*nx, 2);		tmpy = numeric_data((k-1)*nx+1:k*nx, 1);
        [lat1,lon1] = circ_geo(tmpy, tmpx, rng, numeric_data((k-1)*nx+1:k*nx, 6)+180, 1);
        [lat2,lon2] = circ_geo(tmpy(end), tmpx(end), rng(end), numeric_data((k-1)*nx+1, 6), 1);
        lat = [lat1;lat2];      lon = [lon1;lon2];
        hLine(k) = line('XData',lon,'YData',lat,'Parent',handles.axes1,'Color','k','Tag','FaultTrace', 'HitTest','off', ...
			'Vis','off','HandleVisibility','off');		% So that they are not fished in the flederization process
		
        %Lon. Lat. depth slip rake strike dip
        depth{k} = numeric_data((k-1)*nx+1:k*nx,3);		slip{k} = numeric_data((k-1)*nx+1:k*nx,4);
        rake{k}  = numeric_data((k-1)*nx+1:k*nx,5);		strike{k} = numeric_data((k-1)*nx+1:k*nx,6);
        dip{k}   = numeric_data((k-1)*nx+1:k*nx,7);
		z = -(depth{k}(1) * [1 1 1 1 1] + [0 0 [Dy Dy]*sin(dip{k}(1)*D2R) 0]) * fact;
        for (i = 1:nx)		% Draw patches
            rng_p = Dy * cos(dip{k}(i)*D2R) / 6371 / D2R;
            [lat1,lon1] = circ_geo(lat(i),lon(i),rng_p,strike{k}(i)+90,1);
            [lat2,lon2] = circ_geo(lat(i+1),lon(i+1),rng_p,strike{k}(i)+90,1);
            x = [lon(i) lon(i+1) lon2 lon1 lon(i)];    y = [lat(i) lat(i+1) lat2 lat1 lat(i)];
			cor = cmap( round( (slip{k}(i) - minSlip) / (maxSlip - minSlip) * nCores + 1 ), :);
            hp(i) = patch('XData',x,'YData',y,'Parent',handles.axes1,'FaceColor',cor,'Tag','SlipPatch');
			set(hp(i),'UserData',z)		% So that we can Flederize it in 3D
			cmenuHand = uicontextmenu('Parent',handles.figure1);
			set(hp(i), 'UIContextMenu', cmenuHand);
			uimenu(cmenuHand, 'Label', 'Okada', 'Callback', {@calldeform,handles,'okada'});
			uimenu(cmenuHand, 'Label', 'Mansinha', 'Callback', {@calldeform,handles,'mansinha'});
			uimenu(cmenuHand, 'Label', 'Delete', 'Sep', 'on', 'Callback', {@delAll,handles});    
        end
        setappdata(hLine(k),'PatchHand',hp);
	end

	% Save those to accessible by calldeform()
	width = repmat({repmat(Dy, nx, 1)}, 1, ny);		% Width is constant when single segment
	setappdata(handles.axes1,'SlipVars',{hLine, depth, width, strike,slip,dip,rake})

	if (~handles.no_file)
		localHandles = guidata(localHandles.figure1);
		if (strcmp(localHandles.what_next,'okada'))
			deform_okada(handles,hLine,depth,width,strike,slip,dip,rake);
		elseif (strcmp(localHandles.what_next,'mansinha'))
			deform_mansinha(handles,hLine,depth,width,strike,slip,dip,rake);
		end
	end

% ------------------------------------------------------------------------------------------------
function evtag(handles, localHandles, fname, cmap)
% Read in a file in the subfault format
% WARNING, the handles is the Mirone handles
% This function differs from 'subfault()' mainly in what respects the decoding of the
% squizophrenic header of .fsp files. The other difference is that multi-segment faults
% are alowed (it's than that things get real squizo)
%
% A basic assumption that all segments have the same 'nz'
%
% With the complicated multi-segment format each model parameter (e.g. slip)
% is a cell array {nSeg, nz} where nSeg is the number of segments as deffined
% by the format specification (.fsp here) and nz is the number of descretization
% stripes along strike and downdip (that is, rows in the fault plane).
% Then each cell element is a NX column vector with the parameter's value, where
% NX is the segment number of patches along strike.

	D2R = pi / 180;
	nCores = size(cmap,1) - 1;
	[bin,n_column,multi_seg,n_headers] = guess_file(fname,5000,80);
	if (n_headers < 20)
		errordlg('This file is not in the SRCMOD .fsp format.','ERROR')
		return
	end
	[numeric_data,multi_segs_str] = text_read(fname,NaN,0,'%');		% Used to send n_headers but that showed up as a bug!

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

	% Get min/max for eventual background creation or inside-region check
	w = 360;	e = -180;	s = 90;		n = -90;
	for (k = 1:nSeg)
		w = min(min(numeric_data{k}(:,2)), w);		e = max(max(numeric_data{k}(:,2)), e);
		s = min(min(numeric_data{k}(:,1)), s);		n = max(max(numeric_data{k}(:,1)), n);
	end
	
	% I no_file, create a background
	if (handles.no_file)
		offx = Dx / 6371 / D2R;		offy = Dy / 6371 / D2R;
		w = w - offx;		e = e + offx;
		s = s - offy;		n = n + offy;
		mirone('FileNewBgFrame_CB', handles, [w e s n 1])	% Create a background
		hFig = handles.figure1;
		figName = get(hFig,'Name');
		ind = strfind(figName,' @');
		if (~isempty(ind)),		figName = [fname figName(ind(1):end)];
		else					figName = fname;
		end
		set(hFig,'Name',figName)
	elseif (all(~insideRect(handles.head(1:4), [[w; e] [s; n]])))
		errordlg('Time to go to the oculist? The fault model is entirely out of the map region.','ERROR')
		return
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
		warndlg('Where is the RAKE of this file?. Using the average (bad) one.','Warning')
	end
	
	% Guess if the grid is in km or meters. If the guess fails so will the 3D flederization
	fact = 1;
	if (~handles.no_file && handles.validGrid)
		rng = diff(handles.head(5:6));
		if (rng > 12),		fact = 1e3;		end		% Grid is meters (but depth is in km)
	end

	% Pre-allocations
	depth  = cell(nSeg,nz);		slip = cell(nSeg,nz);		rake = cell(nSeg,nz);
	strike = cell(nSeg,nz);		dip = cell(nSeg,nz);		hLine = zeros(nSeg,nz);
	width  = cell(nSeg,nz);

	for (s = 1:nSeg)
		data = numeric_data{s};
		strk = repmat(avgSTRK(s),nx(s),1);
		strike(s,:) = repmat({strk},nz,1);
		dip(s,:)  = repmat({repmat(avgDIP(s),nx(s),1)},nz,1);
		rng = repmat((Dx / 2) / 6371 / D2R, nx(s), 1);
		if (~rake_in_7),		rake(s,:)  = repmat({repmat(avgRAKE,nx(s),1)},nz,1);		end
		maxSlip = max(data(:,6));	minSlip = min(data(:,6));
		
		for (k = 1:nz)			% Loop over downdip patches (that is, rows)
            tmpx = data((k-1)*nx(s)+1:k*nx(s), 2);		tmpy = data((k-1)*nx(s)+1:k*nx(s), 1);
            [lat1,lon1] = circ_geo(tmpy, tmpx, rng, strike{s,k}(1)+180, 1);
            [lat2,lon2] = circ_geo(tmpy(end), tmpx(end), rng(end), strike{s,k}(1), 1);
            lat = [lat1;lat2];      lon = [lon1;lon2];
	        hLine(s,k) = line('XData',lon,'YData',lat,'Parent',handles.axes1,'Color','k','Tag','FaultTrace', 'HitTest','off', ...
				'Vis','off','HandleVisibility','off');		% So that they are not fished in the flederization processe
			
            depth{s,k} = data((k-1)*nx(s)+1:k*nx(s),5);		slip{s,k} = data((k-1)*nx(s)+1:k*nx(s),6);
			width{s,k} = repmat(Dy, nx(s), 1);
			if (rake_in_7),		rake{s,k}  = data((k-1)*nx(s)+1:k*nx(s),7);		end
			hp = zeros(1,nx(s));		% Pre-allocation and resetting from a previous (now old) size
			
            for (i = 1:nx(s))		% Draw patches
                rng_p = Dy * cos(dip{s,k}(i)*D2R) / 6371 / D2R;
                [lat1,lon1] = circ_geo(lat(i),lon(i),rng_p,strike{s,k}(i)+90,1);
                [lat2,lon2] = circ_geo(lat(i+1),lon(i+1),rng_p,strike{s,k}(i)+90,1);
	            x = [lon(i) lon(i+1) lon2 lon1 lon(i)];    y = [lat(i) lat(i+1) lat2 lat1 lat(i)];
				z = -(depth{s,k}(i) * [1 1 1 1 1] + [0 0 [Dy Dy]*sin(dip{s,k}(1)*D2R) 0]) * fact;
				cor = cmap( round( (slip{s,k}(i) - minSlip) / (maxSlip - minSlip) * nCores + 1 ), :);
                hp(i) = patch('XData',x,'YData',y,'Parent',handles.axes1,'FaceColor',cor,'Tag','SlipPatch');
				set(hp(i),'UserData',z)		% So that we can Flederize it in 3D
				cmenuHand = uicontextmenu('Parent',handles.figure1);
				set(hp(i), 'UIContextMenu', cmenuHand);
				uimenu(cmenuHand, 'Label', 'Okada', 'Callback', {@calldeform,handles,'okada'});
				uimenu(cmenuHand, 'Label', 'Mansinha', 'Callback', {@calldeform,handles,'mansinha'});
				uimenu(cmenuHand, 'Label', 'Delete', 'Sep', 'on', 'Callback', {@delAll,handles});    
            end
            setappdata(hLine(s,k),'PatchHand',hp);
		end
	end

	% Save those to accessible by calldeform()
	setappdata(handles.axes1,'SlipVars',{hLine, depth, width, strike,slip,dip,rake})

	if (~handles.no_file)
		localHandles = guidata(localHandles.figure1);
		if (strcmp(localHandles.what_next,'okada'))
			deform_okada(handles,hLine,depth,width,strike,slip,dip,rake);
		elseif (strcmp(localHandles.what_next,'mansinha'))
			deform_mansinha(handles,hLine,depth,width,strike,slip,dip,rake);
		end
	end

% -----------------------------------------------------------------------------------------
function delAll(obj,event,handles)
	delete(findobj(handles.axes1,'Tag','SlipPatch'))
	vars = getappdata(handles.axes1,'SlipVars');
	delete(vars{1});		% Delete the hiden lines
	rmappdata(handles.axes1,'SlipVars')

% -----------------------------------------------------------------------------------------
function calldeform(obj,event, handles, opt)
% Remember, this handles is the one from Mirone
	vars = getappdata(handles.axes1,'SlipVars');
	handles = guidata(handles.figure1);
	if (strcmp(opt,'okada'))
		deform_okada(handles,vars{:});
	else
		deform_mansinha(handles,vars{:});
	end

% --------------------------------------------------------------------
function res = insideRect(rect,pt)
% Check which elements of the  [x y] (Mx2) PT array are inside the rectangle RECT
% RECT = [x_min x_max y_min y_max]
% RES is a logical column vector with length = size(PT,1)
% NO ERROR TESTING
    res = ( pt(:,1) >= rect(1) & pt(:,1) <= rect(2) & pt(:,2) >= rect(3) & pt(:,2) <= rect(4) );

% -----------------------------------------------------------------------------
function figure1_KeyPressFcn(hObject, eventdata)
	if isequal(get(hObject,'CurrentKey'),'escape')
        delete(hObject);
	end

% --- Creates and returns a handle to the GUI figure. 
function [h1, handles] = fault_models_LayoutFcn()

h1 = figure('Visible','off', ...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'KeyPressFcn',@figure1_KeyPressFcn,...
'MenuBar','none',...
'Name','Fault models',...
'NumberTitle','off',...
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
'Call',@fault_models_uiCB,...
'FontWeight','bold',...
'Tooltip','Browse for wanted file',...
'Tag','push_externalFile');

uicontrol('Parent',h1, 'Position',[10 49 170 15],...
'Call',@fault_models_uiCB,...
'Style','radiobutton',...
'String', 'At the end run ''Okada'' tool', ...
'Tooltip','When finished reading call the Okada window',...
'Tag','radio_Okada');

uicontrol('Parent',h1, 'Position',[10 29 190 15],...
'Call',@fault_models_uiCB,...
'Style','radiobutton',...
'String', 'At the end run ''Mansinha'' tool', ...
'Tooltip','When finished reading call the Mansinha window',...
'Tag','radio_Mansinha');

uicontrol('Parent',h1, 'Position',[10 9 160 15],...
'Call',@fault_models_uiCB,...
'Style','radiobutton',...
'Value',1,...
'String', 'Do nothing at the end', ...
'Tooltip','When finished reading ... do nothing else',...
'Tag','radio_nada');

% Create handles
handles = guihandles(h1);

function fault_models_uiCB(hObject, eventdata)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));

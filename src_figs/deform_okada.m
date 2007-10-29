function varargout = deform_okada(varargin)
% Compute Elastic deformations (M-File changed by desGUIDE)

%	Copyright (c) 2004-2006 by J. Luis
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

	if isempty(varargin)
		errordlg('DEFORM OKADA: Wrong number of input args','Error');    return
	end
    
	hObject = figure('Tag','figure1','Visible','off');
	deform_okada_LayoutFcn(hObject);
	handles = guihandles(hObject);
	movegui(hObject,'east');

	handles.input_locations = [];   % May contain ground points positions
	D2R = pi / 180;

	if (length(varargin) >= 8)      % Fault patch collection
		handles = set_all_faults(handles,varargin{:});
		handles.fault_in = 1;
	else                            % "Normal" case
		handles.h_fault = varargin{2};      % Handles to the fault lines (each may have more than one segment)
		handles.FaultStrike = varargin{3};
		handles.fault_in = 0;
	end

	handMir = varargin{1};
	handles.geog = handMir.geog;
	head = handMir.head;
	handles.head = head;
	handles.h_calling_fig = handMir.figure1;     % Handles to the calling figure
	handles.n_faults = length(handles.h_fault);

	if (handles.n_faults > 1)
		s_format = ['%.' num2str(fix(log10(handles.n_faults))+1) 'd'];
		S = cell(handles.n_faults,1);
		for (i=1:handles.n_faults),     S{i} = ['Fault ' sprintf(s_format,i)];   end
		set(handles.popup_fault,'String',S)
		set(handles.h_fault(1),'LineStyle','--');   % set the top fault one with a dashed line type
		refresh(handMir.figure1);             % otherwise, ML BUG
	else
		set(handles.popup_fault,'Visible','off')
		delete(handles.fault_number)
	end

	fault_x = get(handles.h_fault,'XData');     fault_y = get(handles.h_fault,'YData');
	if (handles.n_faults > 1)
        nvert = zeros(1,handles.n_faults);
		for (k=1:handles.n_faults), nvert(k) = size(fault_x{k},2) - 1;  end
	else
		nvert = size(fault_x,2) - 1;
	end

	if (any(nvert > 1))
		set(handles.popup_segment,'Visible','on');
		% Even if we have more than one fault, the segments popup will start with only the first fault's segments
		s_format = ['%.' num2str(fix(log10(max(nvert)))+1) 'd'];
		S = cell(nvert(1),1);
		for (i=1:nvert(1)),     S{i} = ['Segment ' sprintf(s_format,i)];   end
		set(handles.popup_segment,'String',S)
	else
		set(handles.popup_segment,'Visible','off')
		delete(handles.txtFaultSeg)    % Otherwise it would reborn in Pro look
	end

	% Try to guess if we are dealing with other (m or km) than geogs
	handles.is_meters = 0;     handles.is_km = 0;      handles.um_milhao = 1e6;
	if (~handles.geog)      % Try to guess if user units are km or meters
		dx = head(2) - head(1);   dy = head(4) - head(3);
		len = sqrt(dx.*dx + dy.*dy);         % Distance in user unites
		if (len > 1e5 || head(8) >= 10)      % If grid's diagonal > 1e5 || Dx >= 10 consider we have meters
			handles.is_meters = 1;         handles.um_milhao = 1e3;
			set(handles.popup_GridCoords,'Value',2)
		else
			handles.is_km = 1;
			set(handles.popup_GridCoords,'Value',3)
		end
	end

	handles.fault_x = fault_x;
	handles.fault_y = fault_y;
	handles.nvert = nvert;
	handles.hide_planes(1:handles.n_faults) = 0;
	handles.dms_xinc = 0;           handles.dms_yinc = 0;
	handles.txt_Mw_pos = get(handles.h_txt_Mw,'Position');
	handles.Mw(1:handles.n_faults) = 0;
	handles.FaultLength = LineLength(handles.h_fault,handles.geog);
	handles.one_or_zero = ~head(7);
	handles.x_min_or = head(1);			handles.x_max_or = head(2);
	handles.y_min_or = head(3);			handles.y_max_or = head(4);

	if (~handles.fault_in)					% "NORMAL" case (not a fault-patch collection)
		% Make them all cell arrays to simplify logic
		if (~iscell(handles.FaultLength)),  handles.FaultLength = {handles.FaultLength};   end
		if (~iscell(handles.FaultStrike)),  handles.FaultStrike = {handles.FaultStrike};   end
		if (~iscell(handles.fault_x)),      handles.fault_x = {handles.fault_x};    handles.fault_y = {handles.fault_y};   end
		handles.DislocStrike = handles.FaultStrike;

		for (k=1:handles.n_faults)
			handles.FaultDip{k}(1:nvert(k)) = 25;       handles.FaultWidth{k}(1:nvert(k)) = NaN;
			handles.FaultDepth{k}(1:nvert(k)) = NaN;	handles.FaultTopDepth{k}(1:nvert(k)) = 0;
			handles.DislocSlip{k}(1:nvert(k)) = 1;		handles.DislocRake{k}(1:nvert(k)) = 90;
			handles.ux{k}(1:nvert(k)) = 0;              handles.uy{k}(1:nvert(k)) = 1;
			handles.uz{k}(1:nvert(k)) = 0;
		end

		z1 = num2str(handles.FaultLength{1}(1));    z2 = sprintf('%.1f',handles.FaultStrike{1}(1));
		z3 = sprintf('%.1f',handles.FaultDip{1}(1));
		set(handles.edit_FaultLength,'String',z1,'Enable','off')
		set(handles.edit_FaultStrike,'String',z2,'Enable','off')
		set(handles.edit_FaultDip,'String',z3)
		set(handles.edit_DislocStrike,'String',z2)
		set(handles.edit_DislocSlip,'String','1')
		set(handles.edit_DislocRake,'String','90')
		handles.DislocSlip{1}(1) = 1;
		handles.DislocRake{1}(1) = 90;

		% Set a default unit dislocation as a thrust motion
		set(handles.edit_ux,'String',0)
		set(handles.edit_uy,'String',1)
		set(handles.edit_uz,'String',0)
	
		set(handles.edit_FaultTopDepth,'String','0')    % Default the top depth fault to zero
		% If we have one fault only, provide a default Width value
		if (handles.n_faults == 1)
			faultWidth = handles.FaultLength{1}(1) / 4;
			if (handles.is_meters),		faultWidth = round(faultWidth * 1e-3);     end
			handles = edit_FaultWidth_Callback([], faultWidth, handles);    % Compute the rest
	        set(handles.edit_FaultWidth,'String',num2str(faultWidth));
		end
	else
		set(handles.edit_FaultLength,'Enable','off')
		set(handles.edit_FaultStrike,'Enable','off')
        set(handles.edit_FaultLength,'String',num2str(handles.FaultLength{1}(1)));
		% Compute Mag for each fault
		totalM0 = 0;
		for (k = handles.n_faults:-1:1)
			[handles,mag,MO] = compMag(handles, k);
			totalM0 = totalM0 + MO;
		end
		mag = 2/3*(log10(totalM0) - 9.1);
		uicontrol('Parent',hObject,'Enable','inactive','FontSize',10,...
		'HorizontalAlignment','left','Position',[400 190 100 16],...
		'String',['Tot Mw = ' sprintf('%.1f',mag)],'Style','text');

	end
	
	% Set a default view vector (for interferograms s must be different)
	sx = cos((90-handles.FaultStrike{1})*D2R);   handles.sx = sx;
	sy = sin((90-handles.FaultStrike{1})*D2R);   handles.sy = sy;       handles.sz(1:nvert) = 0;
	set(handles.edit_sx,'String',num2str(sx(1)))
	set(handles.edit_sy,'String',num2str(sy(1)))
	set(handles.edit_sz,'String','0')

	%-----------
	% Fill in the grid limits boxes (in case user wants to compute a grid)
	nDigit = round( log10(abs(max(head(1:4)))) );		% Number of digits of the integer part
	frmt = sprintf('%%.%dg',nDigit+8);			% it will be of the type '%.Ng'
	set(handles.edit_x_min,'String',sprintf(frmt,head(1)))
	set(handles.edit_x_max,'String',sprintf(frmt,head(2)))
	set(handles.edit_y_min,'String',sprintf(frmt,head(3)))
	set(handles.edit_y_max,'String',sprintf(frmt,head(4)))
	handles.x_min = head(1);			handles.x_max = head(2);
	handles.y_min = head(3);			handles.y_max = head(4);
	handles.x_inc = head(8);			handles.y_inc = head(9);

	[m,n] = size(getappdata(handMir.figure1,'dem_z'));

	% Fill in the x,y_inc and nrow,ncol boxes
	nDigit = round( log10(abs(max(head(8:9)))) );		% Number of digits of the integer part
	frmt = sprintf('%%.%dg',nDigit+10);		% it will be of the type '%.Ng'
	set(handles.edit_Nrows,'String',m);		set(handles.edit_Ncols,'String',n)
	set(handles.edit_y_inc,'String',sprintf(frmt,head(9)))
	set(handles.edit_x_inc,'String',sprintf(frmt,head(8)))
	%-----------

	% If non-grid use image dims to estimate nRows, nCols
	if (m == 0)
		dim_funs('xInc', handles.edit_x_inc, handles)
		dim_funs('yInc', handles.edit_y_inc, handles)
		handles.nrows = round(str2double(get(handles.edit_Nrows,'String')));
		handles.ncols = round(str2double(get(handles.edit_Ncols,'String')));
		handles = guidata(hObject);		% It was changed inside dim_funs
	else
		handles.nrows = m;      handles.ncols = n;
	end

	%------------ Give a Pro look (3D) to the frame boxes  -------------------------------
	bgcolor = get(0,'DefaultUicontrolBackgroundColor');
	framecolor = max(min(0.65*bgcolor,[1 1 1]),[0 0 0]);
	h_f = findobj(hObject,'Style','Frame');
	for i=1:length(h_f)
		frame_size = get(h_f(i),'Position');
		f_bgc = get(h_f(i),'BackgroundColor');
		usr_d = get(h_f(i),'UserData');
		if abs(f_bgc(1)-bgcolor(1)) > 0.01           % When the frame's background color is not the default's
			frame3D(hObject,frame_size,framecolor,f_bgc,usr_d)
		else
			frame3D(hObject,frame_size,framecolor,'',usr_d)
			delete(h_f(i))
		end
	end
	% Recopy the text fields on top of previously created frames (uistack is to slow)
	h_t = [handles.txtFGeom handles.txtDGeom handles.txtGGeom handles.txtIGPos];
	for i=1:length(h_t)
		usr_d = get(h_t(i),'UserData');
		t_size = get(h_t(i),'Position');   t_str = get(h_t(i),'String');    fw = get(h_t(i),'FontWeight');
		bgc = get (h_t(i),'BackgroundColor');   fgc = get (h_t(i),'ForegroundColor');
		t_just = get(h_t(i),'HorizontalAlignment');     t_tag = get (h_t(i),'Tag');
		uicontrol('Parent',hObject, 'Style','text', 'Position',t_size,'String',t_str,'Tag',t_tag,...
			'BackgroundColor',bgc,'ForegroundColor',fgc,'FontWeight',fw,...
			'UserData',usr_d,'HorizontalAlignment',t_just);
	end
	delete(h_t)
	%------------- END Pro look (3D) -------------------------------------------------------

	guidata(hObject, handles);
	set(hObject,'Visible','on');
	if (nargout),   varargout{1} = hObject;     end

% ------------------------------------------------------------------------------------
function handles = edit_FaultWidth_Callback(hObject, eventdata, handles)
	% Actualize the "FaultWidth" field. EVENTDATA may not be empty
	if (nargout)
		xx = eventdata;
    else
		xx = str2double(get(hObject,'String'));
	end
	[fault,seg] = getFaultSeg(handles);
	if (xx < 0)         % If user tried to give a negative width
		xx = -xx;
		set(hObject,'String',num2str(xx))
	end
	dip = str2double(get(handles.edit_FaultDip,'String'));
	top_d = str2double(get(handles.edit_FaultTopDepth,'String'));
	depth = top_d + xx * cos((90-dip)*pi/180);
	set(handles.edit_FaultDepth,'String',num2str(depth));
	handles.FaultWidth{fault}(seg) = xx;
	handles.FaultDepth{fault}(seg) = depth;

	% Update the patch that represents the surface projection of the fault plane
	xx = [handles.fault_x{fault}(seg); handles.fault_x{fault}(seg+1)];
	yy = [handles.fault_y{fault}(seg); handles.fault_y{fault}(seg+1)];

	D2R = pi / 180;
	off = handles.FaultWidth{fault}(seg) * cos(handles.FaultDip{fault}(seg)*D2R);
	strk = handles.FaultStrike{fault}(seg);

	if (handles.geog)
		rng = off / 6371 / D2R;
		[lat1,lon1] = circ_geo(yy(1),xx(1),rng,strk+90,1);
		[lat2,lon2] = circ_geo(yy(2),xx(2),rng,strk+90,1);
	else
		if (handles.is_meters), off = off * 1e3;    end
		lon1 = xx(1) + off * cos(strk*D2R);     lon2 = xx(2) + off * cos(strk*D2R);
		lat1 = yy(1) - off * sin(strk*D2R);     lat2 = yy(2) - off * sin(strk*D2R);
	end
	x = [xx(1) xx(2) lon2 lon1 xx(1)];    y = [yy(1) yy(2) lat2 lat1 yy(1)];
	hp = getappdata(handles.h_fault(fault),'PatchHand');
	try     set(hp(seg),'XData',x,'YData',y,'FaceColor',[.8 .8 .8],'EdgeColor','k','LineWidth',1);  end
	
	z = -[top_d top_d depth depth top_d];
	if ( diff(handles.head(5:6)) > 10 ),	z = z * 1000;		end		% Assume grid's depth is in meters
	z = z + handles.head(5);				% CRUDE. It should be mean depth along the fault's length
	set(hp, 'UserData', z)					% So that we can Flederize it in 3D 

	handles = compMag(handles, fault);      % Compute and update Fault's Mw magnitude
	guidata(handles.figure1, handles);

% ------------------------------------------------------------------------------------
function edit_FaultStrike_Callback(hObject, eventdata, handles)
	% Cannot be changed

% ------------------------------------------------------------------------------------
function edit_FaultDip_Callback(hObject, eventdata, handles)
	% Actualize the "FaultDip" field
	xx = str2double(get(hObject,'String'));
	[fault,seg] = getFaultSeg(handles);
	top_d = str2double(get(handles.edit_FaultTopDepth,'String'));
	W = str2double(get(handles.edit_FaultWidth,'String'));
	depth = top_d + W * cos((90-xx)*pi/180);
	set(handles.edit_FaultDepth,'String',num2str(depth));
	handles.FaultDip{fault}(seg) = xx;
	handles.FaultDepth{fault}(seg) = depth;

	% Update the patch that represents the surface projection of the fault plane
	xx = [handles.fault_x{fault}(seg); handles.fault_x{fault}(seg+1)];
	yy = [handles.fault_y{fault}(seg); handles.fault_y{fault}(seg+1)];

	D2R = pi / 180;
	off = handles.FaultWidth{fault}(seg) * cos(handles.FaultDip{fault}(seg)*D2R);
	strk = handles.FaultStrike{fault}(seg);

	if (handles.geog)
		rng = off / 6371 / D2R;
		[lat1,lon1] = circ_geo(yy(1),xx(1),rng,strk+90,1);
		[lat2,lon2] = circ_geo(yy(2),xx(2),rng,strk+90,1);
	else
		if (handles.is_meters), off = off * 1e3;    end
		lon1 = xx(1) + off * cos(strk*D2R);     lon2 = xx(2) + off * cos(strk*D2R);
		lat1 = yy(1) - off * sin(strk*D2R);     lat2 = yy(2) - off * sin(strk*D2R);
	end
	x = [xx(1) xx(2) lon2 lon1 xx(1)];    y = [yy(1) yy(2) lat2 lat1 yy(1)];
	hp = getappdata(handles.h_fault(fault),'PatchHand');
	try     set(hp(seg),'XData',x,'YData',y,'FaceColor',[.8 .8 .8],'EdgeColor','k','LineWidth',1);  end

	z = -[top_d top_d depth depth top_d];
	if ( diff(handles.head(5:6)) > 10 ),	z = z * 1000;		end		% Assume grid's depth is in meters
	z = z + handles.head(5);				% CRUDE. It should be mean depth along the fault's length
	set(hp, 'UserData', z)					% So that we can Flederize it in 3D 

	guidata(hObject, handles);

% ------------------------------------------------------------------------------------
function edit_FaultDepth_Callback(hObject, eventdata, handles)
	% Actualize the "FaultTopDepth" field
	xx = str2double(get(hObject,'String'));
	if (xx < 0)         % If user tried to give a negative depth
		xx = -xx;
		set(hObject,'String',num2str(xx))
	end
	W = str2double(get(handles.edit_FaultWidth,'String'));
	dip = str2double(get(handles.edit_FaultDip,'String'));
	top_d = xx - W * cos((90-dip)*pi/180);
	set(handles.edit_FaultTopDepth,'String',num2str(top_d));
	[fault,seg] = getFaultSeg(handles);
	handles.FaultDepth{fault}(seg) = xx;
	handles.FaultTopDepth{fault}(seg) = top_d;
	guidata(hObject, handles);

% ------------------------------------------------------------------------------------
function edit_FaultTopDepth_Callback(hObject, eventdata, handles)
	% Actualize the "FaultDepth" field
	xx = str2double(get(hObject,'String'));
	if (xx < 0)         % If user tried to give a negative depth
		xx = -xx;
		set(hObject,'String',num2str(xx))
	end
	W = str2double(get(handles.edit_FaultWidth,'String'));
	dip = str2double(get(handles.edit_FaultDip,'String'));
	depth = xx + W * cos((90-dip)*pi/180);
	set(handles.edit_FaultDepth,'String',num2str(depth));
	[fault,seg] = getFaultSeg(handles);
	handles.FaultTopDepth{fault}(seg) = xx;
	handles.FaultDepth{fault}(seg) = depth;
	guidata(hObject, handles);

% ------------------------------------------------------------------------------------
function popup_segment_Callback(hObject, eventdata, handles)
	seg = get(hObject,'Value');
	if (handles.n_faults > 1),  fault = get(handles.popup_fault,'Value');
	else                        fault = 1;
	end

	% Fault parameters
	set(handles.edit_FaultLength,'String',num2str(handles.FaultLength{fault}(seg)))
	set(handles.edit_FaultStrike,'String',sprintf('%.1f',handles.FaultStrike{fault}(seg)))
	
	if (isnan(handles.FaultWidth{fault}(seg))),    str = '';
	else    str = num2str(handles.FaultWidth{fault}(seg));
	end
	set(handles.edit_FaultWidth,'String',str)
	
	set(handles.edit_FaultDip,'String',num2str(handles.FaultDip{fault}(seg),'%.1f'))
	set(handles.edit_FaultTopDepth,'String',num2str(handles.FaultTopDepth{fault}(seg)))
	
	if (isnan(handles.FaultDepth{fault}(seg))),    str = '';
	else    str = num2str(handles.FaultDepth{fault}(seg));
	end
	set(handles.edit_FaultDepth,'String',str)
	
	% Dislocation parameters
	set(handles.edit_DislocStrike,'String',sprintf('%.1f',handles.DislocStrike{fault}(seg)))
	if (isnan(handles.DislocSlip{fault}(seg))),    str = '';
	else    str = num2str(handles.DislocSlip{fault}(seg));
	end
	set(handles.edit_DislocSlip,'String',str)
	if (isnan(handles.DislocRake{fault}(seg))),    str = '';
	else    str = sprintf('%.1f',handles.DislocRake{fault}(seg));
	end
	set(handles.edit_DislocRake,'String',str)
	
	set(handles.edit_ux,'String',num2str(handles.ux{fault}(seg)))
	set(handles.edit_uy,'String',num2str(handles.uy{fault}(seg)))
	set(handles.edit_uz,'String',num2str(handles.uz{fault}(seg)))

% -----------------------------------------------------------------------------------------
function popup_fault_Callback(hObject, eventdata, handles)
	fault = get(hObject,'Value');
	S = cell(handles.nvert(fault),1);
	s_format = ['%.' num2str(fix(log10(handles.nvert(fault)))+1) 'd'];
	for (i=1:handles.nvert(fault)),     S{i} = ['Segment ' sprintf(s_format,i)];   end
	set(handles.popup_segment,'String',S,'Value',1)    
	seg = 1;    % Make current the first segment

	% Identify the currently active fault by setting its linestyle to dash
	set(handles.h_fault,'LineStyle','-')
	set(handles.h_fault(fault),'LineStyle','--')

	% Set the hide planes checkbox with the correct value for this fault
	if (handles.hide_planes(fault))
		set(handles.checkbox_hideFaultPlanes,'Value',1)
	else
		set(handles.checkbox_hideFaultPlanes,'Value',0)
	end

	% Fault parameters
	set(handles.edit_FaultLength,'String',num2str(handles.FaultLength{fault}(seg)))
	set(handles.edit_FaultStrike,'String',sprintf('%.1f',handles.FaultStrike{fault}(seg)))

	if (isnan(handles.FaultWidth{fault}(seg))),    str = '';
	else    str = num2str(handles.FaultWidth{fault}(seg));
	end
	set(handles.edit_FaultWidth,'String',str)

	set(handles.edit_FaultDip,'String',sprintf('%.1f',handles.FaultDip{fault}(seg)))
	set(handles.edit_FaultTopDepth,'String',num2str(handles.FaultTopDepth{fault}(seg)))

	if (isnan(handles.FaultDepth{fault}(seg))),    str = '';
	else    str = num2str(handles.FaultDepth{fault}(seg));
	end
	set(handles.edit_FaultDepth,'String',str)

	% Dislocation parameters
	set(handles.edit_DislocStrike,'String',sprintf('%.1f',handles.DislocStrike{fault}(seg)))
	if (isnan(handles.DislocSlip{fault}(seg))),    str = '';
	else    str = num2str(handles.DislocSlip{fault}(seg));
	end
	set(handles.edit_DislocSlip,'String',str)
	if (isnan(handles.DislocRake{fault}(seg))),    str = '';
	else    str = sprintf('%.1f',handles.DislocRake{fault}(seg));
	end
	set(handles.edit_DislocRake,'String',str)
	if (isnan(handles.ux{fault}(seg))),    str = '';
	else    str = num2str(handles.ux{fault}(seg));
	end
	set(handles.edit_ux,'String',str)
	if (isnan(handles.uy{fault}(seg))),    str = '';
	else    str = num2str(handles.uy{fault}(seg));
	end
	set(handles.edit_uy,'String',str)
	if (isnan(handles.uz{fault}(seg))),    str = '0';
	else    str = num2str(handles.uz{fault}(seg));
	end
	set(handles.edit_uz,'String',str)

	if (handles.Mw(fault) > 0)
		txt = ['Mw Magnitude = ' sprintf('%.1f',handles.Mw(fault))];
		set(handles.h_txt_Mw,'String',txt,'Position',handles.txt_Mw_pos + [0 0 30 0])
	else
		set(handles.h_txt_Mw,'String','Mw Magnitude = ','Position',handles.txt_Mw_pos)
	end
	refresh(handles.h_calling_fig);         % otherwise, ML BUG

% ---------------------------------------------------------------
function popup_GridCoords_Callback(hObject, eventdata, handles)
	xx = get(hObject,'Value');
	if (xx == 1),       handles.geog = 1;       handles.is_meters = 0;  handles.is_km = 0;
	elseif (xx == 2),   handles.is_meters = 1;  handles.is_geog = 0;    handles.is_km = 0;
	elseif (xx == 3),   handles.is_km = 1;      handles.is_geog = 0;    handles.is_meters = 0;
	end
	guidata(hObject,handles)

% ------------------------------------------------------------------------------------
function edit_DislocStrike_Callback(hObject, eventdata, handles)
	xx = str2double(get(hObject,'String'));
	[fault,seg] = getFaultSeg(handles);
	if (isnan(xx)),     set(hObject,'String',handles.DislocStrike{fault}(seg));   return,	end
	handles = convGeometry(handles, fault, seg, 'Aki');
	handles.DislocStrike{fault}(seg) = xx;
	guidata(hObject, handles);

% ------------------------------------------------------------------------------------
function edit_DislocRake_Callback(hObject, eventdata, handles)
	xx = str2double(get(hObject,'String'));
	[fault,seg] = getFaultSeg(handles);
	if (isnan(xx)),     set(hObject,'String',handles.DislocRake{fault}(seg));   return;     end
	handles = convGeometry(handles, fault, seg, 'Aki');
	handles.DislocRake{fault}(seg) = xx;
	guidata(hObject, handles);

% ------------------------------------------------------------------------------------
function edit_DislocSlip_Callback(hObject, eventdata, handles)
	xx = str2double(get(hObject,'String'));
	[fault,seg] = getFaultSeg(handles);
	if (isnan(xx)),     set(hObject,'String',handles.DislocSlip{fault}(seg));   return;     end
	handles = convGeometry(handles, fault, seg, 'Aki');
	handles.DislocSlip{fault}(seg) = xx;
	handles = compMag(handles, fault);      % Compute and update Fault's Mw magnitude
	guidata(hObject, handles);

% ------------------------------------------------------------------------------------
function edit_ux_Callback(hObject, eventdata, handles)
	xx = str2double(get(hObject,'String'));
	[fault,seg] = getFaultSeg(handles);
	if (isnan(xx)),     set(hObject,'String',handles.ux{fault}(seg));   return;     end
	handles = convGeometry(handles, fault, seg, 'Us');
	handles.ux{fault}(seg) = xx;
	guidata(hObject, handles);

% ------------------------------------------------------------------------------------
function edit_uy_Callback(hObject, eventdata, handles)
	xx = str2double(get(hObject,'String'));
	[fault,seg] = getFaultSeg(handles);
	if (isnan(xx)),     set(hObject,'String',handles.uy{fault}(seg));   return;     end
	handles = convGeometry(handles, fault, seg, 'Us');
	handles.uy{fault}(seg) = xx;
	guidata(hObject, handles);

% ------------------------------------------------------------------------------------
function edit_uz_Callback(hObject, eventdata, handles)
	xx = str2double(get(hObject,'String'));
	[fault,seg] = getFaultSeg(handles);
	if (isnan(xx)),     set(hObject,'String',handles.uz{fault}(seg));   return;     end
	handles = convGeometry(handles, fault, seg, 'Us');
	handles.uz{fault}(seg) = xx;
	guidata(hObject, handles);

% ------------------------------------------------------------------------------------
function [fault,seg] = getFaultSeg(handles)
	fault = 1;
	if (handles.n_faults > 1),	fault = get(handles.popup_fault,'Value');	end
	seg = get(handles.popup_segment,'Value');

% ------------------------------------------------------------------------------------
function handles = convGeometry(handles, fault, seg, opt)
	% Convert between the Aki & Richards and Ux, Uy, Uz dislocation conventions
    
	D2R = pi / 180;
	if (strcmp(opt, 'Aki'))
		f_strike = str2double(get(handles.edit_FaultStrike,'String'));      % Fault strike
		d_strike = str2double(get(handles.edit_DislocStrike,'String'));     % Dislocation strike
		slip = str2double(get(handles.edit_DislocSlip,'String'));           % Dislocation slip
		rake = str2double(get(handles.edit_DislocRake,'String'));           % Dislocation rake
		ux = slip * cos(rake*D2R) * cos((f_strike - d_strike)*D2R);
		uy = slip * sin(rake*D2R) * cos((f_strike - d_strike)*D2R);
		uz = slip * sin((f_strike - d_strike)*D2R);
		if (abs(ux) < 1e-6),    ux = 0;     end
		if (abs(uy) < 1e-6),    uy = 0;     end
		if (abs(uz) < 1e-6),    uz = 0;     end
		set(handles.edit_ux,'String',num2str(ux))
		set(handles.edit_uy,'String',num2str(uy))
		set(handles.edit_uz,'String',num2str(uz))
		handles.ux{fault}(seg) = ux;	handles.uy{fault}(seg) = uy;	handles.uz{fault}(seg) = uz;
	else
		ux = str2double(get(handles.edit_ux,'String'));
		uy = str2double(get(handles.edit_uy,'String'));
		uz = str2double(get(handles.edit_uz,'String'));
		slip = sqrt(ux^2 + uy^2 + uz^2);
		rake = atan2(uy, ux) / D2R;
		set(handles.edit_DislocSlip,'String',num2str(slip))
		set(handles.edit_DislocRake,'String',sprintf('%.1f',rake))
		handles.DislocSlip{fault}(seg) = slip;
		handles.DislocRake{fault}(seg) = rake;
	end

% -------------------------------------------------------------------------------------
function edit_x_min_Callback(hObject, eventdata, handles)
	dim_funs('xMin', hObject, handles)

% -------------------------------------------------------------------------------------
function edit_x_max_Callback(hObject, eventdata, handles)
	dim_funs('xMax', hObject, handles)

% --------------------------------------------------------------------
function edit_y_min_Callback(hObject, eventdata, handles)
	dim_funs('yMin', hObject, handles)

% --------------------------------------------------------------------
function edit_y_max_Callback(hObject, eventdata, handles)
	dim_funs('yMax', hObject, handles)

% --------------------------------------------------------------------
function edit_x_inc_Callback(hObject, eventdata, handles)
	dim_funs('xInc', hObject, handles)

% --------------------------------------------------------------------
function edit_Ncols_Callback(hObject, eventdata, handles)
	dim_funs('nCols', hObject, handles)

% --------------------------------------------------------------------
function edit_y_inc_Callback(hObject, eventdata, handles)
	dim_funs('yInc', hObject, handles)

% --------------------------------------------------------------------
function edit_Nrows_Callback(hObject, eventdata, handles)
	dim_funs('nRows', hObject, handles)

% ------------------------------------------------------------------------------------
function pushbutton_Help_R_Callback(hObject, eventdata, handles)
message = {'That''s prety obvious to guess what this option does. You select an area,'
    'the grid spacing or the number of rows/columns and the deformation will'
    'be computed at all nodes of that grid.'};
helpdlg(message,'Help on deformation grid');

% ------------------------------------------------------------------------------------
function pushbutton_Help_H_Callback(hObject, eventdata, handles)
	message = {'If you have a file with x,y positions, then the deformation will be computed at those postions'};
	helpdlg(message,'Little Help');

% ------------------------------------------------------------------------------------
function edit_InputFile_Callback(hObject, eventdata, handles)
fname = get(hObject,'String');
if isempty(fname),  handles.input_locations = [];   return;   end
hFig = gcf;
[bin,n_column,multi_seg,n_headers] = guess_file(fname);
% If msgbox exist we have to move it from behind the main window. So get it's handle
hMsgFig = gcf;
if (hFig ~= hMsgFig),       figure(hMsgFig);   end   % If msgbox exists, bring it forward
% If error in reading file
if isempty(bin) && isempty(n_column) && isempty(multi_seg) && isempty(n_headers)
    errordlg(['Error reading file ' fname],'Error');    return
end
if multi_seg ~= 0   % multisegments are not spported
    errordlg('Multisegment files are yet not supported.','Error');   return
end
if (bin == 0)   % ASCII
    if n_column < 2
        errordlg('File error. Your file doesn''t have at least 2 columns','Error'); return
    end
    handles.input_locations = read_xy(fname,n_column,n_headers);
    if (hFig ~= hMsgFig);       figure(hFig);   end     % gain access to the drawing figure
    nr = size(handles.input_locations,1);
    if (nr == 0)
        errordlg('Your file is empty.','Chico Clever');   return
    end
    if (n_headers > 0)      % We have headers in file (ai!, ai!)
        set(handles.checkbox_Option_H,'Value',1)
        set(handles.edit_nHeaders,'String',num2str(n_headers))
    end
else        % BINARY
    errordlg('Sorry, reading binary files is not yet programed','Error');   return
end
guidata(hObject,handles)

% ------------------------------------------------------------------------------------
function pushbutton_InputFile_Callback(hObject, eventdata, handles)
if (~isempty(handles.h_calling_fig))                    % If we know the handle to the calling fig
    cfig_handles = guidata(handles.h_calling_fig);      % get handles of the calling fig
    last_dir = cfig_handles.last_dir;
    home = cfig_handles.home_dir;
else
    last_dir = [];
end

if (~isempty(last_dir)),    cd(last_dir);   end
[FileName,PathName] = uigetfile({'*.dat;*.DAT;*.xy', 'Maregraph location (*.dat,*.DAT,*.xy)';'*.*', 'All Files (*.*)'},'Select Maregraphs position');
if (~isempty(last_dir)),    cd(home);   end
if isequal(FileName,0);     return;     end
fname = [PathName FileName];
hFig = gcf;
[bin,n_column,multi_seg,n_headers] = guess_file(fname);
% If msgbox exist we have to move it from behind the main window. So get it's handle
hMsgFig = gcf;
if (hFig ~= hMsgFig),       figure(hMsgFig);   end   % If msgbox exists, bring it forward
% If error in reading file
if isempty(bin) && isempty(n_column) && isempty(multi_seg) && isempty(n_headers)
    errordlg(['Error reading file ' fname],'Error');    return
end
if multi_seg ~= 0   % multisegments are not spported
    errordlg('Multisegment files are yet not supported.','Error');   return
end
if (bin == 0)   % ASCII
    if (n_column < 2)
        errordlg('File error. Your file doesn''t have at least 2 columns','Error'); return
    end
    handles.input_locations = read_xy(fname,n_column,n_headers);
    if (hFig ~= hMsgFig);       figure(hFig);   end     % gain access to the drawing figure
    nr = size(handles.input_locations,1);
    if (nr == 0)
        errordlg('Your file is empty.','Chico Clever');   return
    end
    if (n_headers > 0)      % We have headers in file (ai!, ai!)
        set(handles.checkbox_Option_H,'Value',1)
        set(handles.edit_nHeaders,'String',num2str(n_headers))
    end
else        % BINARY
    errordlg('Sorry, reading binary files is not yet programed','Error');   return
end
set(handles.edit_InputFile,'String',fname)
guidata(hObject,handles)

% ------------------------------------------------------------------------------------
function edit_sx_Callback(hObject, eventdata, handles)
	if ( isempty(get(hObject,'String')) ),  set(hObject,'String','0');   return;     end

% ------------------------------------------------------------------------------------
function edit_sy_Callback(hObject, eventdata, handles)
	if ( isempty(get(hObject,'String')) ),  set(hObject,'String','0');   return;     end

% ------------------------------------------------------------------------------------
function edit_sz_Callback(hObject, eventdata, handles)
	if ( isempty(get(hObject,'String')) ),  set(hObject,'String','0');   return;     end

% ------------------------------------------------------------------------------------
function pushbutton_compute_Callback(hObject, eventdata, handles)
	% If cartesian coordinates, they must be in meters
	if (any(isnan(cat(1,handles.FaultWidth{:}))))
		errordlg('One or more segments where not set with the fault''s Width','Error');    return
	end
	if (any(isnan(cat(1,handles.FaultDepth{:}))))
		errordlg('One or more segments where not set with the fault''s Depth','Error');    return
	end
	if (any(isnan(cat(1,handles.DislocSlip{:}))))
		errordlg('One or more segments where not set with the movement''s slip','Error');    return
	end
	
	if( all(cat(1,handles.ux{:}) == 0) && all(cat(1,handles.uy{:}) == 0) && all(cat(1,handles.uz{:}) == 0) )
		errordlg('No movement along faults(s), nothing to compute there.','Error')
		return
	end

	sx = str2double(get(handles.edit_sx,'String'));
	sy = str2double(get(handles.edit_sy,'String'));
	sz = str2double(get(handles.edit_sz,'String'));
	s = [sx sy sz];

	if( sx == 0 && sy == 0 && sz == 0 )
		errordlg('Looking vector is looking nowhere.','Error')
		return
	elseif (sqrt(sx^2 + sy^2 + sz^2) < 0.99)
		warndlg('Norm of looking vector is less than 1. It means you are loosing deformation.','Warning')
	end

	% Get grid params
	xmin = str2double(get(handles.edit_x_min,'String'));     xmax = str2double(get(handles.edit_x_max,'String'));
	ymin = str2double(get(handles.edit_y_min,'String'));     ymax = str2double(get(handles.edit_y_max,'String'));
	nrow = str2double(get(handles.edit_Nrows,'String'));    ncol = str2double(get(handles.edit_Ncols,'String'));

	x = handles.fault_x;    y = handles.fault_y;
	if (~iscell(x)),        x = {x};    y = {y};    end
	fig_xlim = [xmin xmax];   fig_ylim = [ymin ymax];
	to_km = 1;      % The conversion from m->km will be done inside range_change

if (handles.geog)    
    opt_M = '-M';
else
    if (handles.is_meters),     to_km = 1000;   end
    for (i=1:handles.n_faults)
        x{i} = x{i} / to_km;   y{i} = y{i} / to_km;
        handles.FaultLength{i} = handles.FaultLength{i} / to_km;
    end
    opt_M = '';
end

for (i=1:handles.n_faults)
    % I have to do fish the patch coords because range_change does not seams to
    % use the fault trace coords but the coordinates of the fault at its depth 
    hp = getappdata(handles.h_fault(i),'PatchHand');
    xp = get(hp,'XData');    yp = get(hp,'YData');
    if (iscell(xp))
        x{i} = [];    y{i} = [];
        for (k=1:length(xp))
            x{i} = [x{i}; xp{k}(4)/to_km];
            y{i} = [y{i}; yp{k}(4)/to_km];
        end
    else
        x{i} = xp(4)/to_km;   y{i} = yp(4)/to_km;
    end
end

if (isempty(handles.input_locations))   % If ground positions were not given, compute a grid
    E = linspace(fig_xlim(1),fig_xlim(2),ncol)/to_km;
    N = linspace(fig_ylim(1),fig_ylim(2),nrow)/to_km;
    N = N(:);               % From the rngchn example, y coords are in a column vector

    % Compute deformation
    %U = range_change(x,y,strike,depth,dip,ux,uy,uz,L,W,E,N,s);
    if (handles.n_faults > 1)           % We have multiple faults
        U = zeros(nrow,ncol);
        h = waitbar(0,'Computing deformation');
        for (i=1:handles.n_faults)
            waitbar(i/handles.n_faults)
            U0 = range_change(x{i}(:),y{i}(:),handles.FaultStrike{i}(:),handles.FaultDepth{i}(:),handles.FaultDip{i}(:),...
                handles.ux{i}(:),handles.uy{i}(:),handles.uz{i}(:),handles.FaultLength{i}(:),handles.FaultWidth{i}(:),...
                E,N,s,opt_M);
            U = U0 + U;
        end
        close(h);    clear U0;
    else                                % We have only one fault
        U = range_change(x{1}(:),y{1}(:),handles.FaultStrike{1}(:),handles.FaultDepth{1}(:),handles.FaultDip{1}(:),...
            handles.ux{1}(:),handles.uy{1}(:),handles.uz{1}(:),handles.FaultLength{1}(:),handles.FaultWidth{1}(:),...
            E,N,s,opt_M);
    end

    z_max = max(U(:));     z_min = min(U(:));
    dx = str2double(get(handles.edit_x_inc,'String'));
    dy = str2double(get(handles.edit_y_inc,'String'));

    if (handles.geog)
        head.head = [xmin xmax ymin ymax z_min z_max 0 dx dy];
        head.X = linspace(xmin,xmax,ncol);
        head.Y = linspace(ymin,ymax,nrow);
    else
        E = E * to_km;   N = N * to_km;     % Convert to grid coords
        head.head = [E(1) E(end) N(1) N(end) z_min z_max 0 E(2)-E(1) N(2)-N(1)];
        head.X = E;     head.Y = N;
    end
    
    % Test if the result seams correct in terms of size
    [m,n] = size(U);        hWarn = [];
    if ( (m ~= nrow) || (n ~= ncol) )
        msg{1} = 'Someting went wrong. Output file has not the required size. Maybe a meters<->kilometers bad guess?';
        if (abs(dx - dy) > 1e-5)
            msg{2} = ' ';
            msg{3} = 'No. Almost likely this was due to the fact that or X and Y spacings are diferent.';
        end
        hWarn = warndlg(msg,'Warning');
    end

	% SHOW WHAT WE HAVE GOT
    if get(handles.radiobutton_deformation,'Value')		% Show deformation 
        mirone(U,head,'Deformation',handles.h_calling_fig);
    else												% Show Interferogram
		val = get(handles.list_Interfero,'Val');
		switch val
			case 1,		cdo = 28.4;				% mm
			case 2,		cdo = 28.4 / 10;		% cm
			case 3,		cdo = 28.4 / 1000;		% meters
		end
        mirone(U,head,'Interfero',cdo);
    end
    if (~isempty(hWarn)),   figure(hWarn);      end
    
else        % Ground positions were given
    E = handles.input_locations(:,1)/to_km;
    N = handles.input_locations(:,2)/to_km;
    if (handles.n_faults > 1)           % We have multiple faults
        U = 0;
        for (i=1:handles.n_faults)
            if (get(handles.checkbox_ToggleXY,'Value'))
                U0 = range_change(x{i}(:),y(:),handles.FaultStrike(:),handles.FaultDepth(:),handles.FaultDip(:),...
                    handles.ux(:),handles.uy(:),handles.uz(:),handles.FaultLength(:),handles.FaultWidth(:),...
                    N,E,s,opt_M);
            else
                U0 = range_change(x(:),y(:),handles.FaultStrike(:),handles.FaultDepth(:),handles.FaultDip(:),...
                    handles.ux(:),handles.uy(:),handles.uz(:),handles.FaultLength(:),handles.FaultWidth(:),...
                    E,N,s,opt_M);
            end
            U = U0 + U;
        end
        clear U0;
    else                                % We have only one fault
        if (get(handles.checkbox_ToggleXY,'Value'))
            U = range_change(x{1}(:),y{1}(:),handles.FaultStrike{1}(:),handles.FaultDepth{1}(:),handles.FaultDip{1}(:),...
                handles.ux(:),handles.uy{1}(:),handles.uz{1}(:),handles.FaultLength{1}(:),handles.FaultWidth{1}(:),...
                N,E,s,opt_M);
        else
            U = range_change(x{1}(:),y{1}(:),handles.FaultStrike{1}(:),handles.FaultDepth{1}(:),handles.FaultDip{1}(:),...
                handles.ux(:),handles.uy{1}(:),handles.uz{1}(:),handles.FaultLength{1}(:),handles.FaultWidth{1}(:),...
                E,N,s,opt_M);
        end
    end
    [FileName,PathName] = uiputfile({'*.dat;*.DAT', 'Data files (*.dat,*.DAT)';'*.*', 'All Files (*.*)'},'Select output deformation file');
    if isequal(FileName,0);     return;      end     % User gave up
    double2ascii([PathName FileName],[U(:,1) U(:,2) U(:,3)],'%f\t%f\t%f');     %NAO SEI SE E ASSIM. FALTA TESTAR
end

% ------------------------------------------------------------------------------------
function [lonlim,zone] = utmorigin(lon)
	% Returns the UTM longitude limits and the Zone. Note that there is no way of telling
	% if it is a North or South zone, but that should be easy by knowing the latitude.
	lons = (-180:6:180)';	
	ind = find(lons <= lon);  lonsidx = ind(max(ind));
		
	if (lonsidx < 1 || lonsidx > 61),    lonsidx = [];
	elseif (lonsidx == 61),             lonsidx = 60;    end

	zone = num2str(lonsidx);

	if (length(zone) == 1)
		lonsidx = str2double(zone);
	elseif (length(zone) == 2)
		num = str2double(zone);
		if isempty(num)
			lonsidx = str2double(zone(1));
		else
			lonsidx = num;
		end
	end
	lonlims = [(-180:6:174)' (-174:6:180)'];
	lonlim = lonlims(lonsidx,:);

% ------------------------------------------------------------------------------------
function radiobutton_deformation_Callback(hObject, eventdata, handles)
	if get(hObject,'Value')
		set(handles.radiobutton_interfero,'Value',0)
		set(handles.edit_sx,'String',num2str(handles.sx))
		set(handles.edit_sy,'String',num2str(handles.sy))
		set(handles.edit_sz,'String',num2str(handles.sz))
		set(handles.list_Interfero,'Vis','off')
	else
		set(hObject,'Value',1)
	end

% ------------------------------------------------------------------------------------
function radiobutton_interfero_Callback(hObject, eventdata, handles)
	% Use ERS looking vector
	if get(hObject,'Value')
		set(handles.radiobutton_deformation,'Value',0)
		set(handles.edit_sx,'String','0.333')
		set(handles.edit_sy,'String','-0.07')
		set(handles.edit_sz,'String','0.94')
		set(handles.list_Interfero,'Vis','on')
	else
		set(hObject,'Value',1)
	end

% ------------------------------------------------------------------------------------
function pushbutton_cancel_Callback(hObject, eventdata, handles)
    figure1_CloseRequestFcn([], [], handles)

% -----------------------------------------------------------------------------------------
function len = LineLength(h,geog)
x = get(h,'XData');     y = get(h,'YData');
D2R = pi/180;	earth_rad = 6371;
len = [];
if (~iscell(x))
	if (geog)
        x = x * D2R;	y = y * D2R;
        lat_i = y(1:length(y)-1);   lat_f = y(2:length(y));     clear y;
        lon_i = x(1:length(x)-1);   lon_f = x(2:length(x));     clear x;
        tmp = sin(lat_i).*sin(lat_f) + cos(lat_i).*cos(lat_f).*cos(lon_f-lon_i);
        clear lat_i lat_f lon_i lon_f;
        len = [len; acos(tmp) * earth_rad];         % Distance in km
	else
        dx = diff(x);   dy = diff(y);
        len = [len; sqrt(dx.*dx + dy.*dy)];         % Distance in user unites
	end
else
	if (geog)
        for (k=1:length(x))
            xx = x{k} * D2R;    yy = y{k} * D2R;
            lat_i = yy(1:length(yy)-1);   lat_f = yy(2:length(yy));
            lon_i = xx(1:length(xx)-1);   lon_f = xx(2:length(xx));
            tmp = sin(lat_i).*sin(lat_f) + cos(lat_i).*cos(lat_f).*cos(lon_f-lon_i);
            len{k} = acos(tmp) * earth_rad;         % Distance in km
        end
	else
        for (k=1:length(x))
            xx = x{k};      yy = y{k};
            dx = diff(xx);  dy = diff(yy);
            len{k} = sqrt(dx.*dx + dy.*dy);         % Distance in user unites
        end
	end
end

% -----------------------------------------------------------------------------------------
function xy = read_xy(file,n_col,n_head)
	% build the format string to read the data n_columns
	xy = [];    format = [];    fid = fopen(file,'r');
	for (i=1:n_col),    format = [format '%f '];    end
	% Jump header lines
	for (i = 1:n_head),    tline = fgetl(fid);  end

	todos = fread(fid,'*char');
	xy = sscanf(todos,format,[n_col inf])';    % After hours strugling agains this FILHO DA PUTA, I may have found
	fclose(fid);

% -----------------------------------------------------------------------------------------
function checkbox_hideFaultPlanes_Callback(hObject, eventdata, handles)
	fault = getFaultSeg(handles);
	hp = getappdata(handles.h_fault(fault),'PatchHand');
	if (get(hObject,'Value'))
		try     set(hp,'Visible','off');    end
		handles.hide_planes(fault) = 1;
	else
		try     set(hp,'Visible','on');     end
		handles.hide_planes(fault) = 0;
	end
	guidata(hObject,handles)

% -----------------------------------------------------------------------------------------
function handles = set_all_faults(handles,varargin)
	% varargin contains a set of parameters of a Slip model transmited by fault_models.m  
	handles.h_calling_fig = varargin{1}.figure1;
	handles.h_fault = varargin{2};
	handles.FaultTopDepth = varargin{3};
	handles.FaultWidth = varargin{4};
	handles.FaultStrike = varargin{5};
	handles.DislocStrike = varargin{5};
	handles.DislocSlip = varargin{6};
	handles.FaultDip  = varargin{7};
	handles.DislocRake = varargin{8};

	nSeg = size(handles.h_fault,1);
	if (nSeg == 1)
		nx = numel(handles.FaultStrike{1});		n_fault = numel(handles.h_fault);
		
		% Pre-allocations
		handles.ux = repmat({zeros(1, nx)}, n_fault, 1);
		handles.uy = repmat({zeros(1, nx)}, n_fault, 1);
		handles.uz = repmat({zeros(1, nx)}, n_fault, 1);
		handles.FaultDepth = repmat({zeros(1, nx)}, n_fault, 1);
		
	else				% Multi-seg Slip model file. Old logic obliges to have {1,nFaultTotal}
		% Remember that multi-segs had declaration in fault_models.m like:
		% hLine = zeros(nSeg,nz); and 	strike = cell(nSeg,nz);
		% where each cell element contains "nPatch" values
		nz = size(handles.h_fault,2);
		handles.ux = [];
		for (k = 1:nSeg)
			nx = numel(handles.FaultStrike{k});
			handles.ux = [handles.ux; repmat({zeros(1, nx)}, nz, 1)];
		end
		handles.uy = handles.ux;
		handles.uz = handles.ux;
		handles.FaultDepth = handles.ux;
		
		% Now reshape them for the old logic
		handles.h_fault = reshape((varargin{2})',1,[]);
		handles.FaultTopDepth = reshape((varargin{3})',1,[]);
		handles.FaultWidth = reshape((varargin{4})',1,[]);
		handles.FaultStrike = reshape((varargin{5})',1,[]);
		handles.DislocStrike = reshape((varargin{5})',1,[]);
		handles.DislocSlip = reshape((varargin{6})',1,[]);
		handles.FaultDip  = reshape((varargin{7})',1,[]);
		handles.DislocRake = reshape((varargin{8})',1,[]);
		n_fault = numel(handles.DislocRake);
	end
	
	for (k=1:n_fault)
		handles.FaultDepth{k} = handles.FaultTopDepth{k} + handles.FaultWidth{k} .* ...
			sin(handles.FaultDip{k} * pi/180);
		handles.ux{k} = handles.DislocSlip{k} .* cos(handles.DislocRake{k} * pi/180);
		handles.uy{k} = handles.DislocSlip{k} .* sin(handles.DislocRake{k} * pi/180);
	end

	set(handles.edit_FaultWidth,'String',handles.FaultWidth{1}(1));
	set(handles.edit_FaultStrike,'String',num2str(handles.FaultStrike{1}(1)));
	set(handles.edit_FaultDip,'String',num2str(handles.FaultDip{1}(1)));
	set(handles.edit_FaultTopDepth,'String',num2str(handles.FaultTopDepth{1}(1)));
	set(handles.edit_FaultDepth,'String',num2str(handles.FaultDepth{1}(1)));
	set(handles.edit_DislocSlip,'String',num2str(handles.DislocSlip{1}(1)));
	set(handles.edit_DislocStrike,'String',num2str(handles.FaultStrike{1}(1)));
	set(handles.edit_DislocRake,'String',num2str(handles.DislocRake{1}(1)));
	set(handles.edit_ux,'String',handles.ux{1}(1));
	set(handles.edit_uy,'String',handles.uy{1}(1));
	set(handles.edit_uz,'String','0');

% ------------------------------------------------------------------------------------
function [handles, mag, M0] = compMag(handles, fault)
	% Compute Moment magnitude
	M0 = 3.0e10 * handles.um_milhao * handles.DislocSlip{fault}(:) .* handles.FaultWidth{fault}(:) .* ...
		handles.FaultLength{fault}(:);
	if (length(M0) > 1),    M0 = sum(M0);   end
	mag = 2/3*(log10(M0) - 9.1);
	if (~isnan(mag))
		txt = ['Mw Magnitude = ' sprintf('%.1f',mag)];
		set(handles.h_txt_Mw,'String',txt,'Position',handles.txt_Mw_pos + [0 0 30 0])
		handles.Mw(fault) = mag;
	end

% ------------------------------------------------------------------------------------
function figure1_CloseRequestFcn(hObject, eventdata, handles)
	if (~handles.fault_in)
		for (i=1:numel(handles.h_fault))
			hp = getappdata(handles.h_fault(i),'PatchHand');
			try     set(hp,'XData', [], 'YData',[],'Visible','on');     end     % Had to use a try (f.. endless errors)
		end
    end
	delete(handles.figure1)

% ------------------------------------------------------------------------------------
% --- Creates and returns a handle to the GUI figure. 
function deform_okada_LayoutFcn(h1)

set(h1,'Units','pixels',...
'PaperUnits',get(0,'defaultfigurePaperUnits'),...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'CloseRequestFcn',{@deform_okada_uicallback,h1,'figure1_CloseRequestFcn'},...
'MenuBar','none',...
'Name','Okada deformation',...
'NumberTitle','off',...
'Position',[520 415 536 365],...
'Renderer',get(0,'defaultfigureRenderer'),...
'RendererMode','manual',...
'Resize','off',...
'Tag','figure1',...
'HandleVisibility','callback');

uicontrol('Parent',h1,'Position',[320 224 211 131],'Style','frame','Tag','frame4');
uicontrol('Parent',h1,'Position',[10 224 181 131],'Style','frame','Tag','frame3');
uicontrol('Parent',h1,'Enable','inactive','Position',[11 15 350 67],'Style','frame','Tag','frame2');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Position',[20 311 71 21],...
'Style','edit',...
'TooltipString','Fault length (km)',...
'Tag','edit_FaultLength');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@deform_okada_uicallback,h1,'edit_FaultWidth_Callback'},...
'Position',[110 311 71 21],...
'Style','edit',...
'TooltipString','Fault width (km)',...
'Tag','edit_FaultWidth');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@deform_okada_uicallback,h1,'edit_FaultStrike_Callback'},...
'Position',[20 271 71 21],...
'Style','edit',...
'TooltipString','Fault strike (degrees)',...
'Tag','edit_FaultStrike');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@deform_okada_uicallback,h1,'edit_FaultDip_Callback'},...
'Position',[110 271 71 21],...
'Style','edit',...
'TooltipString','Fault dip (degrees)',...
'Tag','edit_FaultDip');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@deform_okada_uicallback,h1,'edit_FaultDepth_Callback'},...
'Position',[20 232 71 21],...
'Style','edit',...
'TooltipString','Depth of the base of fault''s plane',...
'Tag','edit_FaultDepth');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@deform_okada_uicallback,h1,'edit_FaultTopDepth_Callback'},...
'Position',[110 231 71 21],...
'Style','edit',...
'TooltipString','Alternatively, give depth to the fault''s top ',...
'Tag','edit_FaultTopDepth');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@deform_okada_uicallback,h1,'edit_DislocStrike_Callback'},...
'Position',[330 311 51 21],...
'Style','edit',...
'Tag','edit_DislocStrike');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@deform_okada_uicallback,h1,'edit_DislocRake_Callback'},...
'Position',[400 311 51 21],...
'Style','edit',...
'TooltipString','Displacement angle clock-wise from horizontal',...
'Tag','edit_DislocRake');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@deform_okada_uicallback,h1,'edit_DislocSlip_Callback'},...
'Position',[470 311 51 21],...
'Style','edit',...
'TooltipString','Total displacement',...
'Tag','edit_DislocSlip');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@deform_okada_uicallback,h1,'edit_ux_Callback'},...
'Position',[330 271 51 21],...
'Style','edit',...
'TooltipString','Left-lateral displacement along the fault plane (along strike)',...
'Tag','edit_ux');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@deform_okada_uicallback,h1,'edit_uy_Callback'},...
'Position',[400 271 51 21],...
'Style','edit',...
'TooltipString','Displacement up-dip the fault plane (across strike)',...
'Tag','edit_uy');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@deform_okada_uicallback,h1,'edit_uz_Callback'},...
'Position',[470 271 51 21],...
'Style','edit',...
'TooltipString','fault tensile slip',...
'Tag','edit_uz');

uicontrol('Parent',h1,'Enable','inactive','Position',[10 109 350 93],...
'String',{''},'Style','frame','Tag','frame1');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@deform_okada_uicallback,h1,'edit_sx_Callback'},...
'Position',[330 231 51 21],...
'Style','edit',...
'TooltipString','Component of unit vector along North coords',...
'Tag','edit_sx');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@deform_okada_uicallback,h1,'edit_sy_Callback'},...
'Position',[400 231 51 21],...
'Style','edit',...
'TooltipString','Component of unit vector along East coords',...
'Tag','edit_sy');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@deform_okada_uicallback,h1,'edit_sz_Callback'},...
'Position',[470 231 51 21],...
'Style','edit',...
'TooltipString','Component of unit vector along Vertical coords',...
'Tag','edit_sz');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@deform_okada_uicallback,h1,'edit_x_min_Callback'},...
'Position',[76 162 71 21],...
'Style','edit',...
'TooltipString','X min value',...
'Tag','edit_x_min');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@deform_okada_uicallback,h1,'edit_x_max_Callback'},...
'Position',[152 162 71 21],...
'Style','edit',...
'TooltipString','X max value',...
'Tag','edit_x_max');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@deform_okada_uicallback,h1,'edit_x_inc_Callback'},...
'Position',[228 162 71 21],...
'Style','edit',...
'TooltipString','DX grid spacing',...
'Tag','edit_x_inc');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@deform_okada_uicallback,h1,'edit_Ncols_Callback'},...
'Position',[304 162 45 21],...
'Style','edit',...
'TooltipString','Number of columns in the grid',...
'Tag','edit_Ncols');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@deform_okada_uicallback,h1,'edit_y_min_Callback'},...
'Position',[76 136 71 21],...
'Style','edit',...
'TooltipString','Y min value',...
'Tag','edit_y_min');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@deform_okada_uicallback,h1,'edit_y_max_Callback'},...
'Position',[152 136 71 21],...
'Style','edit',...
'TooltipString','Y max value',...
'Tag','edit_y_max');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@deform_okada_uicallback,h1,'edit_y_inc_Callback'},...
'Position',[228 136 71 21],...
'Style','edit',...
'TooltipString','DY grid spacing',...
'Tag','edit_y_inc');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@deform_okada_uicallback,h1,'edit_Nrows_Callback'},...
'Position',[304 136 45 21],...
'Style','edit',...
'TooltipString','Number of columns in the grid',...
'Tag','edit_Nrows');

uicontrol('Parent',h1,...
'BackgroundColor',[0.831372559070587 0.815686285495758 0.7843137383461],...
'Callback',{@deform_okada_uicallback,h1,'pushbutton_Help_R_Callback'},...
'FontWeight','bold',...
'ForegroundColor',[0 0 1],...
'Position',[289 114 61 18],...
'String','?',...
'Tag','pushbutton_Help_R');

uicontrol('Parent',h1,...
'Position',[22 54 65 15],...
'String','Headers?',...
'Style','checkbox',...
'TooltipString','Are there any header lines in the input file?',...
'Tag','checkbox_Option_H');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'HorizontalAlignment','left',...
'Position',[171 50 31 20],...
'String','1',...
'Style','edit',...
'TooltipString','How many?',...
'Tag','edit_nHeaders');

uicontrol('Parent',h1,...
'Position',[221 53 75 19],...
'String','Toggle x,y',...
'Style','checkbox',...
'TooltipString','Toggle x and y columns',...
'Tag','checkbox_ToggleXY');

uicontrol('Parent',h1,'Enable','inactive','Position',[18 167 55 15],...
'String','X Direction','Style','text');

uicontrol('Parent',h1,'Enable','inactive','Position',[17 141 55 15],...
'String','Y Direction','Style','text');

uicontrol('Parent',h1,'Enable','inactive','Position',[169 184 41 13],...
'String','Max','Style','text');

uicontrol('Parent',h1,'Enable','inactive','Position',[91 185 41 13],...
'String','Min','Style','text');

uicontrol('Parent',h1,'Enable','inactive','Position',[246 185 41 13],...
'String','Spacing','Style','text');

uicontrol('Parent',h1,'Enable','inactive','Position',[302 185 51 13],...
'String','# of lines','Style','text');

uicontrol('Parent',h1,'Enable','inactive','Position',[30 195 121 15],...
'String','Griding Line Geometry','Style','text','Tag','txtGGeom');

uicontrol('Parent',h1,...
'BackgroundColor',[0.831372559070587 0.815686285495758 0.7843137383461],...
'Callback',{@deform_okada_uicallback,h1,'pushbutton_Help_H_Callback'},...
'FontWeight','bold',...
'ForegroundColor',[0 0 1],...
'Position',[331 50 22 22],...
'String','?',...
'Tag','pushbutton_Help_H');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@deform_okada_uicallback,h1,'edit_InputFile_Callback'},...
'HorizontalAlignment','left',...
'Position',[21 22 310 22],...
'Style','edit',...
'TooltipString','File name with x, y positions where to compute deformation',...
'Tag','edit_InputFile');

uicontrol('Parent',h1,...
'Callback',{@deform_okada_uicallback,h1,'pushbutton_InputFile_Callback'},...
'Position',[331 21 23 23],...
'Tag','pushbutton_InputFile');

uicontrol('Parent',h1,'Enable','inactive','Position',[32 74 145 15],...
'String','Input Ground Positions File','Style','text','Tag','txtIGPos');

uicontrol('Parent',h1,'HorizontalAlignment','left','Position',[102 53 67 15],...
'String','Nº of headers','Style','text','TooltipString','How many?');

uicontrol('Parent',h1,'Enable','inactive','Position',[36 333 41 13],...
'String','Length','Style','text');

uicontrol('Parent',h1,'Enable','inactive','Position',[125 334 41 13],...
'String','Width','Style','text');

uicontrol('Parent',h1,'Enable','inactive','Position',[34 293 41 13],...
'String','Strike','Style','text');

uicontrol('Parent',h1,'Enable','inactive','Position',[124 293 41 13],...
'String','Dip','Style','text');

uicontrol('Parent',h1,'Enable','inactive','Position',[108 252 75 16],...
'String','Depth to Top','Style','text',...
'TooltipString','Depth to the top of the fault (>= 0)',...
'Tag','text14');

uicontrol('Parent',h1,'Enable','inactive','Position',[335 333 41 13],...
'String','Strike','Style','text');

uicontrol('Parent',h1,'Enable','inactive','Position',[404 333 41 13],...
'String','Rake','Style','text');

uicontrol('Parent',h1,'Enable','inactive','Position',[474 333 41 13],...
'String','Slip','Style','text');

uicontrol('Parent',h1, 'Position',[400 105 79 15],...
'Callback',{@deform_okada_uicallback,h1,'radiobutton_deformation_Callback'},...
'String','Deformation',...
'Style','radiobutton',...
'Value',1,...
'Tag','radiobutton_deformation');

uicontrol('Parent',h1, 'Position',[400 85 83 15],...
'Callback',{@deform_okada_uicallback,h1,'radiobutton_interfero_Callback'},...
'String','Interferogram',...
'TooltipString','Will compute an interferogram using the looking vector of the ERS satellite',...
'Style','radiobutton',...
'Tag','radiobutton_interfero');

uicontrol('Parent',h1, 'Position',[485 80 50 51],...
'Background',[1 1 1], ...
'String',{'mm' 'cm' 'm'},...
'Style','listbox',...
'Visible','off',...
'TooltipString','Choose the unites of Slip. A wrong selection leads to a completely erroneous result.',...
'Tag','list_Interfero');

uicontrol('Parent',h1,...
'Callback',{@deform_okada_uicallback,h1,'pushbutton_compute_Callback'},...
'FontWeight','bold',...
'Position',[420 45 71 23],...
'String','Compute',...
'Tag','pushbutton_compute');

uicontrol('Parent',h1,'Enable','inactive','Position',[34 253 41 16],...
'String','Depth','Style','text',...
'TooltipString','Depth to the top of the fault (>= 0)');

uicontrol('Parent',h1,'Enable','inactive','Position',[346 293 21 13],...
'String','u1','Style','text');

uicontrol('Parent',h1,'Enable','inactive','Position',[415 293 21 13],...
'String','u2','Style','text');

uicontrol('Parent',h1,'Enable','inactive','Position',[485 293 21 13],...
'String','u3','Style','text');

uicontrol('Parent',h1,'Enable','inactive','Position',[345 253 21 13],...
'String','Sn','Style','text');

uicontrol('Parent',h1,'Enable','inactive','Position',[414 253 21 13],...
'String','Se','Style','text');

uicontrol('Parent',h1,'Enable','inactive','Position',[484 253 21 13],...
'String','Sz','Style','text');

uicontrol('Parent',h1,...
'Callback',{@deform_okada_uicallback,h1,'pushbutton_cancel_Callback'},...
'FontWeight','bold',...
'Position',[420 15 71 23],...
'String','Cancel',...
'Tag','pushbutton_cancel');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@deform_okada_uicallback,h1,'popup_segment_Callback'},...
'Position',[210 319 91 22],...
'Style','popupmenu',...
'TooltipString','Set parameters with respect to this segment',...
'Value',1,'Tag','popup_segment');

uicontrol('Parent',h1,'Enable','inactive','Position',[225 342 57 15],...
'String','Segment','Style','text','Tag','txtFaultSeg');

uicontrol('Parent',h1,'Enable','inactive','Position',[53 348 85 15],...
'String','Fault Geometry','Style','text','Tag','txtFGeom');

uicontrol('Parent',h1,'Enable','inactive','Position',[373 348 111 15],...
'String','Dislocation Geometry','Style','text','Tag','txtDGeom');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@deform_okada_uicallback,h1,'popup_fault_Callback'},...
'Position',[210 273 91 22],...
'Style','popupmenu',...
'TooltipString','Toggle between faults',...
'Value',1,'Tag','popup_fault');

uicontrol('Parent',h1,'Enable','inactive','Position',[236 295 32 15],...
'String','Faults','Style','text','Tag','fault_number');

uicontrol('Parent',h1,'Enable','inactive','FontSize',10,...
'HorizontalAlignment','left','Position',[400 170 100 16],...
'String','Mw Magnitude =','Style','text','Tag','h_txt_Mw');

uicontrol('Parent',h1,...
'Callback',{@deform_okada_uicallback,h1,'checkbox_hideFaultPlanes_Callback'},...
'Position',[400 140 104 17],...
'String','Hide fault plane',...
'Style','checkbox','Tag','checkbox_hideFaultPlanes');

uicontrol('Parent',h1,'Position',[225 248 55 15],'ForegroundColor',[1 0 0],...
'String','CONFIRM','Style','text');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@deform_okada_uicallback,h1,'popup_GridCoords_Callback'},...
'String', {'Geogs' 'Meters' 'Kilometers'},...
'Position',[210 225 91 22],'Style','popupmenu',...
'TooltipString','GRID COORDINATES: IT IS YOUR RESPONSABILITY THAT THIS IS CORRECT',...
'Value',1,'Tag','popup_GridCoords');

function deform_okada_uicallback(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
feval(callback_name,hObject,[],guidata(h1));

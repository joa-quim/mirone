function varargout = deform_okada(varargin)
% Compute Elastic deformations using the rngchng MEX

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

% $Id: deform_okada.m 7757 2016-01-26 12:28:55Z j $

	if isempty(varargin)
		errordlg('DEFORM OKADA: Wrong number of input args','Error');    return
	end
    
	hObject = figure('Vis','off');
	deform_okada_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hObject,'right');

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
	handles.home_dir = handMir.home_dir;
	handles.work_dir = handMir.work_dir;
	handles.last_dir = handMir.last_dir;
	handles.hCallingFig = handMir.figure1;     % Handles to the calling figure
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
	handles.txt_Mw_pos = get(handles.h_txt_Mw,'Pos');
	handles.Mw(1:handles.n_faults) = 0;
	handles.FaultLength = LineLength(handles.h_fault,handles.geog);
	handles.one_or_zero = ~head(7);
	handles.x_min_or = head(1);			handles.x_max_or = head(2);
	handles.y_min_or = head(3);			handles.y_max_or = head(4);
	handles.mu = 3;							% Shear modulus (x 10^10)

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
			handles = edit_FaultWidth_CB(handles.edit_FaultWidth, handles, faultWidth);    % Compute the rest
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
		uicontrol('Parent',hObject,'Enable','inactive','FontSize',10,'FontName','Helvetica',...
		'HorizontalAlignment','left','Pos',[400 180 100 16],...
		'String',sprintf('Tot Mw = %.1f',mag),'Style','text');

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

	%------------ Give a Pro look (3D) to the frame boxes  --------
	new_frame3D(hObject, [handles.txtFGeom handles.txtDGeom handles.txtGGeom handles.txtIGPos])
	%------------- END Pro look (3D) ------------------------------

	% Make a strike slip image and put it in the focal pushbutton
	ind212 = [1:8 14:19 27:31 40:42 53 54 66 67 79 92 222 235 248 249 261 262 274:276 287:291 300:305 313:320];
	ind255 = [171:181 184:194 197:207 210:220 224:233 237:246 251:259 265:272 279:285 293:298 308:311];
	img = uint8(false(13,25));
	img(ind212) = 212;		img(ind255) = 255;
	img = [img; fliplr(img(12:-1:1,:))];
	img = cat(3, img, img, img);
	set(handles.push_focal,'Cdata',img)

	guidata(hObject, handles);
	set(hObject,'Visible','on');
	if (nargout),   varargout{1} = hObject;     end

% ------------------------------------------------------------------------------------
function handles = edit_FaultWidth_CB(hObject, handles, faultWidth)
% Actualize the "FaultWidth" field. EVENTDATA may not be empty
	if (nargin == 3)
		xx = faultWidth;
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
function edit_FaultStrike_CB(hObject, handles)
	% Cannot be changed

% ------------------------------------------------------------------------------------
function edit_FaultDip_CB(hObject, handles)
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
function edit_FaultDepth_CB(hObject, handles)
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
function edit_FaultTopDepth_CB(hObject, handles)
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
function popup_segment_CB(hObject, handles)
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
function popup_fault_CB(hObject, handles)
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
		set(handles.check_hideFaultPlanes,'Value',1)
	else
		set(handles.check_hideFaultPlanes,'Value',0)
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
		txt = sprintf('Mw Magnitude = %.1f',handles.Mw(fault));
		set(handles.h_txt_Mw,'String',txt,'Pos',handles.txt_Mw_pos + [0 0 30 0])
	else
		set(handles.h_txt_Mw,'String','Mw Magnitude = ','Pos',handles.txt_Mw_pos)
	end
	refresh(handles.hCallingFig);         % otherwise, ML BUG

% ---------------------------------------------------------------
function popup_GridCoords_CB(hObject, handles)
	xx = get(hObject,'Value');
	if (xx == 1),       handles.geog = 1;       handles.is_meters = 0;  handles.is_km = 0;
	elseif (xx == 2),   handles.is_meters = 1;  handles.is_geog = 0;    handles.is_km = 0;
	elseif (xx == 3),   handles.is_km = 1;      handles.is_geog = 0;    handles.is_meters = 0;
	end
	guidata(hObject,handles)

% ------------------------------------------------------------------------------------
function edit_DislocStrike_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	[fault,seg] = getFaultSeg(handles);
	if (isnan(xx)),     set(hObject,'String',handles.DislocStrike{fault}(seg));   return,	end
	handles = convGeometry(handles, fault, seg, 'Aki');
	handles.DislocStrike{fault}(seg) = xx;
	guidata(hObject, handles);

% ------------------------------------------------------------------------------------
function edit_DislocRake_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	[fault,seg] = getFaultSeg(handles);
	if (isnan(xx)),     set(hObject,'String',handles.DislocRake{fault}(seg));   return;     end
	handles = convGeometry(handles, fault, seg, 'Aki');
	handles.DislocRake{fault}(seg) = xx;
	guidata(hObject, handles);

% ------------------------------------------------------------------------------------
function edit_DislocSlip_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	[fault,seg] = getFaultSeg(handles);
	if (isnan(xx)),     set(hObject,'String',handles.DislocSlip{fault}(seg));   return;     end
	handles = convGeometry(handles, fault, seg, 'Aki');
	handles.DislocSlip{fault}(seg) = xx;
	handles = compMag(handles, fault);      % Compute and update Fault's Mw magnitude
	guidata(hObject, handles);

% ------------------------------------------------------------------------------------
function edit_ux_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	[fault,seg] = getFaultSeg(handles);
	if (isnan(xx)),     set(hObject,'String',handles.ux{fault}(seg));   return;     end
	handles = convGeometry(handles, fault, seg, 'Us');
	handles.ux{fault}(seg) = xx;
	guidata(hObject, handles);

% ------------------------------------------------------------------------------------
function edit_uy_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	[fault,seg] = getFaultSeg(handles);
	if (isnan(xx)),     set(hObject,'String',handles.uy{fault}(seg));   return;     end
	handles = convGeometry(handles, fault, seg, 'Us');
	handles.uy{fault}(seg) = xx;
	guidata(hObject, handles);

% ------------------------------------------------------------------------------------
function edit_uz_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	[fault,seg] = getFaultSeg(handles);
	if (isnan(xx)),     set(hObject,'String',handles.uz{fault}(seg));   return;     end
	handles = convGeometry(handles, fault, seg, 'Us');
	handles.uz{fault}(seg) = xx;
	guidata(hObject, handles);

% ------------------------------------------------------------------------------------
function edit_mu_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx)),		set(hObject,'String',handles.mu),	return,		end
	handles.mu = abs(xx);
	fault = getFaultSeg(handles);
	handles = compMag(handles, fault);
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

% --------------------------------------------------------------------
function edit_x_min_CB(hObject, handles)
	dim_funs('xMin', hObject, handles)

% --------------------------------------------------------------------
function edit_x_max_CB(hObject, handles)
	dim_funs('xMax', hObject, handles)

% --------------------------------------------------------------------
function edit_y_min_CB(hObject, handles)
	dim_funs('yMin', hObject, handles)

% --------------------------------------------------------------------
function edit_y_max_CB(hObject, handles)
	dim_funs('yMax', hObject, handles)

% --------------------------------------------------------------------
function edit_x_inc_CB(hObject, handles)
	dim_funs('xInc', hObject, handles)

% --------------------------------------------------------------------
function edit_Ncols_CB(hObject, handles)
	dim_funs('nCols', hObject, handles)

% --------------------------------------------------------------------
function edit_y_inc_CB(hObject, handles)
	dim_funs('yInc', hObject, handles)

% --------------------------------------------------------------------
function edit_Nrows_CB(hObject, handles)
	dim_funs('nRows', hObject, handles)

% ------------------------------------------------------------------------------------
function push_Help_R_CB(hObject, handles)
	message = {'That''s prety obvious to guess what this option does. You select an area,'
		'the grid spacing or the number of rows/columns and the deformation will'
		'be computed at all nodes of that grid.'};
	helpdlg(message,'Help on deformation grid');

% ------------------------------------------------------------------------------------
function push_Help_H_CB(hObject, handles)
	message = {'If you have a file with x,y positions, then the deformation will be computed at those postions'};
	helpdlg(message,'Little Help');

% ------------------------------------------------------------------------------------
function edit_InputFile_CB(hObject, handles, opt)
	if (nargin == 3),		fname = opt;
	else					fname = get(hObject,'String');
	end
	if isempty(fname)
		set(hObject,'Str',''),		handles.input_locations = [];
		guidata(hObject,handles)
		return
	end

	numeric_data = text_read(fname);
	if (isa(numeric_data, 'cell'))
		numeric_data = [numeric_data{:}];
	end
	if (size(numeric_data,2) < 2)
		errordlg('File error. Your file doesn''t have at least 2 columns','Error')
		set(hObject,'Str',''),		handles.input_locations = [];
		guidata(hObject,handles)
		return
	end
	if (~get(handles.check_ToggleXY,'Value'))
		handles.input_locations = numeric_data;
	else
		handles.input_locations = [numeric_data(:,2) numeric_data(:,1)];
	end
	guidata(hObject,handles)

% ------------------------------------------------------------------------------------
function push_InputFile_CB(hObject, handles)
% Just get the file name and sent it to edit_InputFile_CB to do the work
	[FileName,PathName,handles] = put_or_get_file(handles, ...
		{'*.dat;*.DAT;*.xy', 'Maregraph location (*.dat,*.DAT,*.xy)';'*.*', 'All Files (*.*)'},'Select Maregraphs position','get');
	if isequal(FileName,0),		return,		end
	fname = [PathName FileName];
	set(handles.edit_InputFile,'String',fname)
	edit_InputFile_CB(handles.edit_InputFile, handles, fname)

% ------------------------------------------------------------------------------------
function edit_sx_CB(hObject, handles)
	if ( isempty(get(hObject,'String')) ),  set(hObject,'String','0'),	return,		end

% ------------------------------------------------------------------------------------
function edit_sy_CB(hObject, handles)
	if ( isempty(get(hObject,'String')) ),  set(hObject,'String','0'),	return,		end

% ------------------------------------------------------------------------------------
function edit_sz_CB(hObject, handles)
	if ( isempty(get(hObject,'String')) ),  set(hObject,'String','0'),	return,		end

% ------------------------------------------------------------------------------------
function push_focal_CB(hObject, handles)
	strike = str2double(get(handles.edit_DislocStrike,'String'));
	rake = str2double(get(handles.edit_DislocRake,'String'));
	dip = str2double(get(handles.edit_FaultDip,'String'));
	meca_studio(strike, dip, rake)

% ------------------------------------------------------------------------------------
function push_compute_CB(hObject, handles)
	% If cartesian coordinates, they must be in meters
	if (any(isnan(cat(1,handles.FaultWidth{:}))))
		errordlg('One or more segments where not set with the fault''s Width','Error'),		return
	end
	if (any(isnan(cat(1,handles.FaultDepth{:}))))
		errordlg('One or more segments where not set with the fault''s Depth','Error'),		return
	end
	if (any(isnan(cat(1,handles.DislocSlip{:}))))
		errordlg('One or more segments where not set with the movement''s slip','Error'),	return
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
	xmin = str2double(get(handles.edit_x_min,'String'));	xmax = str2double(get(handles.edit_x_max,'String'));
	ymin = str2double(get(handles.edit_y_min,'String'));	ymax = str2double(get(handles.edit_y_max,'String'));
	nrow = str2double(get(handles.edit_Nrows,'String'));	ncol = str2double(get(handles.edit_Ncols,'String'));

	x = handles.fault_x;	y = handles.fault_y;
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

	for (i = 1:handles.n_faults)
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
		if (handles.n_faults > 1)			% We have multiple faults
			U = zeros(nrow,ncol);
			aguentabar(0,'title','Computing deformation','CreateCancelBtn')
			for (i = 1:handles.n_faults)
				U0 = range_change(x{i}(:),y{i}(:),handles.FaultStrike{i}(:),handles.FaultDepth{i}(:),handles.FaultDip{i}(:),...
					handles.ux{i}(:),handles.uy{i}(:),handles.uz{i}(:),handles.FaultLength{i}(:),handles.FaultWidth{i}(:),...
					E,N,s,opt_M);
				U = U0 + U;
				h = aguentabar(i / handles.n_faults);
				if (isnan(h)),	break,	end
			end
			if (isnan(h)),	return,		end
			clear U0;
		else								% We have only one fault
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
		U = single(U);
		if get(handles.radio_deformation,'Value')		% Show deformation 
			mirone(U,head,'Deformation',handles.hCallingFig);
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
			s = [1 0 0; 0 1 0; 0 0 1; s];	% Compute the X,Y,Z components plus the selected projection along unit vector
			U = zeros(numel(E),4);
			for (k = 1:4)					% Loop over the components plus mag
				U_ = 0;
				for (i = 1:handles.n_faults)% Loop over number of faults
					U0 = range_change(x{i}(:),y{i}(:),handles.FaultStrike{i}(:),handles.FaultDepth{i}(:),handles.FaultDip{i}(:),...
						handles.ux{i}(:),handles.uy{i}(:),handles.uz{i}(:),handles.FaultLength{i}(:),handles.FaultWidth{i}(:),...
						E,N,s(k,:),opt_M);
					U_ = U0 + U_;
				end
				U(:,k) = U_;
			end
		else                                % We have only one fault
			s = [1 0 0; 0 1 0; 0 0 1; s];	% Compute the X,Y,Z components plus the selected projection along unit vector
			U = zeros(numel(E),4);
			for (k = 1:4)
				U(:,k) = range_change(x{1}(:),y{1}(:),handles.FaultStrike{1}(:),handles.FaultDepth{1}(:),handles.FaultDip{1}(:),...
					handles.ux{1}(:),handles.uy{1}(:),handles.uz{1}(:),handles.FaultLength{1}(:),handles.FaultWidth{1}(:),...
					E,N,s(k,:),opt_M);
			end
		end

		[FileName,PathName] = put_or_get_file(handles, ...
			{'*.dat;*.DAT', 'Data files (*.dat,*.DAT)';'*.*', 'All Files (*.*)'},'Select output deformation file','put','.dat');
		if isequal(FileName,0),		return,		end
		if (ispc),		fid = fopen([PathName FileName],'wt');
		elseif (isunix)	fid = fopen([PathName FileName],'w');
		else			error('DEFORM_OKADA: Unknown platform.');
		end
		fprintf(fid, '#Lon\t\tLat\t\tNorth\t\tEast\t\tZ\t\tAlong unit vector\n');
		fprintf(fid, '%f\t%f\t%f\t%f\t%f\t%f\n', [E N U]');
		fclose(fid);
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
function radio_deformation_CB(hObject, handles)
	if get(hObject,'Value')
		set(handles.radio_interfero,'Value',0)
		set(handles.edit_sx,'String',num2str(handles.sx))
		set(handles.edit_sy,'String',num2str(handles.sy))
		set(handles.edit_sz,'String',num2str(handles.sz))
		set(handles.list_Interfero,'Vis','off')
	else
		set(hObject,'Value',1)
	end

% ------------------------------------------------------------------------------------
function radio_interfero_CB(hObject, handles)
% Use ERS looking vector
	if get(hObject,'Value')
		set(handles.radio_deformation,'Value',0)
		set(handles.edit_sx,'String','0.333')
		set(handles.edit_sy,'String','-0.07')
		set(handles.edit_sz,'String','0.94')
		set(handles.list_Interfero,'Vis','on')
	else
		set(hObject,'Value',1)
	end

% -----------------------------------------------------------------------------------------
function len = LineLength(h,geog)
	x = get(h,'XData');		y = get(h,'YData');
	D2R = pi/180;			earth_rad = 6371;
	len = [];
	if (~iscell(x))
		if (geog)
			x = x * D2R;	y = y * D2R;
			lat_i = y(1:numel(y)-1);		lat_f = y(2:numel(y));
			lon_i = x(1:numel(x)-1);		lon_f = x(2:numel(x));
			tmp = sin(lat_i).*sin(lat_f) + cos(lat_i).*cos(lat_f).*cos(lon_f-lon_i);
			len = [len; acos(tmp) * earth_rad];			% Distance in km
		else
			dx = diff(x);   dy = diff(y);
			len = [len; sqrt(dx.*dx + dy.*dy)];			% Distance in user unites
		end
	else
		len = cell(1, numel(x));
		if (geog)
			for (k=1:numel(x))
				xx = x{k} * D2R;	yy = y{k} * D2R;
				lat_i = yy(1:length(yy)-1);   lat_f = yy(2:length(yy));
				lon_i = xx(1:length(xx)-1);   lon_f = xx(2:length(xx));
				tmp   = sin(lat_i).*sin(lat_f) + cos(lat_i).*cos(lat_f).*cos(lon_f-lon_i);
				len{k}= acos(tmp) * earth_rad;			% Distance in km
			end
		else
			for (k=1:numel(x))
				xx = x{k};      yy = y{k};
				dx = diff(xx);  dy = diff(yy);
				len{k} = sqrt(dx.*dx + dy.*dy);         % Distance in user unites
			end
		end
	end

% -----------------------------------------------------------------------------------------
function check_hideFaultPlanes_CB(hObject, handles)
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
	handles.hCallingFig = varargin{1}.figure1;
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
	mu = handles.mu * 1e10;
	M0 = mu * handles.um_milhao * handles.DislocSlip{fault}(:) .* handles.FaultWidth{fault}(:) .* ...
		handles.FaultLength{fault}(:);
	if (length(M0) > 1),    M0 = sum(M0);   end
	mag = 2/3*(log10(M0) - 9.1);
	if (~isnan(mag))
		txt = sprintf('Mw Magnitude = %.1f',mag);
		set(handles.h_txt_Mw,'String',txt,'Pos',handles.txt_Mw_pos + [0 0 30 0])
		handles.Mw(fault) = mag;
	end

% ------------------------------------------------------------------------------------
function figure1_CloseRequestFcn(hObject, evt)
	handles = guidata(hObject);
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
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'CloseRequestFcn',@figure1_CloseRequestFcn,...
'MenuBar','none',...
'Name','Okada deformation',...
'NumberTitle','off',...
'Pos',[520 415 536 365],...
'Resize','off',...
'Tag','figure1',...
'HandleVisibility','callback');

uicontrol('Parent',h1,'Pos',[320 224 211 131],'Style','frame');
uicontrol('Parent',h1,'Pos',[10 224 181 131],'Style','frame');
uicontrol('Parent',h1,'Pos',[11 15 350 67],'Style','frame');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Pos',[20 311 71 21],...
'Style','edit',...
'Tooltip','Fault length (km)',...
'Tag','edit_FaultLength');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@deform_okada_uiCB,...
'Pos',[110 311 71 21],...
'Style','edit',...
'Tooltip','Fault width (km)',...
'Tag','edit_FaultWidth');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@deform_okada_uiCB,...
'Pos',[20 271 71 21],...
'Style','edit',...
'Tooltip','Fault strike (degrees)',...
'Tag','edit_FaultStrike');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@deform_okada_uiCB,...
'Pos',[110 271 71 21],...
'Style','edit',...
'Tooltip','Fault dip (degrees)',...
'Tag','edit_FaultDip');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@deform_okada_uiCB,...
'Pos',[20 232 71 21],...
'Style','edit',...
'Tooltip','Depth of the base of fault''s plane',...
'Tag','edit_FaultDepth');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@deform_okada_uiCB,...
'Pos',[110 231 71 21],...
'Style','edit',...
'Tooltip','Alternatively, give depth to the fault''s top ',...
'Tag','edit_FaultTopDepth');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@deform_okada_uiCB,...
'Pos',[330 311 51 21],...
'Style','edit',...
'Tag','edit_DislocStrike');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@deform_okada_uiCB,...
'Pos',[400 311 51 21],...
'Style','edit',...
'Tooltip','Displacement angle clock-wise from horizontal',...
'Tag','edit_DislocRake');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@deform_okada_uiCB,...
'Pos',[470 311 51 21],...
'Style','edit',...
'Tooltip','Total displacement',...
'Tag','edit_DislocSlip');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@deform_okada_uiCB,...
'Pos',[330 271 51 21],...
'Style','edit',...
'Tooltip','Left-lateral displacement along the fault plane (along strike)',...
'Tag','edit_ux');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@deform_okada_uiCB,...
'Pos',[400 271 51 21],...
'Style','edit',...
'Tooltip','Displacement up-dip the fault plane (across strike)',...
'Tag','edit_uy');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@deform_okada_uiCB,...
'Pos',[470 271 51 21],...
'Style','edit',...
'Tooltip','fault tensile slip',...
'Tag','edit_uz');

uicontrol('Parent',h1,'Pos',[10 109 350 93],'Style','frame');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@deform_okada_uiCB,...
'Pos',[330 231 51 21],...
'Style','edit',...
'Tooltip','Component of unit vector along North coords',...
'Tag','edit_sx');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@deform_okada_uiCB,...
'Pos',[400 231 51 21],...
'Style','edit',...
'Tooltip','Component of unit vector along East coords',...
'Tag','edit_sy');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@deform_okada_uiCB,...
'Pos',[470 231 51 21],...
'Style','edit',...
'Tooltip','Component of unit vector along Vertical coords',...
'Tag','edit_sz');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',{@deform_okada_uiCB,h1,'edit_x_min_CB'},...
'Pos',[76 162 71 21],...
'Style','edit',...
'Tooltip','X min value',...
'Tag','edit_x_min');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@deform_okada_uiCB,...
'Pos',[152 162 71 21],...
'Style','edit',...
'Tooltip','X max value',...
'Tag','edit_x_max');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@deform_okada_uiCB,...
'Pos',[228 162 71 21],...
'Style','edit',...
'Tooltip','DX grid spacing',...
'Tag','edit_x_inc');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@deform_okada_uiCB,...
'Pos',[304 162 45 21],...
'Style','edit',...
'Tooltip','Number of columns in the grid',...
'Tag','edit_Ncols');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@deform_okada_uiCB,...
'Pos',[76 136 71 21],...
'Style','edit',...
'Tooltip','Y min value',...
'Tag','edit_y_min');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@deform_okada_uiCB,...
'Pos',[152 136 71 21],...
'Style','edit',...
'Tooltip','Y max value',...
'Tag','edit_y_max');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@deform_okada_uiCB,...
'Pos',[228 136 71 21],...
'Style','edit',...
'Tooltip','DY grid spacing',...
'Tag','edit_y_inc');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@deform_okada_uiCB,...
'Pos',[304 136 45 21],...
'Style','edit',...
'Tooltip','Number of columns in the grid',...
'Tag','edit_Nrows');

uicontrol('Parent',h1,...
'BackgroundColor',[0.831372559070587 0.815686285495758 0.7843137383461],...
'Call',@deform_okada_uiCB,...
'FontWeight','bold',...
'ForegroundColor',[0 0 1],...
'Pos',[289 114 61 18],...
'String','?',...
'Tag','push_Help_R');

uicontrol('Parent',h1,...
'Pos',[22 54 75 19],...
'String','Toggle x,y',...
'Style','checkbox',...
'Tooltip','Toggle x and y columns',...
'Tag','check_ToggleXY');

uicontrol('Parent',h1,'Enable','inactive','Pos',[18 167 55 15],'Str','X Direction','Style','text');
uicontrol('Parent',h1,'Enable','inactive','Pos',[17 141 55 15],'Str','Y Direction','Style','text');
uicontrol('Parent',h1,'Enable','inactive','Pos',[169 184 41 13],'Str','Max','Style','text');
uicontrol('Parent',h1,'Enable','inactive','Pos',[91 185 41 13],'Str','Min','Style','text');
uicontrol('Parent',h1,'Enable','inactive','Pos',[246 185 41 13],'Str','Spacing','Style','text');
uicontrol('Parent',h1,'Enable','inactive','Pos',[302 185 51 13],'Str','# of lines','Style','text');

uicontrol('Parent',h1,'Enable','inactive','Pos',[30 195 121 15],...
'Str','Griding Line Geometry','Style','text','Tag','txtGGeom');

uicontrol('Parent',h1,...
'BackgroundColor',[0.831372559070587 0.815686285495758 0.7843137383461],...
'Call',@deform_okada_uiCB,...
'FontWeight','bold',...
'ForegroundColor',[0 0 1],...
'Pos',[331 50 22 22],...
'String','?',...
'Tag','push_Help_H');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@deform_okada_uiCB,...
'HorizontalAlignment','left',...
'Pos',[21 22 310 22],...
'Style','edit',...
'Tooltip','File name with x, y positions where to compute deformation',...
'Tag','edit_InputFile');

uicontrol('Parent',h1,...
'Call',@deform_okada_uiCB,...
'Pos',[331 21 23 23],...
'Tag','push_InputFile');

uicontrol('Parent',h1,'Enable','inactive','Pos',[32 74 145 15],...
'Str','Input Ground Positions File','Style','text','Tag','txtIGPos');

uicontrol('Parent',h1,'Enable','inactive','Pos',[36 333 41 13],'Str','Length','Style','text');
uicontrol('Parent',h1,'Enable','inactive','Pos',[125 334 41 13],'Str','Width','Style','text');
uicontrol('Parent',h1,'Enable','inactive','Pos',[34 293 41 13],'Str','Strike','Style','text');
uicontrol('Parent',h1,'Enable','inactive','Pos',[124 293 41 13],'Str','Dip','Style','text');

uicontrol('Parent',h1,'Enable','inactive','Pos',[108 252 75 16],...
'Str','Depth to Top','Style','text',...
'Tooltip','Depth to the top of the fault (>= 0)');

uicontrol('Parent',h1,'Enable','inactive','Pos',[335 333 41 13],'Str','Strike','Style','text');
uicontrol('Parent',h1,'Enable','inactive','Pos',[404 333 41 13],'Str','Rake','Style','text');
uicontrol('Parent',h1,'Enable','inactive','Pos',[474 333 41 13],'Str','Slip','Style','text');

uicontrol('Parent',h1, 'Pos',[400 105 90 15],...
'Call',@deform_okada_uiCB,...
'Str','Deformation',...
'Style','radiobutton',...
'Value',1,...
'Tag','radio_deformation');

uicontrol('Parent',h1, 'Pos',[400 85 90 15],...
'Call',@deform_okada_uiCB,...
'String','Interferogram',...
'Tooltip','Will compute an interferogram using the looking vector of the ERS satellite',...
'Style','radiobutton',...
'Tag','radio_interfero');

uicontrol('Parent',h1, 'Pos',[485 80 50 51],...
'Background',[1 1 1], ...
'String',{'mm' 'cm' 'm'},...
'Style','listbox',...
'Visible','off',...
'Tooltip','Choose the unites of Slip. A wrong selection leads to a completely erroneous result.',...
'Tag','list_Interfero');

uicontrol('Parent',h1,'Enable','inactive','Pos',[34 253 41 16],...
'Str','Depth','Style','text',...
'Tooltip','Depth to the top of the fault (>= 0)');

uicontrol('Parent',h1,'Enable','inactive','Pos',[346 293 21 13],'Str','u1','Style','text');
uicontrol('Parent',h1,'Enable','inactive','Pos',[415 293 21 13],'Str','u2','Style','text');
uicontrol('Parent',h1,'Enable','inactive','Pos',[485 293 21 13],'Str','u3','Style','text');
uicontrol('Parent',h1,'Enable','inactive','Pos',[345 253 21 13],'Str','Sn','Style','text');
uicontrol('Parent',h1,'Enable','inactive','Pos',[414 253 21 13],'Str','Se','Style','text');
uicontrol('Parent',h1,'Enable','inactive','Pos',[484 253 21 13],'Str','Sz','Style','text');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@deform_okada_uiCB,...
'Pos',[210 319 91 22],...
'Style','popupmenu',...
'Tooltip','Set parameters with respect to this segment',...
'Value',1,'Tag','popup_segment');

uicontrol('Parent',h1,'Enable','inactive','Pos',[225 342 57 15],'Str','Segment','Style','text','Tag','txtFaultSeg');
uicontrol('Parent',h1,'Enable','inactive','Pos',[53 348 85 15],'Str','Fault Geometry','Style','text','Tag','txtFGeom');
uicontrol('Parent',h1,'Enable','inactive','Pos',[373 348 111 15],'Str','Dislocation Geometry','Style','text','Tag','txtDGeom');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@deform_okada_uiCB,...
'Pos',[210 273 91 22],...
'Style','popupmenu',...
'Tooltip','Toggle between faults',...
'Value',1,'Tag','popup_fault');

uicontrol('Parent',h1,'Enable','inactive','Pos',[236 295 32 15],'Str','Faults','Style','text','Tag','fault_number');
uicontrol('Parent',h1, 'Pos',[419 197 60 18],'Str','Mu (x10^10)','Style','text','Tag','text_mu');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@deform_okada_uiCB,...
'Pos',[481 198 40 21],...
'String','3.0',...
'Style','edit',...
'Tooltip','Shear modulus (for Mw calculation)',...
'Tag','edit_mu');

uicontrol('Parent',h1,'Enable','inactive','FontSize',10,...
'HorizontalAlignment','left','Pos',[400 165 100 16],...
'Str','Mw Magnitude =','Style','text','Tag','h_txt_Mw');

uicontrol('Parent',h1,...
'Call',@deform_okada_uiCB,...
'Pos',[400 140 104 17],...
'Str','Hide fault plane',...
'Style','checkbox','Tag','check_hideFaultPlanes');

uicontrol('Parent',h1,'Pos',[225 248 55 15],'ForegroundColor',[1 0 0],'Str','CONFIRM','Style','text');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@deform_okada_uiCB,...
'String', {'Geogs' 'Meters' 'Kilometers'},...
'Pos',[210 225 91 22],'Style','popupmenu',...
'Tooltip','GRID COORDINATES: IT IS YOUR RESPONSABILITY THAT THIS IS CORRECT',...
'Value',1,'Tag','popup_GridCoords');

uicontrol('Parent',h1,...
'Call',@deform_okada_uiCB,...
'Pos',[375 41 29 29],...
'Tooltip','Show focal mechanism',...
'Tag','push_focal');

uicontrol('Parent',h1,...
'Call',@deform_okada_uiCB,...
'FontWeight','bold',...
'Pos',[420 15 71 21],...
'String','Compute',...
'Tag','push_compute');

function deform_okada_uiCB(hObject, eventdata)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));

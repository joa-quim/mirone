function varargout = deform_mansinha(varargin)
% Compute Elastic deformations

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

	if isempty(varargin)
		errordlg('DEFORM MANSINHA: Wrong number of input args','Error');    return
	end

	hObject = figure('Vis','off');
	deform_mansinha_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hObject,'right');

	if (length(varargin) >= 8)      % Fault patch collection
		handles = set_all_faults(handles,varargin{:});
		handles.fault_in = 1;
	else                            % "Normal" case
		handles.h_fault = varargin{2};			% Handles to the fault lines (each may have more than one segment)
		handles.FaultStrike = varargin{3};
		handles.fault_in = 0;
	end

	handMir = varargin{1};
	handles.last_dir = handMir.last_dir;
	handles.work_dir = handMir.work_dir;
	handles.geog = handMir.geog;
	head = handMir.head;
	handles.head = head;
	handles.hCallingFig = handMir.figure1;		% Handles to the calling figure
	handles.hCallingAxes = handMir.axes1;

	handles.n_faults = numel(handles.h_fault);
	handles.nFaults = handles.n_faults;			% Works as a copy of use with the scc option
	if (handles.n_faults > 1)
		s_format = ['%.' num2str(fix(log10(handles.n_faults))+1) 'd'];
		S = cell(handles.n_faults,1);
		for (i=1:handles.n_faults),     S{i} = ['Fault ' sprintf(s_format,i)];   end
		set(handles.popup_fault,'String',S)
		set(handles.h_fault(1),'LineStyle','--');   % set the top fault one with a dashed line type
		refresh(handMir.figure1);				% otherwise, ML BUG
	else
		set(handles.popup_fault,'Visible','off')
		delete(handles.txtFaultNum)				% Otherwise it would reborn in Pro look
	end

	fault_x = get(handles.h_fault,'XData');     fault_y = get(handles.h_fault,'YData');
	if (handles.n_faults > 1)
		nvert = zeros(1,handles.n_faults);
		for (k = 1:handles.n_faults),  nvert(k) = size(fault_x{k},2) - 1;  end
	else
		nvert = size(fault_x,2) - 1;
	end

	if (any(nvert > 1))
		set(handles.popup_segment,'Visible','on');
		% Even if we have more than one fault, the segments popup will start with only the first fault's segments
		S = cell(nvert(1),1);
		for (i=1:nvert(1)),    S{i} = sprintf('Segment %d',i);   end
		set(handles.popup_segment,'String',S)
	else
		set(handles.popup_segment,'Visible','off')
		delete(handles.txtFaultSeg)    % Otherwise it would reborn in Pro look
	end

	% Try to guess if we are dealing with other (m or km) than geogs
	handles.is_meters = 0;     handles.is_km = 0;   handles.um_milhao = 1e6;
	if (~handles.geog)      % Try to guess if user units are km or meters
		dx = head(2) - head(1);   dy = head(4) - head(3);
		len = sqrt(dx.*dx + dy.*dy);         % Distance in user unites
		if (len > 1e5 || head(8) >= 10)      % If grid's diagonal > 1e5 || Dx >= 10 consider we have meters
			handles.is_meters = 1;     handles.is_km = 0;   handles.um_milhao = 1e3;
			set(handles.popup_GridCoords,'Value',2)
		else
			handles.is_meters = 0;     handles.is_km = 1;
			set(handles.popup_GridCoords,'Value',3)
		end
	end

	handles.fault_x = fault_x;
	handles.fault_y = fault_y;
	handles.nvert = nvert;
	handles.hide_planes(1:handles.n_faults) = 0;
	handles.dms_xinc = 0;           handles.dms_yinc = 0;
	handles.qValue = 0.3;
	handles.patchFatias = [];
	handles.txt_Mw_pos = get(handles.h_txt_Mw,'Position');
	handles.Mw(1:handles.n_faults) = 0;
	handles.FaultLength = LineLength(handles.h_fault,handles.geog);
	handles.one_or_zero = ~head(7);
	handles.x_min_or = head(1);         handles.x_max_or = head(2);
	handles.y_min_or = head(3);         handles.y_max_or = head(4);
	handles.mu = 3;							% Shear modulus (x 10^10)

	if (~iscell(handles.FaultLength)),  handles.FaultLength = {handles.FaultLength};   end

	if (~handles.fault_in)					% "NORMAL" case (not a fault-patch collection)
		% Make them all cell arrays to simplify logic
		if (~iscell(handles.FaultStrike)),  handles.FaultStrike = {handles.FaultStrike};   end
		if (~iscell(handles.fault_x)),      handles.fault_x = {handles.fault_x};    handles.fault_y = {handles.fault_y};   end
		handles.DislocStrike = handles.FaultStrike;

		for (k = 1:handles.n_faults)
			handles.FaultDip{k}(1:nvert(k)) = 25;       handles.FaultWidth{k}(1:nvert(k)) = NaN;
			handles.FaultDepth{k}(1:nvert(k)) = NaN;	handles.FaultTopDepth{k}(1:nvert(k)) = 0;
			handles.DislocSlip{k}(1:nvert(k)) = 1;	    handles.DislocRake{k}(1:nvert(k)) = 90;
		end
		handles.DislocRakeCopy = 90;		% In case of SCC
		handles.DislocSlipCopy = 1;			% 		"

		z2 = sprintf('%.1f',handles.FaultStrike{1}(1));
		set(handles.edit_FaultLength,'String',handles.FaultLength{1}(1))
		set(handles.edit_FaultStrike,'String',z2)
		if (handles.n_faults > 1)				% Otherwise we are allowed to change Fault's length
			set(handles.edit_FaultLength,'Enable','off')
			set(handles.edit_FaultStrike,'Enable','off')
		end
		set(handles.edit_FaultDip,'String',sprintf('%.1f',handles.FaultDip{1}(1)))
		set(handles.edit_DislocStrike,'String',z2)
		set(handles.edit_DislocSlip,'String','1')
		set(handles.edit_DislocRake,'String','90')

		% Default the top depth fault to zero
		set(handles.edit_FaultTopDepth,'String','0')
		if (handles.n_faults == 1)
			faultWidth = handles.FaultLength{1}(1) / 4;
			if (handles.is_meters),		faultWidth = round(faultWidth * 1e-3);     end
			handles = edit_FaultWidth_CB([], handles, faultWidth);    % Compute the rest
            set(handles.edit_FaultWidth,'String',num2str(faultWidth));
		end
	else			% We have a multi-patch slip model
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
		'HorizontalAlignment','left','Position',[460 135 90 16],...
		'String',['Tot = ' sprintf('%.1f',mag)],'Style','text');
	end

	if (handles.n_faults > 1 || ~handles.geog || handles.fault_in)	% In any of these cases there is nothing to save
		set(handles.push_save_subfault,'Vis','off')
	end

	%-----------
	% Fill in the grid limits boxes (in case user wants to compute a grid)
	nDigit = round( log10(abs(max(head(1:4)))) );		% Number of digits of the integer part
	frmt = sprintf('%%.%dg',nDigit+10);					% it will be of the type '%.Ng'
	set(handles.edit_x_min,'String',sprintf(frmt,head(1)))
	set(handles.edit_x_max,'String',sprintf(frmt,head(2)))
	set(handles.edit_y_min,'String',sprintf(frmt,head(3)))
	set(handles.edit_y_max,'String',sprintf(frmt,head(4)))
	handles.x_min = head(1);			handles.x_max = head(2);
	handles.y_min = head(3);			handles.y_max = head(4);
	handles.x_inc = head(8);			handles.y_inc = head(9);

	[m,n] = size(getappdata(handMir.figure1,'dem_z'));

	% Fill in the x,y_inc and nrow,ncol boxes
	nDigit = max(0, round(log10(abs(max(head(8:9))))) );		% Number of digits of the integer part
	frmt = sprintf('%%.%dg',nDigit+10);		% it will be of the type '%.Ng'
	set(handles.edit_Nrows,'String',m);		set(handles.edit_Ncols,'String',n)
	set(handles.edit_y_inc,'String',sprintf(frmt,head(9)))
	set(handles.edit_x_inc,'String',sprintf(frmt,head(8)))

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
	new_frame3D(hObject, [handles.faultGeom handles.dislocGeom handles.gridGeom])
	%------------- END Pro look (3D) ------------------------------

	% Add this figure handle to the carra?as list
	plugedWin = getappdata(handMir.figure1,'dependentFigs');
	plugedWin = [plugedWin hObject];
	setappdata(handMir.figure1,'dependentFigs',plugedWin);

	if (handles.fault_in)
		set([handles.check_SCC handles.edit_qValue handles.edit_nSlices handles.txtqValue handles.txtNfatias],...
			'Visible','off')
	end

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
function edit_FaultStrike_CB(hObject, handles)
% Update the fault's strike to the value entered via this edit box
	xx = str2double(get(hObject,'String'));
	if (isnan(xx))
		set(hObject, 'Str', handles.FaultStrike{1}(1))
		warndlg('Please, pay attention to what you are doing.','Chico Clever')
		return
	end

	[lat2,lon2] = vreckon(handles.fault_y{1}(1), handles.fault_x{1}(1), handles.FaultLength{1}(1)*1e3, xx, 1);
	set(handles.h_fault, 'XData', [handles.fault_x{1}(1) lon2(end)], 'YData', [handles.fault_y{1}(1) lat2(end)])
	handles.fault_x{1}(2) = lon2;		handles.fault_y{1}(2) = lat2;
	handles.FaultStrike{1}(1) = xx;		% We know index is '1' because only single faults can be edited here.
	set(handles.edit_DislocStrike,'String',get(hObject,'Str'))
	edit_FaultWidth_CB(handles.edit_FaultWidth, handles);	% Let this do all the updating work (also uppdates handles)

% ------------------------------------------------------------------------------------
function edit_FaultLength_CB(hObject, handles)
% Update the fault's length to the value entered via this edit box
	xx = str2double(get(hObject,'String'));
	if (isnan(xx) || xx <= 0)
		set(hObject, 'Str', handles.FaultLength{1}(1))
		warndlg('Please, pay attention to what you are doing.','Chico Clever')
		return
	end

	az = str2double(get(handles.edit_FaultStrike,'String'));
	[lat2,lon2] = vreckon(handles.fault_y{1}(1), handles.fault_x{1}(1), xx * 1000, az, 1);
	set(handles.h_fault, 'XData', [handles.fault_x{1}(1) lon2(end)], 'YData', [handles.fault_y{1}(1) lat2(end)])
	handles.fault_x{1}(2) = lon2;		handles.fault_y{1}(2) = lat2;
	handles.FaultLength{1}(1) = xx;		% We know index is '1' because only single faults can be edited here.
	edit_FaultWidth_CB(handles.edit_FaultWidth, handles);	% Let this do all the updating work (also uppdates handles)

% ------------------------------------------------------------------------------------
function handles = edit_FaultWidth_CB(hObject, handles, opt)
% Actualize the "FaultWidth" field. EVENTDATA may not be empty
	if (nargout),		xx = opt;
	else				xx = str2double(get(hObject,'String'));
	end
	if (xx < 0)         % If user tried to give a negative width
        xx = -xx;       set(hObject,'String',xx)
	end
	dip = str2double(get(handles.edit_FaultDip,'String'));
	top_d = str2double(get(handles.edit_FaultTopDepth,'String'));
	depth = top_d + xx * cos((90-dip)*pi/180);
	set(handles.edit_FaultDepth,'String',depth);
	[fault,seg] = getFaultSeg(handles);
	handles.FaultWidth{fault}(seg) = xx;
	handles.FaultDepth{fault}(seg) = depth;
	handles.FaultWidthCopy = xx;        % To use (if) in scc. It is a scalar because one-fault-only
	
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
	try     set(hp(seg), 'XData',x, 'YData',y, 'FaceColor',[.8 .8 .8],'EdgeColor','k','LineWidth',1);  end
	
	z = -[top_d top_d depth depth top_d];
	if ( diff(handles.head(5:6)) > 10 ),	z = z * 1000;		end		% Assume grid's depth is in meters
	z = z + handles.head(5);				% CRUDE. It should be mean depth along the fault's length
	set(hp, 'UserData', z)					% So that we can Flederize it in 3D 

	handles = compMag(handles, fault);      % Compute and update Fault's Mw magnitude
	guidata(handles.figure1,handles)

% ------------------------------------------------------------------------------------
function edit_FaultDip_CB(hObject, handles)
% Actualize the "FaultDip" field
	xx = str2double(get(hObject,'String'));
	top_d = str2double(get(handles.edit_FaultTopDepth,'String'));
	W = str2double(get(handles.edit_FaultWidth,'String'));
	depth = top_d + W * cos((90-xx)*pi/180);
	set(handles.edit_FaultDepth,'String',depth);
	[fault,seg] = getFaultSeg(handles);
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
	if (diff(handles.head(5:6)) > 10),		z = z * 1000;		end		% Assume grid's depth is in meters
	z = z + handles.head(5);				% CRUDE. It should be mean depth along the fault's length
	set(hp, 'UserData', z)					% So that we can Flederize it in 3D 

	guidata(handles.figure1,handles)

% ------------------------------------------------------------------------------------
function edit_FaultDepth_CB(hObject, handles)
% Actualize the "FaultTopDepth" field
	xx = str2double(get(hObject,'String'));
	if (xx < 0)         % If user tried to give a negative depth
        xx = -xx;       set(hObject,'String',xx)
	end
	W = str2double(get(handles.edit_FaultWidth,'String'));
	dip = str2double(get(handles.edit_FaultDip,'String'));
	top_d = xx - W * cos((90-dip)*pi/180);
	set(handles.edit_FaultTopDepth,'String',top_d);
	[fault,seg] = getFaultSeg(handles);
	handles.FaultDepth{fault}(seg) = xx;
	handles.FaultTopDepth{fault}(seg) = top_d;
	guidata(handles.figure1,handles)

% ------------------------------------------------------------------------------------
function edit_FaultTopDepth_CB(hObject, handles)
% Actualize the "FaultDepth" field
	xx = str2double(get(hObject,'String'));
	if (xx < 0)         % If user tried to give a negative depth
        xx = -xx;       set(hObject,'String',xx)
	end
	W = str2double(get(handles.edit_FaultWidth,'String'));
	dip = str2double(get(handles.edit_FaultDip,'String'));
	depth = xx + W * cos((90-dip)*pi/180);
	set(handles.edit_FaultDepth,'String',depth);
	[fault,seg] = getFaultSeg(handles);
	handles.FaultTopDepth{fault}(seg) = xx;
	handles.FaultDepth{fault}(seg) = depth;
	guidata(handles.figure1,handles)

% ------------------------------------------------------------------------------------
function popup_segment_CB(hObject, handles)
	seg = get(hObject,'Value');
	fault = getFaultSeg(handles);

	% Fault parameters
	set(handles.edit_FaultLength,'String',handles.FaultLength{fault}(seg))
	set(handles.edit_FaultStrike,'String',sprintf('%.1f',handles.FaultStrike{fault}(seg)))

	if (isnan(handles.FaultWidth{fault}(seg))),    str = '';
	else	str = num2str(handles.FaultWidth{fault}(seg));
	end
	set(handles.edit_FaultWidth,'String',str)

	set(handles.edit_FaultDip,'String',sprintf('%.1f',handles.FaultDip{fault}(seg)))
	set(handles.edit_FaultTopDepth,'String',num2str(handles.FaultTopDepth{fault}(seg)))

	if (isnan(handles.FaultDepth{fault}(seg))),    str = '';
	else	str = num2str(handles.FaultDepth{fault}(seg));
	end
	set(handles.edit_FaultDepth,'String',str)

	% Dislocation parameters
	set(handles.edit_DislocStrike,'String',sprintf('%.1f',handles.DislocStrike{fault}(seg)))
	if (isnan(handles.DislocSlip{fault}(seg))),    str = '';
	else	str = num2str(handles.DislocSlip{fault}(seg));
	end
	set(handles.edit_DislocSlip,'String',str)
	if (isnan(handles.DislocRake{fault}(seg))),    str = '';
	else	str = sprintf('%.1f',handles.DislocRake{fault}(seg));
	end
	set(handles.edit_DislocRake,'String',str)

% -----------------------------------------------------------------------------------------
function popup_fault_CB(hObject, handles)
	fault = get(hObject,'Value');
	S = cell(handles.nvert(fault),1);
	for (i=1:handles.nvert(fault)),    S{i} = ['Segment ' sprintf('%d',i)];   end
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
	set(handles.edit_FaultLength,'String',handles.FaultLength{fault}(seg))
	set(handles.edit_FaultStrike,'String',sprintf('%.1f',handles.FaultStrike{fault}(seg)))

	if (isnan(handles.FaultWidth{fault}(seg))),    str = '';
	else	str = num2str(handles.FaultWidth{fault}(seg));
	end
	set(handles.edit_FaultWidth,'String',str)

	set(handles.edit_FaultDip,'String',sprintf('%.1f',handles.FaultDip{fault}(seg)))
	set(handles.edit_FaultTopDepth,'String',handles.FaultTopDepth{fault}(seg))

	if (isnan(handles.FaultDepth{fault}(seg))),    str = '';
	else	str = num2str(handles.FaultDepth{fault}(seg));
	end
	set(handles.edit_FaultDepth,'String',str)

	% Dislocation parameters
	set(handles.edit_DislocStrike,'String',sprintf('%.1f',handles.DislocStrike{fault}(seg)))
	if (isnan(handles.DislocSlip{fault}(seg))),    str = '';
	else	str = num2str(handles.DislocSlip{fault}(seg));
	end
	set(handles.edit_DislocSlip,'String',str)

	if (isnan(handles.DislocRake{fault}(seg))),    str = '';
	else    str = sprintf('%.1f',handles.DislocRake{fault}(seg));
	end
	set(handles.edit_DislocRake,'String',str)
	if (handles.Mw(fault) > 0)
        txt = ['Mw Magnitude = ' sprintf('%.1f',handles.Mw(fault))];
        set(handles.h_txt_Mw,'String',txt,'Position',handles.txt_Mw_pos + [0 0 30 0])
	else
        set(handles.h_txt_Mw,'String','Mw Magnitude = ','Position',handles.txt_Mw_pos)
	end
	refresh(handles.hCallingFig);         % otherwise, ML BUG

% ------------------------------------------------------------------------------------
function edit_DislocStrike_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	[fault,seg] = getFaultSeg(handles);
	if (isnan(xx)),     set(hObject,'String',handles.DislocStrike{fault}(seg));   return;     end
	handles.DislocStrike{fault}(seg) = xx;
	guidata(handles.figure1,handles)

% ------------------------------------------------------------------------------------
function edit_DislocRake_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	[fault,seg] = getFaultSeg(handles);
	if (isnan(xx)),     set(hObject,'String',handles.DislocRake{fault}(seg));   return;     end
	handles.DislocRake{fault}(seg) = xx;
	handles.DislocRakeCopy = xx;            % To be of use in the SCC cases
	if (get(handles.check_SCC,'Val'))
        handles.DislocRake = repmat({xx},1,handles.n_faults);
	end
	guidata(handles.figure1,handles)

% ------------------------------------------------------------------------------------
function edit_DislocSlip_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	[fault,seg] = getFaultSeg(handles);
	if (isnan(xx))
		set(hObject,'String','')
		handles.DislocSlip{fault}(seg) = NaN;
		return
	else
		handles.DislocSlip{fault}(seg) = xx;
	end
	handles.DislocSlipCopy = xx;            % To be of use in the SCC cases

	handles = compMag(handles, fault);      % Compute and update Fault's Mw magnitude
	guidata(handles.figure1,handles)

% -------------------------------------------------------------------------------------
function edit_x_min_CB(hObject, handles)
	dim_funs('xMin', hObject, handles)

% -------------------------------------------------------------------------------------
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
function push_focal_CB(hObject, handles)
	strike = str2double(get(handles.edit_DislocStrike,'String'));
	rake = str2double(get(handles.edit_DislocRake,'String'));
	dip = str2double(get(handles.edit_FaultDip,'String'));
	meca_studio(strike, dip, rake)

% ------------------------------------------------------------------------------------
function push_compute_CB(hObject, handles)
	% If cartesian coordinates, they must be in meters
	if (any(isnan(cat(1,handles.FaultWidth{:}))))
		errordlg('One or more segments where not set with the fault''s Width','Error');    return
	end
	if (any(isnan(cat(1,handles.FaultDepth{:}))))
		errordlg('One or more segments where not set with the fault''s Depth','Error');    return
	end
	if (any(isnan(cat(1,handles.DislocRake{:}))))
		errordlg('One or more segments where not set with the movement''s rake','Error');    return
	end
	if (any(isnan(cat(1,handles.DislocSlip{:}))))
		errordlg('One or more segments where not set with the movement''s slip','Error');    return
	end

	% Get grid params
	xmin = str2double(get(handles.edit_x_min,'String'));
	xmax = str2double(get(handles.edit_x_max,'String'));
	ymin = str2double(get(handles.edit_y_min,'String'));
	ymax = str2double(get(handles.edit_y_max,'String'));
	xinc = str2double(get(handles.edit_x_inc,'String'));
	yinc = str2double(get(handles.edit_y_inc,'String'));
	nrow = str2double(get(handles.edit_Nrows,'String'));
	ncol = str2double(get(handles.edit_Ncols,'String'));

	if (handles.is_km)      % Than we must convert those to meters
		xmin = xmin * 1e3;    xmax = xmax * 1e3;
		ymin = ymin * 1e3;    ymax = ymax * 1e3;
		xinc = xinc * 1e3;    yinc = yinc * 1e3;
	end

	opt_R = sprintf('-R%.12g/%.12g/%.12g/%.12g',xmin, xmax, ymin, ymax);
	opt_I = sprintf('-I%.12g/%.12g',xinc, yinc);

	n_seg = sum(handles.nvert);     % Total number of fault segments
	x = handles.fault_x;            y = handles.fault_y;
	if (~iscell(x)),                x = {x};    y = {y};    end
	kk = 1;
	opt_F = cell(1, n_seg * handles.n_faults);		opt_A = opt_F;		opt_E = opt_F;
	for (i = 1:handles.n_faults)
		for (k = 1:handles.nvert(i))		% Loop over number of segments of this fault
			if (handles.is_meters)			% Fault's length must be given in km to mansinha_m
				handles.FaultLength{i}(k) = handles.FaultLength{i}(k) / 1000;
			elseif (handles.is_km)			% This is a messy case. -E & -I must also be in meters
				x{i}(k) = x{i}(k) * 1e3;    y{i}(k) = y{i}(k) * 1e3;
			end
			opt_F{kk} = ['-F' num2str(handles.FaultLength{i}(k)) '/' num2str(handles.FaultWidth{i}(k)) '/' ...
					num2str(handles.FaultTopDepth{i}(k))];
			dip = handles.FaultDip{i}(k);
			if (dip == 90.0),	dip = 89.9999;	end		% There is bug in mex that screws when dip = 90.0
			opt_A{kk} = ['-A' num2str(dip) '/' num2str(handles.FaultStrike{i}(k)) '/' ...
					num2str(handles.DislocRake{i}(k)) '/' num2str(handles.DislocSlip{i}(k))];
			opt_E{kk} = sprintf('-E%.5f/%.5f',x{i}(k), y{i}(k));
			kk = kk + 1;
		end
	end

	if (handles.geog),  opt_M = '-M';
	else                opt_M = '';
	end

	% Compute deformation
	if (n_seg > 1)
		U = zeros(nrow,ncol);
		aguentabar(0,'title','Computing deformation','CreateCancelBtn')
		for (k = 1:n_seg)
			U0 = double(mansinha_m(opt_R, opt_I, opt_A{k}, opt_F{k}, opt_E{k}, opt_M));
			U = U0 + U;
			h = aguentabar(k / n_seg);
			if (isnan(h)),	break,	end
		end
		if (isnan(h)),	return,		end
		clear U0;
	else
		U = mansinha_m(opt_R, opt_I, opt_A{1}, opt_F{1}, opt_E{1}, opt_M);
	end

	z_max = double(max(U(:)));		z_min = double(min(U(:)));
	dx = str2double(get(handles.edit_x_inc,'String'));
	dy = str2double(get(handles.edit_y_inc,'String'));
    
    [m,n] = size(U);		hWarn = [];
    if ( (m ~= nrow) || (n ~= ncol) )
        msg{1} = 'Someting went wrong. Output file has not the required size. Maybe a meters<->kilometers bad guess?';
        if (abs(dx - dy) > 1e-5)
            msg{2} = ' ';
            msg{3} = 'No. Almost likely this was due to the fact that or X and Y spacings are diferent.';
        end
        hWarn = warndlg(msg,'Warning');
    end

	head.head = [xmin xmax ymin ymax z_min z_max 0 dx dy];
	head.X = linspace(xmin,xmax,ncol);
	head.Y = linspace(ymin,ymax,nrow);

    U = single(U);
	mirone(U,head,'Deformation',handles.hCallingFig);
    if (~isempty(hWarn)),   figure(hWarn);      end

% -----------------------------------------------------------------------------------------
function len = LineLength(h,geog)
	x = get(h,'XData');     y = get(h,'YData');
	len = [];
	if (~iscell(x))
		if (geog)
			ll = draw_funs(h, 'show_LineLength', [], [], h, 'k');
			len = ll.seg_len;
		else
			dx = diff(x);   dy = diff(y);
			len = [len; sqrt(dx.*dx + dy.*dy)];         % Distance in user unites
		end
	else
		len = cell(1, numel(x));
		if (geog)
			for (k = 1:numel(x))
				ll = draw_funs(h(k), 'show_LineLength', [], [], h(k), 'k');
				len{k} = ll.seg_len;
			end
		else
			for (k = 1:numel(x))
				xx = x{k};      yy = y{k};
				dx = diff(xx);  dy = diff(yy);
				len{k} = sqrt(dx.*dx + dy.*dy);         % Distance in user unites
			end
		end
	end

% -----------------------------------------------------------------------------------------
function popup_GridCoords_CB(hObject, handles)
	xx = get(hObject,'Value');
	if (xx == 1),       handles.geog = 1;       handles.is_meters = 0;  handles.is_km = 0;  handles.um_milhao = 1e6;
	elseif (xx == 2),   handles.is_meters = 1;  handles.is_geog = 0;    handles.is_km = 0;  handles.um_milhao = 1e3;
	elseif (xx == 3),   handles.is_km = 1;      handles.is_geog = 0;    handles.is_meters = 0;  handles.um_milhao = 1e6;
	end
    guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------------------
function check_hideFaultPlanes_CB(hObject, handles)
	fault = getFaultSeg(handles);
	hp = getappdata(handles.h_fault(fault),'PatchHand');
	if (get(hObject,'Value'))
        try    set(hp,'Visible','off');    end
        handles.hide_planes(fault) = 1;
	else
        try    set(hp,'Visible','on');     end
        handles.hide_planes(fault) = 0;
	end
    guidata(handles.figure1,handles)

% ------------------------------------------------------------------------------------
function edit_mu_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx)),		set(hObject,'String',handles.mu),	return,		end
	handles.mu = abs(xx);
	fault = getFaultSeg(handles);
	handles = compMag(handles, fault);
	guidata(hObject, handles);

% -----------------------------------------------------------------------------------------
function check_SCC_CB(hObject, handles)
% Activate (or de-...) the variable slip mode.
	if (~handles.geog)
		warndlg('Sorry but this works only with geographic grids','Warning')
		set(hObject,'Value', 0)
		return
	end

	if (get(hObject,'Value'))
		msg = [];
		if (handles.nFaults ~= 1),              msg = 'Currently only one fault is allowed';    end
		if (numel(handles.FaultWidth{1}) ~= 1), msg = 'Currently only a fault with one segment is allowed';    end
		if (~isempty(msg)),     errordlg(msg,'Error');      set(hObject,'Value',0),		return,		end
		if (isnan(handles.FaultDepth{1})),      msg = 'Must set the fault''s Depth first';          end
		if (isnan(handles.FaultWidth{1})),      msg = 'Must set the total fault''s Width first';    end
		if (~isempty(msg))
			errordlg(msg,'Error');
			set(hObject,'Value',0);
			return
		end
		handles.restoreOldPlane = ~handles.hide_planes;
		handles.nvert_back = handles.nvert;		% We will need the original if unsetting the SCC
		do_scc(handles);
	else								% Remove the scc stripes and restore previous const slip state
		deleteFatias([],[], handles)
		if (handles.restoreOldPlane)
			hp = getappdata(handles.h_fault(1),'PatchHand');
			try     set(hp,'Visible','on');     end
			set(handles.check_hideFaultPlanes,'Val',0)
			handles.hide_planes(1) = 0;
		else
			set(handles.check_hideFaultPlanes,'Val',1)
		end

		handles.patchFatias = [];       % Reset it so a new cicle may begin
		if (isappdata(handles.figure1,'handBak'))
			handBak = getappdata(handles.figure1,'handBak');
			handles.n_faults = handBak.n_faults;        handles.FaultTopDepth = handBak.FaultTopDepth;
			handles.nvert = handBak.nvert;              handles.FaultLength = handBak.FaultLength;
			handles.FaultDip = handBak.FaultDip;        handles.FaultWidth = handBak.FaultWidth;
			handles.DislocRake = handBak.DislocRake;    handles.DislocSlip = handBak.DislocSlip;
			handles.fault_x = handBak.fault_x;          handles.fault_y = handBak.fault_y;
			handles.FaultStrike = handBak.FaultStrike(1);
		end
		
		guidata(handles.figure1,handles)
	end

% -----------------------------------------------------------------------------------------
function edit_nSlices_CB(hObject, handles)
    nFatias = str2double(get(hObject,'String'));
    if (isnan(nFatias) || nFatias < 5)
        warndlg('Number of slices requested is either nonsense or too low. Reseting','Warning')
        set(handles.edit_nSlices,'String','20')
    end
    if (get(handles.check_SCC,'Val'))	% We are already on the Fatias building mode. So rebuild them
        do_scc(handles);
    end

% -----------------------------------------------------------------------------------------
function edit_qValue_CB(hObject, handles)
    qValue = str2double(get(hObject,'String'));
    if (isnan(qValue) || qValue < 1e-2 || qValue >= 1)
        set(hObject,'String','0.3');
        qValue = 0.3;
    end
    handles.qValue = qValue;
    if (get(handles.check_SCC,'Val'))	% We are already on the Fatias building mode. So rebuild them
        do_scc(handles);				% It also saves handles
    end

% ------------------------------------------------------------------------------------
function push_save_subfault_CB(hObject, handles)
% Save the current solution in the sub-fault format.
% This function, however, may only be called when geogs and single fault, single segment

	str = {'*.dat;*.DAT;*.txt;*.TXT', 'Data file (*.dat,*.DAT,*.txt,*.TXT)'; '*.*', 'All Files (*.*)'};
	[FileName,PathName] = put_or_get_file(handles,str,'Select file','put','.dat');
	if isequal(FileName,0),		return,		end
	f_name = [PathName FileName];

	fid = fopen(f_name, 'wt');
	if (fid < 0)
		errordlg('Could not open output file.','Error'),	return
	end
	D2R = pi / 180;
	rng = handles.FaultWidth{1}(1) / 6371 / D2R;
	strk = handles.FaultStrike{1}(1);

	% Write header
	fprintf(fid, '#Total number of fault_segments=     1\n');
	fprintf(fid, '#Fault_segment =   1 nx(Along-strike)=   1 Dx= %.2fkm ny(downdip)=   1 Dy= %.2fkm\n', ...
		handles.FaultLength{1}, handles.FaultWidth{1});
	% Write BB
	fprintf(fid, '#Boundary of Fault_segment     1\n');
	fprintf(fid, '#Lon.  Lat.  Depth\n');
	fprintf(fid, '%.5f\t%.5f\t%.5f\n', handles.fault_x{1}(1), handles.fault_y{1}(1), handles.FaultTopDepth{1}(1));
	fprintf(fid, '%.5f\t%.5f\t%.5f\n', handles.fault_x{1}(2), handles.fault_y{1}(2), handles.FaultTopDepth{1}(1));
	[lat1,lon1] = circ_geo(handles.fault_y{1}(1), handles.fault_x{1}(1), rng, strk+90, 1);
	[lat2,lon2] = circ_geo(handles.fault_y{1}(2), handles.fault_x{1}(2), rng, strk+90, 1);
	fprintf(fid, '%.5f\t%.5f\t%.5f\n', lon2, lat2, handles.FaultDepth{1}(1));
	fprintf(fid, '%.5f\t%.5f\t%.5f\n', lon1, lat1, handles.FaultDepth{1}(1));
	fprintf(fid, '%.5f\t%.5f\t%.5f\n', handles.fault_x{1}(1), handles.fault_y{1}(1), handles.FaultTopDepth{1}(1));
	% Write the patch (single one for time being)
	fprintf(fid, '#Lat. Lon. depth slip rake strike dip\n');

	% The shit here is that we need the coordinates of the MIDDLE of the fault trace
	% So we will the same thing as will be done in fault_models/subfault(), but in reverse order.
	rng = (handles.FaultLength{1}(1) / 2) / 6371 / D2R;
	[lat1,lon1] = circ_geo(handles.fault_y{1}(1), handles.fault_x{1}(1), rng, strk, 1);
	fprintf(fid, '%.4f\t%.4f\t%.3f\t%.2f\t%.1f\t%.1f\t%.1f\n', lat1, lon1, handles.FaultDepth{1}(1), ...
		handles.DislocSlip{1}(1)*100, handles.DislocRake{1}(1), handles.DislocStrike{1}(1), handles.FaultDip{1}(1));
	fclose(fid);

% -----------------------------------------------------------------------------------------
function stripes = do_scc(handles)
% Take a constant slip fault plane and divide it it into slices of variable slip
    
    % This is a trick to use the callabck to hide the const slip patches
    if (~get(handles.check_hideFaultPlanes,'Val'))
        set(handles.check_hideFaultPlanes,'Val',1)
        check_hideFaultPlanes_CB(handles.check_hideFaultPlanes, handles)
        set(handles.check_hideFaultPlanes,'Val',0)
    end
    
    % Backup those. We may want to revert to these values
    if (~isappdata(handles.figure1,'handBak'))
        handBak.n_faults = handles.n_faults;        handBak.FaultTopDepth = handles.FaultTopDepth;
        handBak.nvert = handles.nvert;              handBak.FaultLength = handles.FaultLength;
        handBak.FaultDip = handles.FaultDip;        handBak.FaultWidth = handles.FaultWidth;
        handBak.DislocRake = handles.DislocRake;    handBak.DislocSlip = handles.DislocSlip;
        handBak.fault_x = handles.fault_x;          handBak.fault_y = handles.fault_y;
        handBak.FaultStrike = handles.FaultStrike(1);
        setappdata(handles.figure1,'handBak',handBak)
    end

    if (~isempty(handles.patchFatias))          % The SCC slices already exist. We are recomputing. Take care
        handBak = getappdata(handles.figure1,'handBak');
        handles.FaultTopDepth = handBak.FaultTopDepth;   handles.FaultWidth = handBak.FaultWidth;
        handles.FaultDip = handBak.FaultDip;             handles.FaultStrike = handBak.FaultStrike;
        handles.fault_x = handBak.fault_x;               handles.fault_y = handBak.fault_y;
        delete(handles.patchFatias)             % Delete existing slices patches before drawing new ones below
    end
    [stripes,varSlip,sliceWidth] = comp_varSlip(handles);   % <== Compute the variable slip per slice of the source
    handles.varSlip = varSlip;                  % e.g. for ploting
    handles.n_faults = numel(stripes);
    handles.FaultTopDepth = cell(handles.n_faults,1);
    handles.patchFatias = zeros(1,handles.n_faults);
    stripeVis = 'on';           % Control whether stripes patches are drawn visible or ... not
    if (get(handles.check_hideFaultPlanes,'Val')),        stripeVis = 'off';    end
    cmenuHand = zeros(1,handles.n_faults);

    for (i=1:handles.n_faults)
        handles.patchFatias(i) = patch('XData',stripes{i}(1,:),'YData',stripes{i}(2,:),...
            'Parent',handles.hCallingAxes,'FaceColor',rand(1,3),'Visible',stripeVis,'Tag','Fatia');
        handles.FaultTopDepth{i} = stripes{i}(3,1:2);
        cmenuHand(i) = uicontextmenu('Parent',handles.hCallingFig);
        set(handles.patchFatias(i), 'UIContextMenu', cmenuHand(i))
        uimenu(cmenuHand(i), 'Label', 'Delete', 'Call',{@deleteFatias,handles});
        uimenu(cmenuHand(i), 'Label', 'Show slip profile', 'Call',{@showFatias,handles});
    end
    
    handles.nvert = ones(1,handles.n_faults);       % Though 'nvert' it means number of segments of each fault, which is 1 here
    handles.FaultLength = repmat(handles.FaultLength(1),1,handles.n_faults);  % Replicate the (cell array) fault length
    handles.FaultDip = repmat(handles.FaultDip(1),1,handles.n_faults);
    handles.FaultStrike = repmat(handles.FaultStrike(1),1,handles.n_faults);
    handles.FaultWidth = repmat({sliceWidth},1,handles.n_faults);
    
    handles.fault_x = cell(handles.n_faults,1);    handles.fault_y = handles.fault_x;
    for (i=1:handles.n_faults)
        handles.fault_x{i} = stripes{i}(1,1:2);
        handles.fault_y{i} = stripes{i}(2,1:2);
    end
    
    if (~isnan(handles.DislocRake{1}))              % We'll allow this to be set after
        handles.DislocRake = repmat({handles.DislocRakeCopy},1,handles.n_faults);
    end
    if (~isnan(handles.DislocSlip{1}))              % Idem
        handles.DislocSlip = num2cell(handles.DislocSlipCopy * 2 * varSlip);
    end
    
    guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------------------
function [stripes,varSlip,sliceWidth] = comp_varSlip(handles)
% Compute the descretized slip according to eq 10 of Geist & Demowska, PAGEOPH, 485-512, 1999
    nFatias = round(str2double(get(handles.edit_nSlices,'String')));
    q = handles.qValue;
    D2R = pi / 180;
    gama  = (0:nFatias-1) / nFatias;
    D1 = 12/q^3*(q*gama - gama.^2);                 % Two branches of the variable slip equation
    D2 = 12/(1-q)^3*(gama.^2 -gama*(1+q)+q);
    D3 = [D1(gama < q) D2(gama >= q)];              % Join them
    varSlip = cumsum(D3);                           % Integrate to get the slip
    varSlip(varSlip < 0) = 0;                       % At the end we normally have 1 or 2 negative values
    varSlip = varSlip / max(varSlip);               % Normalize
    
    fac = 1 / 6371 * 180 / pi;
    sliceWidth = handles.FaultWidthCopy / nFatias;
    rng = sliceWidth * cos(handles.FaultDip{1} * D2R) * fac;
    stripes = cell(nFatias,1);
    lon_i = handles.fault_x{1}(1);          lon_f = handles.fault_x{1}(2);
    lat_i = handles.fault_y{1}(1);          lat_f = handles.fault_y{1}(2);
    z_i   = handles.FaultTopDepth{1}(1);    z_f = z_i;
    for (i=1:nFatias)
        [lat_a,lon_a] = circ_geo(lat_i,lon_i,rng,handles.FaultStrike{1}+90,1);
        [lat_b,lon_b] = circ_geo(lat_f,lon_f,rng,handles.FaultStrike{1}+90,1);
        z_a = z_i + sliceWidth * sin(handles.FaultDip{1} * D2R);    z_b = z_a;
        stripes{i} = [lon_i lon_f lon_b lon_a; lat_i lat_f lat_b lat_a; z_i z_f z_b z_a];    % 4 corners stripe coordinates
        lon_i = lon_a;      lon_f = lon_b;        % a-b side of this stripe will be i-f of the next
        lat_i = lat_a;      lat_f = lat_b;
        z_i = z_a;          z_f = z_b;
    end

    %   i               f
    %   |---------------|
    %   |               |           stripe{i} has the coordinates in this order i->f->b->a
    %   |---------------|
    %   a               b

% -----------------------------------------------------------------------------------------
function deleteFatias(hObj, evt, handles)
	handles = guidata(handles.figure1);     % Update
	delete(handles.patchFatias)
	handles.patchFatias = [];
	guidata(handles.figure1, handles);

% -----------------------------------------------------------------------------------------
function showFatias(hObj, evt, handles)
	q = handles.qValue;
	gama = (0:100)/100;
	D1 = 12/q^3*(q*gama - gama.^2);     D2 = 12/(1-q)^3*(gama.^2 -gama*(1+q)+q);
	ind = find(gama < q);
	D3 = [D1(1:ind(end)) D2(ind(end)+1:end)];
	D3 = cumsum(D3);                    D3(D3 < 0) = 0;
	D3 = D3 / max(D3);
	hFig = figure('Name','Normalized slip profile','NumberTitle','off');    
	plot(gama,D3)
	hAxes = get(hFig,'Currentaxes');

	nFatias = numel(handles.patchFatias);
	gama  = (0:nFatias-1)/nFatias;
	D3 = [handles.varSlip(1:end); handles.varSlip(1:end)];
	D3 = D3(:)';
	x = [gama(2:end); gama(2:end)];
	x = [0 x(:)' gama(end)+diff(gama(1:2))];
	line('XData',x,'YData',D3,'Parent',hAxes)
	area = cumsum(handles.varSlip) * diff(gama(1:2));
	lost = (1 - area(end) / 0.5)*100;
	text(.25,.1,sprintf('%s%.1f%s','Lost in descretization = ',abs(lost),'%'))

% ------------------------------------------------------------------------------------
function [fault,seg] = getFaultSeg(handles)
	fault = 1;
	if (handles.n_faults > 1),	fault = get(handles.popup_fault,'Value');	end
	seg = get(handles.popup_segment,'Value');

% -------------------------------------------------------------------------------------
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

		% Pre-allocation
		handles.FaultDepth = repmat({zeros(1, nx)}, n_fault, 1);
		
	else				% Multi-seg Slip model file. Old logic obliges to have {1,nFaultTotal}
		% Remember that multi-segs had declaration in fault_models.m like:
		% hLine = zeros(nSeg,nz); and 	strike = cell(nSeg,nz);
		% where each cell element contains "nPatch" values
		nz = size(handles.h_fault,2);
		handles.FaultDepth = [];
		for (k = 1:nSeg)
			nx = numel(handles.FaultStrike{k});
			handles.FaultDepth = [handles.FaultDepth; repmat({zeros(1, nx)}, nz, 1)];
		end
		
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
	end

	set(handles.edit_FaultWidth,'String',handles.FaultWidth{1}(1));
	set(handles.edit_FaultStrike,'String',num2str(handles.FaultStrike{1}(1)));
	set(handles.edit_FaultDip,'String',num2str(handles.FaultDip{1}(1)));
	set(handles.edit_FaultTopDepth,'String',num2str(handles.FaultTopDepth{1}(1)));
	set(handles.edit_FaultDepth,'String',num2str(handles.FaultDepth{1}(1)));
	set(handles.edit_DislocSlip,'String',num2str(handles.DislocSlip{1}(1)));
	set(handles.edit_DislocStrike,'String',num2str(handles.FaultStrike{1}(1)));
	set(handles.edit_DislocRake,'String',num2str(handles.DislocRake{1}(1)));
	
% ------------------------------------------------------------------------------------
function [handles, mag, M0] = compMag(handles, fault)
% Compute Moment magnitude
	mu = handles.mu * 1e10;
	M0 = mu * handles.um_milhao * handles.DislocSlip{fault}(:) .* handles.FaultWidth{fault}(:) .* ...
		handles.FaultLength{fault}(:);
	if (numel(M0) > 1),    M0 = sum(M0);   end
	mag = 2/3*(log10(M0) - 9.1);
	if (~isnan(mag))
		txt = ['Mw Magnitude = ' sprintf('%.1f',mag)];
		set(handles.h_txt_Mw,'String',txt,'Position',handles.txt_Mw_pos + [0 0 30 0])
		handles.Mw(fault) = mag;
	end

% ------------------------------------------------------------------------------------
function figure1_CloseRequestFcn(hObject, evt)
% ...
	handles = guidata(hObject);
	if (get(handles.check_SCC,'Val'))
		deleteFatias([],[], handles)
	elseif (~handles.fault_in)
		for (i=1:numel(handles.h_fault))
			hp = getappdata(handles.h_fault(i),'PatchHand');
			try     set(hp,'XData', [], 'YData',[], 'Visible','on');     end     % Had to use a try (f.. endless errors)
		end
	end
	delete(handles.figure1)

% -----------------------------------------------------------------------------------------
% ---------------- Creates and returns a handle to the GUI figure. 
function deform_mansinha_LayoutFcn(h1)

set(h1 ,...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'CloseRequestFcn',@figure1_CloseRequestFcn,...
'MenuBar','none',...
'Name','Vertical elastic deformation',...
'NumberTitle','off',...
'Position',[520 529 540 265],...
'Resize','off',...
'Tag','figure1',...
'HandleVisibility','callback');

uicontrol('Parent',h1,'Position',[10 126 181 131],'Style','frame');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@deform_mansinha_uiCB,...
'Position',[20 213 71 21],...
'Style','edit',...
'Tooltip','Fault length (km)',...
'Tag','edit_FaultLength');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@deform_mansinha_uiCB,...
'Position',[110 213 71 21],...
'Style','edit',...
'Tooltip','Fault width (km)',...
'Tag','edit_FaultWidth');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@deform_mansinha_uiCB,...
'Position',[20 173 71 21],...
'Style','edit',...
'Tooltip','Fault strike (degrees)',...
'Tag','edit_FaultStrike');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@deform_mansinha_uiCB,...
'Position',[110 173 71 21],...
'Style','edit',...
'Tooltip','Fault dip (degrees)',...
'Tag','edit_FaultDip');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@deform_mansinha_uiCB,...
'Position',[20 134 71 21],...
'Style','edit',...
'Tooltip','Depth of the base of fault''s plane',...
'Tag','edit_FaultDepth');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@deform_mansinha_uiCB,...
'Position',[110 133 71 21],...
'Style','edit',...
'Tooltip','Alternatively, give depth to the fault''s top ',...
'Tag','edit_FaultTopDepth');

uicontrol('Parent',h1,'Enable','inactive','Position',[36 235 41 13],'String','Length','Style','text');
uicontrol('Parent',h1,'Enable','inactive','Position',[125 236 41 13],'String','Width','Style','text');
uicontrol('Parent',h1,'Enable','inactive','Position',[34 195 41 13],'String','Strike','Style','text');
uicontrol('Parent',h1,'Enable','inactive','Position',[124 195 41 13],'String','Dip','Style','text');

uicontrol('Parent',h1,'Enable','inactive','Position',[108 154 75 16],...
'String','Depth to Top','Style','text',...
'Tooltip','Depth to the top of the fault (>= 0)');

uicontrol('Parent',h1,'Enable','inactive','Position',[34 155 41 16],...
'String','Depth','Style','text',...
'Tooltip','Depth to the top of the fault (>= 0)');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@deform_mansinha_uiCB,...
'Position',[210 216 91 22],...
'Style','popupmenu',...
'Tooltip','Set parameters with respect to this segment',...
'Value',1,...
'Tag','popup_segment');

uicontrol('Parent',h1,'Enable','inactive','Position',[229 238 48 16],...
'String','Segments','Style','text','Tag','txtFaultSeg');

uicontrol('Parent',h1,'Enable','inactive','Position',[53 250 85 15],...
'String','Fault Geometry','Style','text','Tag','faultGeom');

uicontrol('Parent',h1,'Position',[320 126 211 131],'Style','frame');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@deform_mansinha_uiCB,...
'Position',[330 213 51 21],...
'Style','edit',...
'Tag','edit_DislocStrike');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@deform_mansinha_uiCB,...
'Position',[400 213 51 21],...
'Style','edit',...
'Tooltip','Displacement angle anti-clockwise from horizontal',...
'Tag','edit_DislocRake');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@deform_mansinha_uiCB,...
'Position',[470 213 51 21],...
'Style','edit',...
'Tooltip','Total displacement',...
'Tag','edit_DislocSlip');

uicontrol('Parent',h1,'Enable','inactive','Position',[335 235 41 13],'String','Strike','Style','text');
uicontrol('Parent',h1,'Enable','inactive','Position',[404 235 41 13],'String','Rake','Style','text');
uicontrol('Parent',h1,'Enable','inactive','Position',[474 235 41 13],'String','Slip','Style','text');

uicontrol('Parent',h1,'Enable','inactive','Position',[373 250 111 15],...
'String','Dislocation Geometry','Style','text','Tag','dislocGeom');

uicontrol('Parent',h1,'Position',[10 11 350 93],'Style','frame');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@deform_mansinha_uiCB,...
'HorizontalAlignment','left',...
'Position',[76 64 71 21],...
'Style','edit',...
'Tooltip','X min value',...
'Tag','edit_x_min');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@deform_mansinha_uiCB,...
'HorizontalAlignment','left',...
'Position',[152 64 71 21],...
'Style','edit',...
'Tooltip','X max value',...
'Tag','edit_x_max');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@deform_mansinha_uiCB,...
'Position',[228 64 71 21],...
'Style','edit',...
'Tooltip','DX grid spacing',...
'Tag','edit_x_inc');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',{@deform_mansinha_uiCB,h1,'edit_Ncols_CB'},...
'Position',[304 64 45 21],...
'Style','edit',...
'Tooltip','Number of columns in the grid',...
'Tag','edit_Ncols');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@deform_mansinha_uiCB,...
'HorizontalAlignment','left',...
'Position',[76 38 71 21],...
'Style','edit',...
'Tooltip','Y min value',...
'Tag','edit_y_min');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@deform_mansinha_uiCB,...
'HorizontalAlignment','left',...
'Position',[152 38 71 21],...
'Style','edit',...
'Tooltip','Y max value',...
'Tag','edit_y_max');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@deform_mansinha_uiCB,...
'Position',[228 38 71 21],...
'Style','edit',...
'Tooltip','DY grid spacing',...
'Tag','edit_y_inc');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@deform_mansinha_uiCB,...
'Position',[304 38 45 21],...
'Style','edit',...
'Tooltip','Number of columns in the grid',...
'Tag','edit_Nrows');

uicontrol('Parent',h1,...
'BackgroundColor',[0.83137256 0.81568628 0.78431374],...
'Call',@deform_mansinha_uiCB,...
'FontWeight','bold',...
'ForegroundColor',[0 0 1],...
'Position',[289 16 61 18],...
'String','?',...
'Tag','push_Help_R');

uicontrol('Parent',h1,'Enable','inactive','Position',[18 69 55 15],'String','X Direction','Style','text');
uicontrol('Parent',h1,'Enable','inactive','Position',[17 43 55 15],'String','Y Direction','Style','text');
uicontrol('Parent',h1,'Enable','inactive','Position',[169 86 41 13],'String','Max','Style','text');
uicontrol('Parent',h1,'Enable','inactive','Position',[91 87 41 13],'String','Min','Style','text');
uicontrol('Parent',h1,'Enable','inactive','Position',[246 87 41 13],'String','Spacing','Style','text');
uicontrol('Parent',h1,'Enable','inactive','Position',[302 87 51 13],'String','# of lines','Style','text');

uicontrol('Parent',h1,'Enable','inactive','Position',[30 97 121 15],...
'String','Griding Line Geometry','Style','text','Tag','gridGeom');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@deform_mansinha_uiCB,...
'Position',[209 171 91 22],'Style','popupmenu',...
'Tooltip','Toggle between the different faults',...
'Value',1,'Tag','popup_fault');

uicontrol('Parent',h1,'Enable','inactive','Position',[229 194 42 15],...
'String','Faults','Style','text','Tag','txtFaultNum');

uicontrol('Parent',h1,...
'Call',@deform_mansinha_uiCB,...
'Position',[330 170 105 15],'String','Hide fault planes',...
'Style','checkbox','Tag','check_hideFaultPlanes');

uicontrol('Parent',h1,'Enable','inactive','FontSize',9,...
'HorizontalAlignment','left','Position',[327 135 105 16],...
'String','Mw Magnitude =','Style','text','Tag','h_txt_Mw');

uicontrol('Parent',h1,...
'Call',@deform_mansinha_uiCB,...
'Position',[450 170 50 15],'String','SCC',...
'Tooltip','Use variable slip',...
'Style','checkbox','Tag','check_SCC');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@deform_mansinha_uiCB,...
'Position',[450 135 31 21],...
'Style','edit','String','20',...
'Tooltip','Number of discretization slices',...
'Tag','edit_nSlices');

uicontrol('Parent',h1,'Enable','inactive','Position',[460 156 10 12],...
'String','N','Style','text','Tag','txtNfatias');

uicontrol('Parent',h1,'Enable','inactive','Position',[495 157 10 13],...
'String','q','Style','text','Tag','txtqValue');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@deform_mansinha_uiCB,...
'Position',[485 135 40 21],...
'Style','edit','String','0.3',...
'Tooltip','Slip skewness control (0 < q < 1)',...
'Tag','edit_qValue');

uicontrol('Parent',h1,'Position',[224 150 60 15],'ForegroundColor',[1 0 0], 'String','CONFIRM','Style','text');

uicontrol('Parent',h1, 'Position',[410 100 60 18],...
'String','Mu (x10^10)',...
'Style','text',...
'Tag','text_mu');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@deform_mansinha_uiCB,...
'Position',[485 100 40 21],...
'String','3.0',...
'Style','edit',...
'Tooltip','Shear modulus (for Mw calculation)',...
'Tag','edit_mu');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@deform_mansinha_uiCB,...
'String', {'Geogs' 'Meters' 'Kilometers'},...
'Position',[209 127 91 22],'Style','popupmenu',...
'Tooltip','GRID COORDINATES: IT IS YOUR RESPONSABILITY THAT THIS IS CORRECT',...
'Value',1,'Tag','popup_GridCoords');

uicontrol('Parent',h1,...
'Call',@deform_mansinha_uiCB,...
'Position',[380 51 29 29],...
'Tooltip','Show focal mechanism',...
'Tag','push_focal');

uicontrol('Parent',h1,...
'Call',@deform_mansinha_uiCB,...
'FontWeight','bold',...
'Position',[370 15 71 21],...
'Tooltip','Save current solution in sub-fault format',...
'String','Save fault',...
'Tag','push_save_subfault');

uicontrol('Parent',h1,...
'Call',@deform_mansinha_uiCB,...
'FontWeight','bold',...
'Position',[460 15 71 21],...
'String','Compute',...
'Tag','push_compute');

function deform_mansinha_uiCB(hObject, eventdata)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));

function varargout = microlev(varargin)
% Helper window to perform microlleveling filter on magnetic grids

%	Copyright (c) 2004-2015 by J. Luis
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

	if (isempty(varargin) || ~ishandle( varargin{1}))
		errordlg('MICROLEV: input argument must be a Mirone figure handle.','Error'),	return
	end

	if ( strcmp(get(varargin{1},'type'), 'line') )	% The handle is actualy a line (rectangle) handle
		hMirFig = get(get(varargin{1},'Parent'),'Parent');
		handMir = guidata(hMirFig);
		gotROI = true;
	else
		hMirFig = varargin{1};
		handMir = guidata(hMirFig);
		gotROI = false;
	end

	if (handMir.no_file)
		errordlg('You didn''t even load a file. What are you expecting then?','ERROR'),	return
	end
	if (~handMir.validGrid)
		errordlg('This operation can only be performed on grids.','Error'),	return
	end

	hObject = figure('Tag','figure1','Visible','off');
	microlev_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hMirFig,hObject)

	handles.hMirFig = hMirFig;
	[handles.X, handles.Y, handles.Z, head] = load_grd(handMir);
	handles.head = head;
	handles.have_nans = handMir.have_nans;
	handles.home_dir = handMir.home_dir;
	handles.work_dir = handMir.work_dir;
	handles.last_dir = handMir.last_dir;

	handles.Z_filt = [];
	handles.cmboMask = [];
	handles.resMask = [];
	handles.trackMask = [];
	handles.micro = [];
	handles.nWin = 3;
	handles.rmsThresh = 6;
	handles.lineThick = 2;
	handles.command{1} = '-Fc';
	handles.gotROI = gotROI;
	if (handMir.geog)
		set(handles.popup_Option_D, 'Val', 2)
		handles.command{3} = '-D1';
	else
		handles.command{3} = '-D0';
	end

	if (gotROI)
		pos = get(handles.push_computeSolution, 'Pos');
		set(handles.push_computeSolution, 'Pos', [pos(1) pos(2)+70 pos(3:4)])	% Move this button a bit up
		set(handles.push_burnSolution, 'Vis', 'on')
		handles.hRect = varargin{1};
		x = get(varargin{1},'XData');		y = get(varargin{1},'YData');
		rect = [x(1) y(1) (x(3)-x(2)) (y(2)-y(1))];
		[handles.Z_rect, r_c] = cropimg(head(1:2),head(3:4),handles.Z, rect, 'out_grid');
		handles.head_orig = head;
		
		X = linspace( head(1) + (r_c(3)-1)*head(8), head(1) + (r_c(4)-1)*head(8), r_c(4) - r_c(3) + 1 );
		Y = linspace( head(3) + (r_c(1)-1)*head(9), head(3) + (r_c(2)-1)*head(9), r_c(2) - r_c(1) + 1 );
		head(1) = X(1);		head(2) = X(end);		head(3) = Y(1);		head(4) = Y(end);
		
		handles.head = head;		% min/max is wrong
		handles.X = X;
		handles.Y = Y;
		handles.row_col = r_c;
	end

	%------------ Give a Pro look (3D) to the frame boxes  -------------------------------
	new_frame3D(hObject, [handles.txt_step1 handles.txt_step2 handles.txt_step3 handles.txt_step4])
	%------------- END Pro look (3D) -----------------------------------------------------

	% ------------ Create a Semaforo for track lines ops -------------
	[semaforo, pal] = aux_funs('semaforo_green');
	handles.hSemaforo = image(semaforo,'Parent', handles.axes1);
	set(handles.axes1, 'XTick',[], 'YTick', [])
	set(handles.figure1, 'Colormap', pal),		drawnow

	set(hObject,'Visible','on');

	guidata(hObject, handles);
	if (nargout),	varargout{1} = 	hObject;	end

% -------------------------------------------------------------------------------------
function popup_FilterType_CB(hObject, handles)
	val = get(hObject,'Value');		str = get(hObject, 'String');
	switch str{val};
		case 'boxcar',				handles.command{1} = '-Fb';
		case 'cosine arch',			handles.command{1} = '-Fc';
		case 'gaussian',			handles.command{1} = '-Fg';
	end
	guidata(hObject, handles);

% -------------------------------------------------------------------------------------
function edit_FilterWidth_CB(hObject, handles)
	xx = get(hObject,'String');
	if ~isempty(xx)		handles.command{2} = xx;
	else				handles.command{2} = [];
	end
	guidata(hObject,handles)

% -------------------------------------------------------------------------------------
function push_HelpFilterType_CB(hObject, handles)
	helpdlg('Choose one only of for boxcar, cosine arch or gaussian filter and specify full width.','Help');

% -------------------------------------------------------------------------------------
function popup_Option_D_CB(hObject, handles)
	val = get(hObject,'Value');     str = get(hObject, 'String');
	switch str{val};
		case '0',        handles.command{3} = '-D0';
		case '1',        handles.command{3} = '-D1';
		case '2',        handles.command{3} = '-D2';
		case '3',        handles.command{3} = '-D3';
		case '4',        handles.command{3} = '-D4'; 
	end
	guidata(handles.figure1, handles);

% -------------------------------------------------------------------------------------
function push_HelpOptionD_CB(hObject, handles)
message = {'Distance flag tells how grid (x,y) relates to filter width as follows:'
			' '
			'flag = 0: grid (x,y) same units as width, Cartesian distances.'
			'flag = 1: grid (x,y) in degrees, width  in  kilometers, Cartesian distances.'
			'flag  =  2:  grid (x,y) in degrees, width in km, dx scaled by cos(middle y), Cartesian distances.'
			' '
			'The above options are fastest  because  they  allow'
			'weight  matrix  to be computed only once.  The next'
			'two  options  are  slower  because  they  recompute'
			'weights for each East-West scan line.'
			' '
			'flag  =  3:  grid (x,y) in degrees, width in km, dx scaled by cosine(y), Cartesian distance calculation.'
			'flag  =  4:  grid  (x,y)  in  degrees, width in km, Spherical distance calculation.'};
helpdlg(message,'Help Distance flag');

% -------------------------------------------------------------------------------------
function push_applyStep1_CB(hObject, handles)
% STEP 1. Compute filtered grid

	if isempty(handles.command{2})
		errordlg('Must specify a Filter width','Error');    return
	end

	opt_F = [handles.command{1} handles.command{2}];
	opt_D = handles.command{3};

	set(handles.figure1,'pointer','watch')
	if (~handles.gotROI)
		[handles.Z_filt, head] = c_grdfilter(handles.Z,handles.head,opt_F,opt_D);
	else
		[handles.Z_filt, head] = c_grdfilter(handles.Z_rect,handles.head,opt_F,opt_D);
	end
	update_report(handles, 1)		% Update the listbox info
	set(handles.figure1,'pointer','arrow')
	guidata(handles.figure1, handles);

	if (get(handles.check_show1,'Val'))
		zMinMax = grdutils(handles.Z_filt,'-L');
		tmp.head = [head(1:4) zMinMax(1:2)' head(7:9)];
		tmp.X = handles.X;			tmp.Y = handles.Y;		tmp.name = 'Filtered grid';
		mirone(handles.Z_filt, tmp)
	end

% -------------------------------------------------------------------------------------
function push_loadFiltered_CB(hObject, handles)
% STEP 1. OR Load a grid previously filtered grid
	[FileName,PathName] = put_or_get_file(handles, ...
		{'*.grd;*.nc', 'Grid files (*.grd,*.nc)';'*.*', 'All Files (*.*)'},'Select grid','get');
	if (isequal(FileName,0))		return,		end
	fname = [PathName FileName];
	set(handles.edit_loadFiltered, 'String', fname)

	% Let the edit_loadFiltered_CB do the rest of the work
	edit_loadFiltered_CB(handles.edit_loadFiltered, handles)

% -------------------------------------------------------------------------------------
function edit_loadFiltered_CB(hObject, handles)
% STEP 1. 
	fname = get(hObject, 'String');
	if isempty(fname),		return,		end

	[handles, X, Y, handles.Z_filt, head_filt] = read_gmt_type_grids(handles,fname);
	if (isempty(X)),    return,		end
	
	if ( ~isequal(size(handles.Z), size(handles.Z_filt)) )
		errordlg('Error: main and filtered grid have different sizes','Error'),	return
	end
	if ( any( abs(handles.head(1:4) - head_filt(1:4)) > 1e-5 ) )
		errordlg('Error: main and filtered grid do not cover the same region','Error'),	return
	end

	update_report(handles, 1)		% Update the listbox info
	guidata(handles.figure1,handles)

	if (get(handles.check_show1,'Val'))
		tmp.head = head_filt;		tmp.name = 'Filtered grid';
		tmp.X = handles.X;			tmp.Y = handles.Y;
		mirone(handles.Z_filt, tmp)
	end

% -------------------------------------------------------------------------------------
function push_applyStep2_CB(hObject, handles)
% STEP 2. Compute the highpass or residue grid. That is, original - filtered

	if (~handles.gotROI)
		handles.Z_resid = cvlib_mex('sub', handles.Z, handles.Z_filt);
	else
		handles.Z_resid = cvlib_mex('sub', handles.Z_rect, handles.Z_filt);
	end
	update_report(handles, 2)		% Update the listbox info

	if (get(handles.check_show2,'Val'))
		tmp.X = handles.X;			tmp.Y = handles.Y;
		tmp.head = handles.head;	tmp.name = 'Highpass grid';		% Wrong min/max
		mirone(handles.Z_resid, tmp)
	end
	guidata(handles.figure1, handles);

% -------------------------------------------------------------------------------------
function edit_nWin_CB(hObject, handles)
% STEP 3. Set the window (neighborwood) size
	xx = round( str2double(get(hObject,'String')) );
	if (isnan(xx)),		set(hObject,'String',handles.nWin),	return,		end
	if (rem(xx,2) == 0)
		errordlg('Window dimension must be an odd number. Fixing it','Error')
		xx = xx + 1;
	end
	handles.nWin = xx;
	if (handles.nWin < 3)       % Minimum allowed is 3
		set(hObject,'String',3)
		handles.nWin = 3;
	end
	guidata(handles.figure1,handles)

% -------------------------------------------------------------------------------------
function edit_rmsThreshold_CB(hObject, handles)
% STEP 3. Set the RMS threshold
	xx = str2double(get(hObject,'String'));
	if (isnan(xx) || xx < 0),		set(hObject,'String',handles.rmsThresh),	return,		end
	handles.rmsThresh = xx;
	guidata(handles.figure1,handles)

% -------------------------------------------------------------------------------------
function push_applyStep3_CB(hObject, handles)
% STEP 3. Compute a mask from clipping a windowed RMS solution

	set(handles.figure1,'pointer','watch')
	opt_W = sprintf('-W%d', handles.nWin);
	opt_N = sprintf('-N%d', handles.have_nans);
	rms = mirblock(handles.Z_resid, handles.head, '-A8', opt_W, opt_N);
	handles.resMask = (rms >= handles.rmsThresh);	clear rms
	set(handles.figure1,'pointer','arrow')

	update_report(handles, 3)		% Update the listbox info
	guidata(handles.figure1,handles)

	if (get(handles.check_show3,'Val'))
		tmp.X = [handles.X(1) handles.X(end)];		tmp.Y = [handles.Y(1) handles.Y(end)];
		tmp.head = handles.head;					tmp.head(5:6) = [0 1];		tmp.name = 'Clip mask';
		mirone(handles.resMask, tmp)
	end

% -------------------------------------------------------------------------------------
function push_loadMask1_CB(hObject, handles)
% STEP 3. OR Load a grid previously computed mask grid
	[FileName,PathName] = put_or_get_file(handles, ...
		{'*.grd;*.nc', 'Grid files (*.grd,*.nc)';'*.*', 'All Files (*.*)'},'Select grid','get');
	if (isequal(FileName,0))		return,		end
	fname = [PathName FileName];
	set(handles.edit_loadMask1, 'String', fname)

	% Let the edit_loadMask1_CB do the rest of the work
	edit_loadMask1_CB(handles.edit_loadMask1, handles)

% -------------------------------------------------------------------------------------
function edit_loadMask1_CB(hObject, handles)
% STEP 3. 
	fname = get(hObject, 'String');
	if isempty(fname),		return,		end

	have_nans = handles.have_nans;		% backup this, which will be reset in next call
	[handles, X, Y, handles.resMask, head_mask] = read_gmt_type_grids(handles,fname);
	if (isempty(X)),    return,		end
	
	if ( ~isequal(size(handles.Z), size(handles.resMask)) )
		errordlg('Error: main and Mask grid have different sizes','Error'),	return
	end
	if ( any( abs(handles.head(1:4) - head_mask(1:4)) > 1e-5 ) )
		errordlg('Error: main and Mask grid do not cover the same region','Error'),	return
	end
	if (head_mask(5) ~= 0 || head_mask(6) ~= 1)
		warndlg('This grid isn''t really a Mask since it has values different from either 0 or 1.','Warning')
	end

	if (handles.have_nans)			% We cannot try to convert NaNs to logical
		handles.resMask(isnan(handles.resMask)) = 0;
	end
	handles.resMask = logical(handles.resMask);
	handles.have_nans = have_nans;	% Reset original value

	update_report(handles, 3)		% Update the listbox info
	guidata(handles.figure1,handles)

	if (get(handles.check_show3,'Val'))
		tmp.X = [handles.X(1) handles.X(end)];		tmp.Y = [handles.Y(1) handles.Y(end)];
		tmp.head = handles.head;					tmp.head(5:6) = [0 1];		tmp.name = 'Clip mask';
		mirone(handles.resMask, tmp)
	end

% -------------------------------------------------------------------------------------
function edit_lineThick_CB(hObject, handles)
% STEP 4 (optional).
	xx = str2double(get(hObject,'String'));
	if (isnan(xx) || xx < 1),		set(hObject,'String',handles.lineThick),	return,		end
	handles.lineThick = xx;
	guidata(handles.figure1,handles)

% -------------------------------------------------------------------------------------
function push_loadXY_CB(hObject, handles)
% STEP 4 (optional).
	[FileName,PathName] = put_or_get_file(handles, ...
		{'*.dat;*.txt', 'Grid files (*.dat,*.txt)';'*.*', 'All Files (*.*)'},'Select grid','get');
	if (isequal(FileName,0))		return,		end
	fname = [PathName FileName];
	set(handles.edit_loadXY, 'String', fname)

% -------------------------------------------------------------------------------------
function push_maskFromLines_CB(hObject, handles)
% STEP 4 (optional). Make a mask from track lines. This is done via a screen capture
% so the linewidth matter a lot

	fnameXY = get(handles.edit_loadXY, 'String');
	if ( isempty(fnameXY) && ~get(handles.check_useMGG, 'Val') )
		errordlg('You must choose someting to do. Look better to this section optios.','Error')
		return
	end
	
	[semaforo, pal] = aux_funs('semaforo_red');
	set(handles.hSemaforo,'CData',semaforo)
	set(handles.figure1, 'Colormap', pal),		drawnow

	% Either case we have a use to these, so create them right away
	tmp.X = [handles.X(1) handles.X(end)];		tmp.Y = [handles.Y(1) handles.Y(end)];
	tmp.head(5:6) = [0 1];						tmp.head = handles.head;	tmp.name = 'Nikles';
	if (~handles.gotROI)
		mask = true(size(handles.Z));	% Create a fake base image where to plot the tracks
	else
		mask = true(size(handles.Z_rect));
	end

	hMirFig = mirone(mask,tmp);		set(hMirFig,'Vis', 'off'),		pause(0.01)
	%set(findobj(hMirFig,'type','image'), 'Vis', 'off')

	if (get(handles.check_useMGG, 'Val'))
		x = [tmp.head(1) tmp.head(1) tmp.head(2)];
		y = [tmp.head(3) tmp.head(4) tmp.head(4)];
		tracks = deal_opts([], 'get_MGGtracks', [],[],x, y);	% 2 [] because get_MGGtracks is a @fun
		mirone('GeophysicsImportGmtFile_CB',guidata(hMirFig), tracks)
	else
		try			load_xyz(guidata(hMirFig), fnameXY)
		catch,		errordlg(lasterr,'Error'),	return
		end
	end

	hL = findobj(hMirFig,'Type','line');	% Fish the line handles
	set(hL,'LineWidth',handles.lineThick)
	mask = imcapture(findobj(hMirFig,'type','axes'), 'img', 0);
	delete(hMirFig)							% Delete the hiden Mirone figure
	mask = mask(:,:,1);
	handles.trackMask = flipud(mask ~= 255);

	[semaforo, pal] = aux_funs('semaforo_green');
	set(handles.hSemaforo,'CData',semaforo)
	set(handles.figure1, 'Colormap', pal),		drawnow

	update_report(handles, 5)		% Update the listbox info (5 because 4rth entry is '-----')
	guidata(handles.figure1,handles)

	if (get(handles.check_show4a,'Val'))
		tmp.name = 'Tracks mask';
		mirone(handles.trackMask, tmp)
	end

% -------------------------------------------------------------------------------------
function push_combineMasks_CB(hObject, handles)
% STEP 4 (optional). Combine (logical &) clipp and tracks masks
	status = update_report(handles, 3);
	if (~status)
		errordlg('Mask from STEP 3 is unknown. You have to compute it first.','Error')
		return
	end

	status = update_report(handles, 5);
	if (~status)
		errordlg('Mask from STEP 4a is unknown. To use this option you have to compute it first.','Error')
		return
	end

	% Now combine the Masks
	handles.comboMask = handles.resMask & handles.trackMask;
	update_report(handles, 6)
	guidata(handles.figure1,handles)

	if (get(handles.check_show4b,'Val'))
		tmp.X = [handles.X(1) handles.X(end)];		tmp.Y = [handles.Y(1) handles.Y(end)];
		tmp.head = handles.head;					tmp.head(5:6) = [0 1];		tmp.name = 'Combined mask';
		mirone(handles.comboMask, tmp)
	end

% -------------------------------------------------------------------------------------
function push_computeSolution_CB(hObject, handles)
% STEP 5. Final

	msg = cell(4,1);	one = false;	c = true(1,4);
	for (k = 1:2)		% After a user request, STEP3 is no longer mandatory. But 3 or 4 must exist
		if (~update_report(handles, k))
			msg{k} = sprintf('STEP %d is not done', k);
			one = true;		c(k) = false;
		end
	end
	if (one)		% Have at least one mandatory step not done
		msg{1} = 'ERROR';
		c(1) = false;		msg(c) = [];	% Remove empty cells that bother poor errordlg
		errordlg(msg, 'ERROR'),		return
	end

	if ( update_report(handles, 5) && ~update_report(handles, 6))
		warndlg('You did STEP 4a but not STEP 4b (mask combining)');
	end

	if ( isempty(handles.comboMask) && isempty(handles.resMask) && isempty(handles.trackMask) )
		errordlg('At least one of the STEP 3 or STEP 4 must be computed.', 'ERROR'),	return
	end

	% Find out wich mask to use
	if (~isempty(handles.comboMask)),		mask = handles.comboMask;
	elseif (~isempty(handles.resMask))		mask = handles.resMask;
	else									mask = handles.trackMask;
	end
	mask = ~mask;		% Negate because we want to add back the non-high order variations

	micro = handles.Z_filt;
	micro(mask) = cvlib_mex('add', micro(mask), handles.Z_resid(mask));

	zMinMax = grdutils(handles.Z_filt,'-L');
	tmp.head = [handles.head(1:4) zMinMax(1:2)' handles.head(7:9)];
	tmp.X = handles.X;			tmp.Y = handles.Y;		tmp.name = 'MicroLev grid';
	mirone(micro, tmp)

	if (handles.gotROI)				% Need to save the solution
		handles.micro = micro;
		guidata(handles.figure1, handles)
	end

% -------------------------------------------------------------------------------------
function push_burnSolution_CB(hObject, handles)
% Apply the STEP 5 into the big grid
	if (isempty(handles.micro))
		errordlg('You did not finish STEP 5','Error'),		return
	end
	handles.Z(handles.row_col(1):handles.row_col(2), handles.row_col(3):handles.row_col(4)) = handles.micro;
	setappdata(handles.hMirFig,'dem_z', handles.Z);			% Update the main grid (still need Min/Max)

% -------------------------------------------------------------------------------------
function resp = update_report(handles, n)
% Set or get status of N step.
% Set status if called with naragout == 0. Otherwise, get status.
	str = get(handles.listbox, 'String');
	ind = strfind(str{n},'(');
	if (nargout)
		if (strcmp(str{n}(ind+1:end-1), 'not done'))
			resp = false;
		else
			resp = true;
		end
	else
		str{n} = [str{n}(1:ind) 'done)'];
		set(handles.listbox, 'String', str, 'Val', n)		
	end

% -------------------------------------------------------------------------------------
function microlev_LayoutFcn(h1)

set(h1, 'Position',[520 239 470 561],...
'Color', get(0,'factoryUicontrolBackgroundColor'),...
'MenuBar','none',...
'Name','Microleveling',...
'NumberTitle','off',...
'Resize','off',...
'HandleVisibility','Call');

uicontrol('Parent',h1,'Position',[10 440 451 111], 'Style','frame');
uicontrol('Parent',h1,'Position',[10 380 451 51],  'Style','frame');
uicontrol('Parent',h1,'Position',[10 270 451 101], 'Style','frame');
uicontrol('Parent',h1,'Position',[10 130 452 131], 'Style','frame');

uicontrol('Parent',h1,'Position',[22 501 122 22],...
'BackgroundColor',[1 1 1],...
'Call',@microlev_uiCB,...
'HorizontalAlignment','right',...
'String',{'boxcar'; 'cosine arch'; 'gaussian' },...
'Style','popupmenu',...
'Value',2,...
'Tag','popup_FilterType');

uicontrol('Parent',h1,'Position',[145 502 47 21],...
'BackgroundColor',[1 0.49 0.49],...
'Call',@microlev_uiCB,...
'HorizontalAlignment','center',...
'String','',...
'Style','edit',...
'TooltipString','specify filter full width',...
'Tag','edit_FilterWidth');

uicontrol('Parent',h1,'Position',[193 501 21 23],...
'Call',@microlev_uiCB,...
'FontWeight','bold',...
'ForegroundColor',[0 0 1],...
'String','?',...
'Tag','push_HelpFilterType');

uicontrol('Parent',h1,'Position',[231 503 72 22],...
'BackgroundColor',[1 1 1],...
'Call',@microlev_uiCB,...
'String',{'0'; '1'; '2'; '3'; '4' },...
'Style','popupmenu',...
'TooltipString','You better read the help availabe on box aside',...
'Value',1,...
'Tag','popup_Option_D');

uicontrol('Parent',h1,'Position',[305 503 21 23],...
'Call',@microlev_uiCB,...
'FontWeight','bold',...
'ForegroundColor',[0 0 1],...
'String','?',...
'Tag','push_HelpOptionD');

uicontrol('Parent',h1,'Position',[55 525 47 15],...
'HorizontalAlignment','left',...
'String','Filter type',...
'Style','text');

uicontrol('Parent',h1,'Position',[145 524 53 15],...
'HorizontalAlignment','left',...
'String','Filter width',...
'Style','text');

uicontrol('Parent',h1,'Position',[234 526 63 15],...
'HorizontalAlignment','left',...
'String','Distance flag',...
'Style','text');

uicontrol('Parent',h1,'Position',[20 451 311 21],...
'BackgroundColor',[1 1 1],...
'Call',@microlev_uiCB,...
'HorizontalAlignment','left',...
'String','',...
'Style','edit',...
'TooltipString','Load a previously filtered grid',...
'Tag','edit_loadFiltered');

uicontrol('Parent',h1,'Position',[330 450 21 23],...
'Call',@microlev_uiCB,...
'FontSize',9,...
'FontWeight','bold',...
'String','...',...
'Tag','push_loadFiltered');

uicontrol('Parent',h1,'Position',[102 540 150 21],...
'FontSize',9,...
'FontWeight','bold',...
'HorizontalAlignment','left',...
'String',' STEP 1 -- Filtering',...
'Style','text',...
'Tag','txt_step1');

uicontrol('Parent',h1,'Position',[80 417 220 21],...
'FontSize',9,...
'FontWeight','bold',...
'HorizontalAlignment','left',...
'String',' STEP 2 -- High Pass (residues)',...
'Style','text',...
'Tag','txt_step2');

uicontrol('Parent',h1,'Position',[80 357 180 21],...
'FontSize',9,...
'FontWeight','bold',...
'HorizontalAlignment','left',...
'String',' STEP 3 -- Clip residues',...
'Style','text',...
'Tag','txt_step3');

uicontrol('Parent',h1,'Position',[80 246 280 22],...
'FontSize',9,...
'FontWeight','bold',...
'HorizontalAlignment','left',...
'String',' STEP 4 -- Create second mask (optional)',...
'Style','text',...
'Tag','txt_step4');

uicontrol('Parent',h1,'Position',[381 506 69 21],...
'Call',@microlev_uiCB,...
'String','Apply',...
'Tag','push_applyStep1');

uicontrol('Parent',h1,'Position',[370 468 89 23],...
'String','Show result',...
'Style','checkbox',...
'TooltipString','Create a new window showing the result of this step',...
'Tag','check_show1');

uicontrol('Parent',h1,'Position',[161 391 90 21],...
'Call',@microlev_uiCB,...
'String','Apply',...
'TooltipString','Compute difference between original and filtered grid',...
'Tag','push_applyStep2');

uicontrol('Parent',h1,'Position',[370 391 89 23],...
'String','Show result',...
'Style','checkbox',...
'TooltipString','Create a new window showing the result of this step',...
'Tag','check_show2');

uicontrol('Parent',h1,'Position',[109 330 31 21],...
'BackgroundColor',[1 1 1],...
'Call',@microlev_uiCB,...
'String','3',...
'Style','edit',...
'TooltipString','Width of the rectangular neighborhood (MUST bo an odd number)',...
'Tag','edit_nWin');

uicontrol('Parent',h1,'Position',[30 332 78 17],...
'FontName','Helvetica',...
'FontSize',9,...
'HorizontalAlignment','right',...
'String','Window size',...
'Style','text',...
'TooltipString','Width of the rectangular neighborhood (MUST bo an odd number)');

uicontrol('Parent',h1,'Position',[259 330 31 21],...
'BackgroundColor',[1 1 1],...
'Call',@microlev_uiCB,...
'String','6',...
'Style','edit',...
'TooltipString','RMS values greater or higher than this are set to 1, and others to 0',...
'Tag','edit_rmsThreshold');

uicontrol('Parent',h1,'Position',[160 330 98 19],...
'FontName','Helvetica',...
'FontSize',9,...
'HorizontalAlignment','right',...
'String','RMS threshold',...
'Style','text');

uicontrol('Parent',h1,'Position',[20 280 311 21],...
'BackgroundColor',[1 1 1],...
'Call',@microlev_uiCB,...
'HorizontalAlignment','left',...
'String','',...
'Style','edit',...
'TooltipString','Load a previously created mask with 1s at positions where "High Pass" will be deleted',...
'Tag','edit_loadMask1');

uicontrol('Parent',h1,'Position',[330 279 21 23],...
'Call',@microlev_uiCB,...
'FontSize',9,...
'FontWeight','bold',...
'String','...',...
'Tag','push_loadMask1');

uicontrol('Parent',h1,'Position',[381 330 69 21],...
'Call',@microlev_uiCB,...
'String','Apply',...
'TooltipString','Compute a mask grid with condition set by "RMS threshold"',...
'Tag','push_applyStep3');

uicontrol('Parent',h1,'Position',[370 296 89 23],...
'String','Show result',...
'Style','checkbox',...
'TooltipString','Create a new window showing the result of this step',...
'Tag','check_show3');

uicontrol('Parent',h1,'Position',[30 218 110 24],...
'String','Use MGG tracks',...
'Style','checkbox',...
'TooltipString','If you have mag data in mgd77+ netCDF files and are using GMT''s x2sys supplement',...
'Tag','check_useMGG');

uicontrol('Parent',h1,'Position',[223 219 31 21],...
'BackgroundColor',[1 1 1],...
'Call',@microlev_uiCB,...
'String','2',...
'Style','edit',...
'TooltipString','Line thickness in points used when drawing lines (not seen). Will impact on mask raster result.',...
'Tag','edit_lineThick');

uicontrol('Parent',h1,'Position',[298 220 147 21],...
'Call',@microlev_uiCB,...
'String','Make mask from track lines',...
'TooltipString','Create a mask grid with 1s at track lines position',...
'Tag','push_maskFromLines');

uicontrol('Parent',h1,'Position',[20 170 311 21],...
'BackgroundColor',[1 1 1],...
'HorizontalAlignment','left',...
'String','',...
'Style','edit',...
'TooltipString','Load a (x,y) multi-segment file with the track lines navigation',...
'Tag','edit_loadXY');

uicontrol('Parent',h1,'Position',[330 169 21 23],...
'Call',@microlev_uiCB,...
'FontSize',9,...
'FontWeight','bold',...
'String','...',...
'Tag','push_loadXY');

uicontrol('Parent',h1,'Position',[370 190 89 21],...
'String','Show result',...
'Style','checkbox',...
'TooltipString','Show result of mask from lines',...
'Tag','check_show4a');

uicontrol('Parent',h1,'Position',[370 138 89 23],...
'String','Show result',...
'Style','checkbox',...
'TooltipString','Show result of mask combination',...
'Tag','check_show4b');

uicontrol('Parent',h1,'Position',[90 136 221 21],...
'Call',@microlev_uiCB,...
'String','Combine Masks (MaskA & MaskB)',...
'Tag','push_combineMasks');

axes('Parent',h1,'Units','pixels','Position',[446 212 14 46],...
'XTick', [], 'YTick', [], ...
'Visible', 'off', ...
'Tag','axes1');

uicontrol('Parent',h1,'Position',[238 55 222 21],...
'FontSize',11,...
'FontWeight','bold',...
'HorizontalAlignment','left',...
'String','STEP 5 -- Compute solution',...
'Style','text');

uicontrol('Parent',h1,'Position',[237 20 221 21],...
'Call',@microlev_uiCB,...
'String','Apply Mask and compute solution',...
'Tag','push_computeSolution');

uicontrol('Parent',h1,'Position',[237 20 221 21],...
'Call',@microlev_uiCB,...
'String','Burn solution in Big GRID',...
'Vis', 'off',...
'Tag','push_burnSolution');

uicontrol('Parent',h1,'Position',[145 474 85 16],...
'FontSize',9,...
'String','OR, Load file',...
'Style','text');

uicontrol('Parent',h1,'Position',[149 301 120 17],...
'FontSize',9,...
'String','OR, Load a Mask file',...
'Style','text');

uicontrol('Parent',h1,'Position',[145 191 120 17],...
'FontSize',9,...
'String','OR, Load a (x,y) file',...
'Style','text');

uicontrol('Parent',h1,'Position',[39 10 171 101],...
'BackgroundColor',[1 1 1],...
'String',{'STEP 1 (not done)'; 'STEP 2 (not done)'; 'STEP 3 (not done)'; '----- Optional --------'; 'STEP 4a (not done)'; 'STEP 4b (not done)' },...
'Style','listbox',...
'Value',1,...
'Tag','listbox');

uicontrol('Parent',h1,'Position',[156 220 65 17],...
'FontName','Helvetica',...
'FontSize',9,...
'HorizontalAlignment','right',...
'String','Line width',...
'Style','text',...
'TooltipString','Line thickness in points');

function microlev_uiCB(hObject, evt)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));

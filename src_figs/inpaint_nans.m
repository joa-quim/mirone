function varargout = inpaint_nans(varargin)
% Fill holes in 2D arrays
 
	hObject = figure('Tag','figure1','Visible','off');
	inpaint_nans_LayoutFcn(hObject);
	handles = guihandles(hObject);
	movegui(hObject, 'center')
 
	if (numel(varargin) > 0)
		handles.handMir = varargin{1};
	else
		delete(hObject)
		return
	end

	if (~handles.handMir.have_nans)
		warndlg('Yeah, fill your head. No NaNs in this grid.','pseudo-Warning')
		delete(hObject),		return
	end

	handles.nCells = 10;			% Default max number of "Paint small" option

	%------------ Give a Pro look (3D) to the frame boxes  -------------------------------
	bgcolor = get(0,'DefaultUicontrolBackgroundColor');
	framecolor = max(min(0.65*bgcolor,[1 1 1]),[0 0 0]);
	frames = [handles.frame1 handles.frame2];
	for (k = 1:numel(frames))
    	frame_size = get(frames(k),'Position');
		frame3D(hObject,frame_size,framecolor,'',get(frames(k),'UserData'))
		delete(frames(k))
	end

	% Recopy the text fields on top of previously created frames
	h_t = [handles.text_what; handles.text_method];
	% This is a destilled minimalist uistack(H,'top')
	Children = findobj(allchild(hObject), 'flat', 'type', 'uicontrol');
	HandleLoc = find(ismember(Children,h_t));
	Children(HandleLoc) = [];  
	NewOrder = [h_t; Children];		% 'top'
	AllChildren = allchild(hObject);
	AllChildren(find(ismember(AllChildren,NewOrder))) = NewOrder;
	set(hObject,'Children',AllChildren);
	%------------- END Pro look (3D) -------------------------------------------------------

	% Add this figure handle to the carraças list
	plugedWin = getappdata(handles.handMir.figure1,'dependentFigs');
	plugedWin = [plugedWin hObject];
	setappdata(handles.handMir.figure1,'dependentFigs',plugedWin);

	guidata(hObject, handles);
	set(hObject,'Visible','on');
	if (nargout),   varargout{1} = hObject;     end

% -----------------------------------------------------------------------------------------
function radio_paintAll_Callback(hObject, eventdata, handles)
	if (~get(hObject,'Value')),		set(hObject,'Value',1),		return,		end
	set(handles.radio_paintSmall,'Val',0)
	set(handles.edit_nCells,'Enable','on')

% -----------------------------------------------------------------------------------------
function radio_paintSmall_Callback(hObject, eventdata, handles)
	if (~get(hObject,'Value')),		set(hObject,'Value',1),		return,		end
	set(handles.radio_paintAll,'Val',0)
	set(handles.edit_nCells,'Enable','off')

% -----------------------------------------------------------------------------------------
function edit_nCells_Callback(hObject, eventdata, handles)
	xx = abs( round(str2double(get(hObject,'String'))) );
	if (isnan(xx))
		set(hObject,'String',handles.nCells)
		return
	end
	handles.nCells = xx;
	guidata(handles.figure1, handles)

% -----------------------------------------------------------------------------------------
function radio_surface_Callback(hObject, eventdata, handles)
	if (~get(hObject,'Value')),		set(hObject,'Value',1),		return,		end
	set([handles.radio_bicubic handles.radio_bilinear],'Val',0)

% -----------------------------------------------------------------------------------------
function radio_bicubic_Callback(hObject, eventdata, handles)
	if (~get(hObject,'Value')),		set(hObject,'Value',1),		return,		end
	set([handles.radio_surface handles.radio_bilinear],'Val',0)

% -----------------------------------------------------------------------------------------
function radio_bilinear_Callback(hObject, eventdata, handles)
	if (~get(hObject,'Value')),		set(hObject,'Value',1),		return,		end
	set([handles.radio_surface handles.radio_bicubic],'Val',0)

% -----------------------------------------------------------------------------------------
function radio_inPlace_Callback(hObject, eventdata, handles)
	if (~get(hObject,'Value')),		set(hObject,'Value',1),		return,		end
	set(handles.radio_newWin,'Val',0)

% -----------------------------------------------------------------------------------------
function radio_newWin_Callback(hObject, eventdata, handles)
	if (~get(hObject,'Value')),		set(hObject,'Value',1),		return,		end
	set(handles.radio_inPlace,'Val',0)

% % -----------------------------------------------------------------------------------------
% function push_help_Callback(hObject, eventdata, handles)

% -----------------------------------------------------------------------------------------
function push_OK_Callback(hObject, eventdata, handles)

	[X,Y,Z,head,rows,cols] = load_grd(handles.handMir);
	if isempty(Z),		return,		end		% An error message was already issued
	if (get(handles.radio_newWin,'Val'))
		hdr.X = X;		hdr.Y = Y;		hdr.head = head;		% Save this for use at this function's end
	end

	bw = isnan(Z);

	if (get(handles.radio_paintSmall,'Val'))		% Retain only <= handles.nCells sized of connected groups
		bw2 = img_fun('bwareaopen', bw, handles.nCells);
		bw = xor(bw, bw2);
		clear bw2;
	end
	
	B = img_fun('find_holes',bw);
	
	if (get(handles.radio_surface,'Val'))
		opt_I = sprintf('-I%.10f/%.10f',head(8),head(9));
	else
		opt_I = ' ';
	end

	aguentabar(0,'title','Tapa buracos (Filling holes)','CreateCancelBtn')

	if (get(handles.radio_paintSmall,'Val'))
		pad = 4;		n_buracos = numel(B);
		for (i = 1:n_buracos)
			% Get rectangles arround each hole
			x_min = min(B{i}(:,2));			x_max = max(B{i}(:,2));
			y_min = min(B{i}(:,1));			y_max = max(B{i}(:,1));
			x_min = max(1,x_min-pad);		x_max = min(x_max+pad,cols);
			y_min = max(1,y_min-pad);		y_max = min(y_max+pad,rows);
			x_min = head(1) + (x_min-1)*head(8);    x_max = head(1) + (x_max-1)*head(8);
			y_min = head(3) + (y_min-1)*head(9);    y_max = head(3) + (y_max-1)*head(9);

			rect_crop = [x_min y_min (x_max-x_min) (y_max-y_min)];
			[Z_rect, r_c]  = cropimg(head(1:2),head(3:4),Z,rect_crop,'out_grid');
			bw_rect = cropimg(head(1:2),head(3:4),bw,rect_crop,'out_grid');
			Z_rect = double(Z_rect);      % It has to be (GHRRRRRRRRRRRRR)

			X = x_min:head(8):x_max;	Y = y_min:head(9):y_max;
			[XX,YY] = meshgrid(X,Y);
			XX(bw_rect) = [];			YY(bw_rect) = [];		Z_rect(bw_rect) = [];

			if (get(handles.radio_surface,'Val'))
				opt_R = sprintf('-R%.10f/%.10f/%.10f/%.10f', X(1), X(end), Y(1), Y(end));
				Z_rect = surface_m( XX(:), YY(:), Z_rect(:), opt_R, opt_I, '-T.25' );
			elseif (get(handles.radio_bicubic,'Val'))
				Z_rect = griddata_j(XX(:), YY(:), Z_rect(:), X, Y, 'cubic');
			else
				Z_rect = griddata_j(XX(:), YY(:), Z_rect(:), X, Y, 'linear');
			end

			% Inprint the processed rectangle back into orig array
			if (isa(Z,'single')),		Z(r_c(1):r_c(2),r_c(3):r_c(4)) = single(Z_rect);
			elseif (isa(Z,'int16')),	Z(r_c(1):r_c(2),r_c(3):r_c(4)) = int16(Z_rect);
			elseif (isa(Z,'uint16')),	Z(r_c(1):r_c(2),r_c(3):r_c(4)) = uint16(Z_rect);
			else						Z(r_c(1):r_c(2),r_c(3):r_c(4)) = single(Z_rect);
			end

			h = aguentabar(i/n_buracos);
			if (isnan(h)),	break,	end
		end
		if (isnan(h)),	return,		end			% User hit cancel button

	else
		% Fill all holes. We'll do that by tilling to speed up and lower memory consumption. (Side effects ?)
		[Z,h] = interpbytiles(handles, Z, head, rows, cols, opt_I, 3);
		if (isnan(h)),	return,		end			% User hit cancel button
		
		have_nans = grdutils(Z,'-N');		% See if we still have NaNs (in almost cases we shouldn't)
		if (have_nans)
			aguentabar(0,'title','Second Run (Filling holes)','CreateCancelBtn')
			[Z,h] = interpbytiles(handles, Z, head, rows, cols, opt_I, 2);		% Recomputing with larger tiles
			if (isnan(h)),	return,		end			% User hit cancel button
			
			% One last check. I hope it will never be true but ...
			have_nans = grdutils(Z,'-N');
			if (have_nans)
				aguentabar(0,'title','Thirth Run (Filling holes)','CreateCancelBtn')
				[Z,h] = interpbytiles(handles, Z, head, rows, cols, opt_I, 1);	% Recomputing with whole grid (shit)
				if (isnan(h)),	return,		end			% User hit cancel button
			end
		end

	end

	zz = grdutils(Z,'-L');		z_min = zz(1);		z_max = zz(2);

	if ( get(handles.radio_newWin,'Val') )		% Output in new window
		hdr.geog = handles.handMir.geog;	hdr.name = 'Filled holes';
		hdr.head(5:6) = [z_min z_max];
		mirone(Z, hdr)
		figure(handles.figure1)		% Don't let it hiden behind the Mirone figure
	else										% Update current figure
		img = scaleto8(Z);
		set(handles.handMir.hImg, 'CData', img)
		setappdata(handles.handMir.figure1,'dem_z',Z);		% Update grid so that coursor display correct values
		have_nans = grdutils(Z,'-N');
		if (handles.handMir.Illumin_type > 0 || ~have_nans)
			handles.handMir.Illumin_type == 0;
			handles.handMir.have_nans = have_nans;
			guidata(handles.handMir.figure1, handles.handMir)	% Update Mirone handles
		end
		delete(handles.figure1)		% Force a restart if user wants more. Otherwise what should be the contents of "Z" ?
	end

% -----------------------------------------------------------------------------------------
function [Z,h] = interpbytiles(handles, Z, head, rows, cols, opt_I, n_tiles)
% Interpolate the holes by tilling the Z array into n_tiles x n_tiles
	um_part_row = fix(rows/n_tiles);		um_part_col = fix(cols/n_tiles);
	pad = 10;
	for (m = 1:n_tiles)			% loop over 1/n_tiles rows
		m_this = [max(1, um_part_row * (m-1) - pad) 	min(um_part_row * m + pad, rows)];
		for (n = 1:n_tiles)		% loop over 1/n_tiles cols
			n_this = [max(1, um_part_col * (n-1) - pad) 	min(um_part_col * n + pad, cols)];
			Z_rect = Z(m_this(1):m_this(2), n_this(1):n_this(2));
			have_nans = grdutils(Z_rect,'-N');
			if ( ~have_nans ),		continue,	end		% This tile has no NaNs so goto next one
			Z_rect = double(Z_rect);			% It has to be double (GHRRRRRRRRRRRRR)
			bw = isnan(Z_rect);			% Get the NaNs mask

			x_min = head(1) + (n_this(1)-1)*head(8);    x_max = head(1) + (n_this(2)-1)*head(8);
			y_min = head(3) + (m_this(1)-1)*head(9);    y_max = head(3) + (m_this(2)-1)*head(9);

			X = x_min:head(8):x_max;	Y = y_min:head(9):y_max;
			[XX,YY] = meshgrid(X,Y);
			XX(bw) = [];			YY(bw) = [];		Z_rect(bw) = [];
			
			if (get(handles.radio_surface,'Val'))
				opt_R = sprintf('-R%.10f/%.10f/%.10f/%.10f', X(1), X(end), Y(1), Y(end));
				Z_rect = surface_m( XX(:), YY(:), Z_rect(:), opt_R, opt_I, '-T.25' );
			elseif (get(handles.radio_bicubic,'Val'))
				Z_rect = griddata_j(XX(:), YY(:), Z_rect(:), X, Y, 'cubic');
			else
				Z_rect = griddata_j(XX(:), YY(:), Z_rect(:), X, Y, 'linear');
			end

			% Inprint the processed rectangle back into orig array
			if (isa(Z,'single')),		Z( m_this(1):m_this(2),n_this(1):n_this(2) ) = single(Z_rect);
			elseif (isa(Z,'int16')),	Z( m_this(1):m_this(2),n_this(1):n_this(2) ) = int16(Z_rect);
			elseif (isa(Z,'uint16')),	Z( m_this(1):m_this(2),n_this(1):n_this(2) ) = uint16(Z_rect);
			else						Z( m_this(1):m_this(2),n_this(1):n_this(2) ) = single(Z_rect);
			end

			h = aguentabar((n*m)/9);		% ...
			if (isnan(h)),	break,	end
		end
	end

% -----------------------------------------------------------------------------------------
% --- Creates and returns a handle to the GUI figure. 
function inpaint_nans_LayoutFcn(h1);

set(h1, 'Position',[520 696 271 120],...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'MenuBar','none',...
'Name','Inpaint NaNs',...
'NumberTitle','off',...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1');

uicontrol('Parent',h1,'Position',[110 29 81 85],'Style','frame','Tag','frame2');
uicontrol('Parent',h1,'Position',[10 29 91 85],'Style','frame','Tag','frame1');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'HorizontalAlignment','left',...
'Position',[15 41 35 15],...
'String','N cells',...
'Style','text',...
'Tag','text1');

uicontrol('Parent',h1,...
'Callback',{@inpaint_nans_uicallback,h1,'radio_paintAll_Callback'},...
'FontName','Helvetica',...
'Position',[20 86 60 15],...
'String','Paint all',...
'Style','radiobutton',...
'TooltipString','Fill all NaN zones ',...
'Tag','radio_paintAll');

uicontrol('Parent',h1,...
'Callback',{@inpaint_nans_uicallback,h1,'radio_paintSmall_Callback'},...
'FontName','Helvetica',...
'Position',[20 64 75 15],...
'String','Paint small',...
'Style','radiobutton',...
'TooltipString','Fill NaN zones smaler than N cells',...
'Value',1,...
'Tag','radio_paintSmall');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@inpaint_nans_uicallback,h1,'edit_nCells_Callback'},...
'Position',[48 37 45 21],...
'String','10',...
'Style','edit',...
'TooltipString','Interpolate only NaN regions with at maximum this number of grid cells',...
'Tag','edit_nCells');

uicontrol('Parent',h1,...
'Callback',{@inpaint_nans_uicallback,h1,'radio_surface_Callback'},...
'FontName','Helvetica',...
'Position',[123 86 60 15],...
'String','Surface',...
'Style','radiobutton',...
'TooltipString','Minimum curvature interpolation',...
'Value',1,...
'Tag','radio_surface');

uicontrol('Parent',h1,...
'Callback',{@inpaint_nans_uicallback,h1,'radio_bicubic_Callback'},...
'FontName','Helvetica',...
'Position',[123 64 60 15],...
'String','Bicubic',...
'Style','radiobutton',...
'TooltipString','Bicubic interpolation',...
'Tag','radio_bicubic');

uicontrol('Parent',h1,...
'Callback',{@inpaint_nans_uicallback,h1,'radio_bilinear_Callback'},...
'FontName','Helvetica',...
'Position',[123 41 60 15],...
'String','Bilinear',...
'Style','radiobutton',...
'TooltipString','Bilinear interpolation',...
'Tag','radio_bilinear');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'Position',[127 106 45 15],...
'String','Method',...
'Style','text',...
'Tag','text_method');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'Position',[23 106 60 15],...
'String','What to do',...
'Style','text',...
'Tag','text_what');

% uicontrol('Parent',h1,...
% 'Callback',{@inpaint_nans_uicallback,h1,'push_help_Callback'},...
% 'FontName','Helvetica',...
% 'FontSize',11,...
% 'FontWeight','bold',...
% 'ForegroundColor',[0 0 1],...
% 'Position',[222 84 23 23],...
% 'String','?',...
% 'Tag','push_help');

uicontrol('Parent',h1,...
'Callback',{@inpaint_nans_uicallback,h1,'radio_inPlace_Callback'},...
'FontName','Helvetica',...
'Position',[118 6 60 15],...
'String','In place',...
'Style','radiobutton',...
'TooltipString','Apply changes to current window. If it screws you have to restart',...
'Tag','radio_inPlace');

uicontrol('Parent',h1,...
'Callback',{@inpaint_nans_uicallback,h1,'radio_newWin_Callback'},...
'FontName','Helvetica',...
'Position',[11 6 86 16],...
'String','New window',...
'Style','radiobutton',...
'TooltipString','Put the results in a new Mirone window',...
'Value',1,...
'Tag','radio_newWin');

uicontrol('Parent',h1,...
'Callback',{@inpaint_nans_uicallback,h1,'push_OK_Callback'},...
'FontName','Helvetica',...
'FontSize',9,...
'Position',[200 38 66 23],...
'String','Compute',...
'Tag','push_OK');

uicontrol('Parent',h1,...
'Callback','delete(gcf)',...
'FontName','Helvetica',...
'FontSize',9,...
'Position',[200 8 66 23],...
'String','Cancel',...
'Tag','push_cancel');


function inpaint_nans_uicallback(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
feval(callback_name,hObject,[],guidata(h1));

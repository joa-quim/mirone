function varargout = grdsample_mir(varargin)
% Helper window to interface with grdsample MEX 

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

% $Id: grdsample_mir.m 9860 2016-10-11 23:39:58Z j $

	if isempty(varargin)
		errordlg('GRDSAMPLE: wrong number of input arguments.','Error'),	return
	end
 
	hObject = figure('Vis','off');
	grdsample_mir_LayoutFcn(hObject);
	handles = guihandles(hObject);

	handMir  = varargin{1};
	handles.Z = getappdata(handMir.figure1,'dem_z');
	move2side(handMir.figure1, hObject,'right')
    
	if (handMir.no_file)
		errordlg('GRDSAMPLE: You didn''t even load a file. What are you expecting then?','ERROR')
        delete(hObject);    return
	end
	if (~handMir.validGrid)
        errordlg('GRDSAMPLE: This operation is deffined only for images derived from DEM grids.','ERROR')
        delete(hObject);    return
	end
	if (isempty(handles.Z))
        errordlg('GRDSAMPLE: Grid was not saved in memory. Increase "Grid max size" and start over.','ERROR')
        delete(hObject);    return
	end

	handles.x_min = [];			handles.x_max = [];
	handles.y_min = [];			handles.y_max = [];
	handles.x_inc = [];			handles.y_inc = [];
	handles.dms_xinc = 0;		handles.dms_yinc = 0;
	handles.bd_cond = [];
    handles.nr_or = size(handles.Z,1);
    handles.nc_or = size(handles.Z,2);
    handles.hMirFig = handMir.figure1;

	%-----------
	% Fill in the grid limits boxes with calling fig values and save some limiting value
    head = handMir.head;
	set(handles.edit_x_min,'String',sprintf('%.8g',head(1)))
	set(handles.edit_x_max,'String',sprintf('%.8g',head(2)))
	set(handles.edit_y_min,'String',sprintf('%.8g',head(3)))
	set(handles.edit_y_max,'String',sprintf('%.8g',head(4)))
	handles.x_min = head(1);            handles.x_max = head(2);
	handles.y_min = head(3);            handles.y_max = head(4);
	handles.x_min_or = head(1);         handles.x_max_or = head(2);
	handles.y_min_or = head(3);         handles.y_max_or = head(4);
	handles.one_or_zero = ~head(7);

	% Fill in the x,y_inc and nrow,ncol boxes
	set(handles.edit_Nrows,'String',sprintf('%d',handles.nr_or))
	set(handles.edit_Ncols,'String',sprintf('%d',handles.nc_or))
	set(handles.edit_y_inc,'String',sprintf('%.10g',head(9)))
	set(handles.edit_x_inc,'String',sprintf('%.10g',head(8)))
	handles.x_inc = head(8);    handles.y_inc = head(9);
	handles.head = head;
	%----------------

	%------------ Give a Pro look (3D) to the frame boxes  -------------------------------
	new_frame3D(hObject, handles.GLG)
	%------------- END Pro look (3D) -----------------------------------------------------

	% Add this figure handle to the carraças list
	plugedWin = getappdata(handles.hMirFig,'dependentFigs');
	plugedWin = [plugedWin hObject];
	setappdata(handles.hMirFig,'dependentFigs',plugedWin);

	guidata(hObject, handles);

	set(hObject,'Visible','on');
	if (nargout),   varargout{1} = hObject;     end

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

%--------------------------------------------------------------------------------
function edit_refGrid_CB(hObject, handles)
	fname = get(hObject,'String');
	if (isempty(fname)),	return,		end
	push_refGrid_CB(handles.push_refGrid, handles, fname)

%--------------------------------------------------------------------------------
function push_refGrid_CB(hObject, handles, opt)
	if (nargin == 2)
		str1 = {'*.grd;*.GRD', 'Grid files (*.grd,*.GRD)';'*.*', 'All Files (*.*)'};
		[FileName, PathName, handles] = put_or_get_file(handles,str1,'Select grid','get');
		if isequal(FileName,0),		return,		end
		fname = [PathName FileName];
		set(handles.edit_refGrid,'Str',fname)
	else
		fname = opt;
	end

	% Now fill the edit boxes after the -R -I of this grid
	hdr = c_grdinfo(fname,'no_struct');
	hdr = hdr.data;
	set(handles.edit_x_min,'Str', hdr(1)),		set(handles.edit_x_max,'Str', hdr(2))
	set(handles.edit_y_min,'Str', hdr(3)),		set(handles.edit_y_max,'Str', hdr(4))
	set(handles.edit_x_inc,'Str', hdr(7)),		set(handles.edit_y_inc,'Str', hdr(8))
	set(handles.edit_Ncols,'Str', hdr(9)),		set(handles.edit_Nrows,'Str', hdr(10))

% --------------------------------------------------------------------
function push_Help_R_F_T_CB(hObject, handles)
	message = {'Min and Max, of "X Direction" and "Y Direction" specify the Region of'
        'interest. To specify boundaries in degrees and minutes [and seconds],'
        'use the dd:mm[:ss.xx] format.'
        '"Spacing" sets the grid size for grid output. You may choose different'
        'spacings for X and Y. Also here you can use the dd:mm[:ss.xx] format.'
        'In "#of lines" it is offered the easyeast way of controling the grid'
        'dimensions (lines & columns).'};
	helpdlg(message,'Help on Grid Line Geometry');

% --------------------------------------------------------------------
function popup_BoundaryCondition_CB(hObject, handles)
	val = get(hObject,'Value');     str = get(hObject, 'String');
	switch str{val};
        case ' ',        handles.bd_cond = [];
        case '',         handles.bd_cond = [];
        case 'x',        handles.bd_cond = '-Lx';
        case 'y',        handles.bd_cond = '-Ly';
        case 'xy',       handles.bd_cond = '-Lxy'; 
        case 'g',        handles.bd_cond = '-Lg';
	end
	guidata(hObject, handles);

% --------------------------------------------------------------------
function push_Help_L_CB(hObject, handles)
	message = {'Boundary condition flag may be "x" or "y" or "xy" indicating data is periodic'
               'in range of x or y or both set by the grids limits in the above boxes,'
               'or flag may be "g" indicating geographical conditions (x and y may be'
               'lon and lat). [Default is no boundary conditions].'};
	helpdlg(message,'Help -L option');

% --------------------------------------------------------------------
function push_OK_CB(hObject, handles)
	
	global gmt_ver

	opt_R = ' ';     opt_N = ' ';     opt_Q = ' ';     opt_L = ' ';
    n_set = 0;
	x_min = get(handles.edit_x_min,'String');   x_max = get(handles.edit_x_max,'String');
	y_min = get(handles.edit_y_min,'String');   y_max = get(handles.edit_y_max,'String');
	if isempty(x_min) || isempty(x_max) || isempty(y_min) || isempty(y_max)
        errordlg('One or more grid limits are empty. Try with your yes open.','Error');    return
	end

	nx = str2double(get(handles.edit_Ncols,'String'));
	ny = str2double(get(handles.edit_Nrows,'String'));
	if (isnan(nx) || isnan(ny))      % I think this was already tested, but ...
        errordlg('One (or two) of the grid dimensions are not valid. Do your best.','Error');   return
	end

	if (nx ~= handles.nc_or || ny ~= handles.nr_or)
        opt_N = ['-N' get(handles.edit_Ncols,'String') '/' get(handles.edit_Nrows,'String')];
        n_set = 1;
	end

	if (get(handles.checkbox_Option_Q,'Value')),    opt_Q = '-Q';   end
	if (~isempty(handles.bd_cond)),     opt_L = handles.bd_cond;    end

	% See if grid limits were changed
	if ((abs(handles.x_min-handles.x_min_or) > 1e-5) || (abs(handles.x_max-handles.x_max_or) > 1e-5) || ...
        (abs(handles.y_min-handles.y_min_or) > 1e-5) || (abs(handles.y_max-handles.y_max_or) > 1e-5))

		opt_R = sprintf('-R%.12g/%.12g/%.12g/%.12g',handles.x_min, handles.x_max, handles.y_min, handles.y_max);
        x_min = handles.x_min;    x_max = handles.x_max;
        y_min = handles.y_min;    y_max = handles.y_max;
        if (~n_set)     % Only limits had changed, but we need also to return the -N option
            opt_N = ['-N' get(handles.edit_Ncols,'String') '/' get(handles.edit_Nrows,'String')];
            n_set = 1;
        end
	end

	if (~n_set)
        errordlg('You haven''t select anything usefull to do.','Chico Clever');   return
	end

	if (gmt_ver == 5)				% Than use the GMT5 syntax
		opt_N(2) = 'I';
		opt_N = strrep(opt_N, '/', '+/');
		opt_N(end+1) = '+';
	end

	set(handles.figure1,'pointer','watch');	    set(handles.hMirFig,'pointer','watch')
	newZ = c_grdsample(handles.Z, handles.head, opt_R, opt_N, opt_Q, opt_L);
	zMinMax = grdutils(newZ,'-L');	    [ny,nx] = size(newZ);
	set(handles.figure1,'pointer','arrow');	    set(handles.hMirFig,'pointer','arrow')
    new_head = [handles.head(1:4) zMinMax(1:2)' handles.head(7)];
	if (~strcmp(opt_R,' '))       % Grid limits did change
        new_head(1:4) = [x_min x_max y_min y_max];    
	end
    x_inc = (new_head(2) - new_head(1)) / (nx - ~handles.head(7));
    y_inc = (new_head(4) - new_head(3)) / (ny - ~handles.head(7));
    new_head(8:9) = [x_inc y_inc];
    tmp.X = linspace(new_head(1),new_head(2),nx);       tmp.Y = linspace(new_head(3),new_head(4),ny);
    tmp.head = new_head;
    tmp.name = 'Resampled grid';
	prjInfoStruc = aux_funs('getFigProjInfo', handles);
	if (~isempty(prjInfoStruc.projWKT))
		tmp.srsWKT = prjInfoStruc.projWKT;
	elseif (~isempty(prjInfoStruc.proj4))
		tmp.srsWKT = prjInfoStruc.projWKT;
	end
    mirone(newZ,tmp);
	if (~isempty(prjInfoStruc.projGMT))		% If we have this rare one it was not applied above so do it now.
		setappdata(handles.figure1, 'ProjGMT', prjInfoStruc.projGMT);
	end
    figure(handles.figure1)         % Don't let this figure forgotten behind the newly created one

% --------------------------------------------------------------------
% --- Executes on key press over figure1 with no controls selected.
function figure1_KeyPressFcn(hObject, handles)
	if isequal(get(hObject,'CurrentKey'),'escape')      % Check for "escape"
		delete(handles.figure1);
	end

% --- Creates and returns a handle to the GUI figure. 
function grdsample_mir_LayoutFcn(h1)

yoff = 30;
set(h1, 'Position',[266 399 441 158+yoff],...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'KeyPressFcn',@figure1_KeyPressFcn,...
'MenuBar','none',...
'Name','Grdsample',...
'NumberTitle','off',...
'RendererMode','manual',...
'Resize','off',...
'Tag','figure1');

uicontrol('Parent',h1,'Enable','inactive','Position',[10 75+yoff 421 75],'Style','frame');

uicontrol('Parent',h1, 'Position',[30 142+yoff 125 15],...
'Enable','inactive',...
'String','Griding Line Geometry',...
'Style','text',...
'Tag','GLG');

uicontrol('Parent',h1, 'Position',[77 111+yoff 80 21],...
'BackgroundColor',[1 1 1],...
'Call',@grdsample_mir_uiCB,...
'HorizontalAlignment','left',...
'Style','edit',...
'Tag','edit_x_min');

uicontrol('Parent',h1, 'Position',[163 111+yoff 80 21],...
'BackgroundColor',[1 1 1],...
'Call',@grdsample_mir_uiCB,...
'HorizontalAlignment','left',...
'Style','edit',...
'Tag','edit_x_max');

uicontrol('Parent',h1, 'Position',[77 85+yoff 80 21],...
'BackgroundColor',[1 1 1],...
'Call',@grdsample_mir_uiCB,...
'HorizontalAlignment','left',...
'Style','edit',...
'Tag','edit_y_min');

uicontrol('Parent',h1, 'Position',[163 85+yoff 80 21],...
'BackgroundColor',[1 1 1],...
'Call',@grdsample_mir_uiCB,...
'HorizontalAlignment','left',...
'Style','edit',...
'Tag','edit_y_max');

uicontrol('Parent',h1, 'Position',[22 115+yoff 55 15],...
'Enable','inactive',...
'HorizontalAlignment','left',...
'String','X Direction',...
'Style','text');

uicontrol('Parent',h1, 'Position',[21 89+yoff 55 15],...
'Enable','inactive',...
'HorizontalAlignment','left',...
'String','Y Direction',...
'Style','text',...
'Tag','text3');

uicontrol('Parent',h1, 'Position',[183 132+yoff 41 13],...
'Enable','inactive',...
'String','Max',...
'Style','text');

uicontrol('Parent',h1, 'Position',[99 132+yoff 41 13],...
'Enable','inactive',...
'String','Min',...
'Style','text');

uicontrol('Parent',h1, 'Position',[248 111+yoff 71 21],...
'BackgroundColor',[1 1 1],...
'Call',@grdsample_mir_uiCB,...
'HorizontalAlignment','left',...
'Style','edit',...
'TooltipString','DX grid spacing',...
'Tag','edit_x_inc');

uicontrol('Parent',h1, 'Position',[248 85+yoff 71 21],...
'BackgroundColor',[1 1 1],...
'Call',@grdsample_mir_uiCB,...
'HorizontalAlignment','left',...
'Style','edit',...
'TooltipString','DY grid spacing',...
'Tag','edit_y_inc');

uicontrol('Parent',h1, 'Position',[324 111+yoff 65 21],...
'BackgroundColor',[1 1 1],...
'Call',@grdsample_mir_uiCB,...
'HorizontalAlignment','center',...
'Style','edit',...
'TooltipString','Number of columns in the grid',...
'Tag','edit_Ncols');

uicontrol('Parent',h1, 'Position',[324 85+yoff 65 21],...
'BackgroundColor',[1 1 1],...
'Call',@grdsample_mir_uiCB,...
'HorizontalAlignment','center',...
'Style','edit',...
'Tooltip','Number of columns in the grid',...
'Tag','edit_Nrows');

uicontrol('Parent',h1, 'Position',[265 133+yoff 41 13],...
'Enable','inactive',...
'String','Spacing',...
'Style','text');

uicontrol('Parent',h1, 'Position',[332 133+yoff 51 13],...
'Enable','inactive',...
'String','# of lines',...
'Style','text');

uicontrol('Parent',h1, 'Position',[400 84+yoff 21 48],...
'BackgroundColor',[0.8313725591 0.815686285495758 0.7843137383461],...
'Call',@grdsample_mir_uiCB,...
'FontWeight','bold',...
'ForegroundColor',[0 0 1],...
'String','?',...
'Tag','push_Help_R_F_T');

%...
uicontrol('Parent',h1,'Enable','inactive','Position',[10 30+yoff 421 40],'Style','frame');

uicontrol('Parent',h1,'Position',[17 42+yoff 60 16],...
'HorizontalAlignment','left',...
'String','OR Ref grid',...
'Style','text');

uicontrol('Parent',h1, 'Position',[77 40+yoff 289 21],...
'Call',@grdsample_mir_uiCB,...
'BackgroundColor',[1 1 1],...
'HorizontalAlignment','left',...
'Tooltip','Alternatively, give a grid name to use its Region and Increment',...
'Style','edit',...
'Tag','edit_refGrid');

uicontrol('Parent',h1, 'Position',[365 39+yoff 23 23],...
'Call',@grdsample_mir_uiCB,...
'String','...', 'FontWeight','bold',...
'Tag','push_refGrid');

uicontrol('Parent',h1,'Position',[17 31 95 16],...
'Enable','inactive',...
'HorizontalAlignment','left',...
'String','Boundary condition',...
'Style','text');

uicontrol('Parent',h1, 'Position',[114 29 47 21],...
'BackgroundColor',[1 1 1],...
'Call',@grdsample_mir_uiCB,...
'HorizontalAlignment','right',...
'String',{' '; 'x'; 'y'; 'xy'; 'g' },...
'Style','popupmenu',...
'Value',1,...
'Tag','popup_BoundaryCondition');

uicontrol('Parent',h1, 'Position',[172 28 21 23],...
'Call',@grdsample_mir_uiCB,...
'FontWeight','bold',...
'ForegroundColor',[0 0 1],...
'String','?',...
'Tag','push_Help_L');

uicontrol('Parent',h1, 'Position',[221 32 140 15],...
'String','Bilinear interpolation',...
'Style','checkbox',...
'Tooltip','Use bilinear rather than bicubic interpolation',...
'Tag','checkbox_Option_Q');

uicontrol('Parent',h1, 'Position',[365 8 66 21],...
'Call',@grdsample_mir_uiCB,...
'String','OK',...
'Tag','push_OK');

function grdsample_mir_uiCB(hObject, eventdata)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));

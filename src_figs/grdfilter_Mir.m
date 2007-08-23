function varargout = grdfilter_Mir(varargin)
% M-File changed by desGUIDE 
% varargin   command line arguments to grdfilter_Mir (see VARARGIN) 

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
 
	hObject = figure('Tag','figure1','Visible','off');
	handles = guihandles(hObject);
	guidata(hObject, handles);
	grdfilter_Mir_LayoutFcn(hObject,handles);
	handles = guihandles(hObject);
	movegui(hObject,'center')

	if ~isempty(varargin)
		handMir  = varargin{1};
		handles.Z = getappdata(handMir.figure1,'dem_z');
	else
        errordlg('GRDSAMPLE: wrong number of arguments.','Error')
        delete(hObject);    return
	end

	if (~handles.ValidGrid)
        errordlg('GRDSAMPLE: This operation is deffined only for images derived from DEM grids.','ERROR')
        delete(hObject);    return
	end
	if (isempty(handles.Z))
        errordlg('GRDSAMPLE: Grid was not saved in memory. Increase "Grid max size" and start over.','ERROR')
        delete(hObject);    return
	end

	handles.command = cell(5,1);
	handles.command{1} = '-Fb';
	if (handMir.geog)
		handles.command{4} = '-D1';
		set(handles.popup_Option_D,'Value',2)
	else
		handles.command{4} = '-D0';
	end

	handles.x_min = [];             handles.x_max = [];
	handles.y_min = [];             handles.y_max = [];
	handles.x_inc = [];             handles.y_inc = [];
	handles.dms_xinc = 0;           handles.dms_yinc = 0;

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
	handles.one_or_zero = head(7);
	handles.head = head;
	[m,n] = size(handles.Z);
	handles.nr_or = m;                  handles.nc_or = n;;

	% Fill in the x,y_inc and nrow,ncol boxes
	set(handles.edit_Nrows,'String',sprintf('%d',handles.nr_or))
	set(handles.edit_Ncols,'String',sprintf('%d',handles.nc_or))
	set(handles.edit_y_inc,'String',sprintf('%.10g',head(9)))
	set(handles.edit_x_inc,'String',sprintf('%.10g',head(8)))
	handles.x_inc = head(8);        handles.y_inc = head(9);
	handles.x_inc_or = head(8);     handles.y_inc_or = head(9);
	%----------------

	% Give a Pro look (3D) to the frame boxes 
	bgcolor = get(0,'DefaultUicontrolBackgroundColor');
	framecolor = max(min(0.65*bgcolor,[1 1 1]),[0 0 0]);
	h_f = findobj(hObject,'Style','Frame');
	frame_size = get(h_f,'Position');
	f_bgc = get(h_f,'BackgroundColor');
	usr_d = get(h_f,'UserData');
	if abs(f_bgc(1)-bgcolor(1)) > 0.01           % When the frame's background color is not the default's
		frame3D(hObject,frame_size,framecolor,f_bgc,usr_d)
	else
		frame3D(hObject,frame_size,framecolor,'',usr_d)
		delete(h_f)
	end

	% Recopy the text fields on top of previously created frames (uistack is to slow)
	h_t = handles.text8;
	usr_d = get(h_t,'UserData');
	t_size = get(h_t,'Position');   t_str = get(h_t,'String');    fw = get(h_t,'FontWeight');
	bgc = get (h_t,'BackgroundColor');   fgc = get (h_t,'ForegroundColor');
	uicontrol('Parent',hObject, 'Style','text', 'Position',t_size,'String',t_str, ...
		'BackgroundColor',bgc,'ForegroundColor',fgc,'FontWeight',fw,'UserData',usr_d);
	delete(h_t)
	
	guidata(hObject, handles);
	
	set(hObject,'Visible','on');
	if (nargout),   varargout{1} = hObject;     end

% -------------------------------------------------------------------------------------
function edit_x_min_Callback(hObject, eventdata, handles)
x_min_or = handles.x_min_or;
xx = get(hObject,'String');     val = test_dms(xx);
if ~isempty(val)            % when dd:mm or dd:mm:ss was given
    x_min = 0;
    if str2double(val{1}) > 0
        for i = 1:length(val)   x_min = x_min + str2double(val{i}) / (60^(i-1));    end
    else
        for i = 1:length(val)   x_min = x_min - abs(str2double(val{i})) / (60^(i-1));   end
    end
    if x_min < x_min_or; set(hObject,'String',x_min_or);    return;     end; 
    handles.x_min = x_min;
    if ~isempty(handles.x_max) & x_min >= handles.x_max
        errordlg('West Longitude >= East Longitude ','Error in Longitude limits')
        set(hObject,'String','');   guidata(hObject, handles);  return
    else
        guidata(hObject, handles);
    end
    nc = get(handles.edit_Ncols,'String');
    if ~isempty(handles.x_max) & ~isempty(nc)       % x_max and ncols boxes are filled
        % Compute Ncols, but first must recompute x_inc
        x_inc = ivan_the_terrible((handles.x_max - x_min),round(abs(str2double(nc))),1);
        xx = floor((handles.x_max - str2double(xx)) / (str2double(get(handles.edit_x_inc,'String')))+0.5) + handles.one_or_zero;
        set(handles.edit_x_inc,'String',num2str(x_inc,8))
        guidata(hObject, handles);
    elseif ~isempty(handles.x_max)      % x_max box is filled but ncol is not, so put to the default (100)
        x_inc = ivan_the_terrible((handles.x_max - x_min),100,1);
        set(handles.edit_x_inc,'String',num2str(x_inc,8))
        set(handles.edit_Ncols,'String','100')
        guidata(hObject, handles);
    end
else                % box is empty, so clear also x_inc and ncols
    set(handles.edit_x_inc,'String','');     set(handles.edit_Ncols,'String','');
    set(hObject,'String','');   guidata(hObject, handles);
end

% -------------------------------------------------------------------------------------
function edit_x_max_Callback(hObject, eventdata, handles)
x_max_or = handles.x_max_or;
xx = get(hObject,'String');     val = test_dms(xx);
if ~isempty(val)
    x_max = 0;
    if str2double(val{1}) > 0
        for i = 1:length(val)   x_max = x_max + str2double(val{i}) / (60^(i-1));    end
    else
        for i = 1:length(val)   x_max = x_max - abs(str2double(val{i})) / (60^(i-1));   end
    end
    if x_max > x_max_or; set(hObject,'String',x_max_or);    return;     end; 
    handles.x_max = x_max;
    if ~isempty(handles.x_min) & x_max <= handles.x_min 
        errordlg('East Longitude <= West Longitude','Error in Longitude limits')
        set(hObject,'String','');   guidata(hObject, handles);  return
    else
        guidata(hObject, handles);
    end
    nc = get(handles.edit_Ncols,'String');
    if ~isempty(handles.x_min) & ~isempty(nc)       % x_max and ncols boxes are filled
        % Compute Ncols, but first must recompute x_inc
        x_inc = ivan_the_terrible((x_max - handles.x_min),round(abs(str2double(nc))),1);
        xx = floor((handles.x_min - str2double(xx)) / (str2double(get(handles.edit_x_inc,'String')))+0.5) + handles.one_or_zero;
        set(handles.edit_x_inc,'String',num2str(x_inc,8))
        guidata(hObject, handles);    
    elseif ~isempty(handles.x_min)      % x_min box is filled but ncol is not, so put to the default (100)
        x_inc = ivan_the_terrible((x_max - handles.x_min),100,1);
        set(handles.edit_x_inc,'String',num2str(x_inc,8))
        set(handles.edit_Ncols,'String','100')
        guidata(hObject, handles);
    end
else                % box is empty, so clear also x_inc and ncols
    set(handles.edit_x_inc,'String','');     set(handles.edit_Ncols,'String','');
    set(hObject,'String','');   guidata(hObject, handles);
end

% -------------------------------------------------------------------------------------
function edit_y_min_Callback(hObject, eventdata, handles)
% Read value either in decimal or in the dd:mm or dd_mm:ss formats and do some tests
y_min_or = handles.y_min_or;
xx = get(hObject,'String');     val = test_dms(xx);
if ~isempty(val)
    y_min = 0;
    if str2double(val{1}) > 0
        for i = 1:length(val)   y_min = y_min + str2double(val{i}) / (60^(i-1));    end
    else
        for i = 1:length(val)   y_min = y_min - abs(str2double(val{i})) / (60^(i-1));   end
    end
    if y_min < y_min_or; set(hObject,'String',y_min_or);    return;     end; 
    handles.y_min = y_min;
    if ~isempty(handles.y_max) & y_min >= handles.y_max
        errordlg('South Latitude >= North Latitude','Error in Latitude limits')
        set(hObject,'String','');   guidata(hObject, handles);  return
    else
        guidata(hObject, handles);
    end
    nr = get(handles.edit_Nrows,'String');
    if ~isempty(handles.y_max) & ~isempty(nr)       % y_max and nrows boxes are filled
        % Compute Nrowss, but first must recompute y_inc
        y_inc = ivan_the_terrible((handles.y_max - y_min),round(abs(str2double(nr))),1);
        xx = floor((handles.y_max - str2double(xx)) / (str2double(get(handles.edit_y_inc,'String')))+0.5) + handles.one_or_zero;
        set(handles.edit_y_inc,'String',num2str(y_inc,8))
        guidata(hObject, handles);
    elseif ~isempty(handles.y_max)      % y_max box is filled but nrows is not, so put to the default (100)
        y_inc = ivan_the_terrible((handles.y_max - y_min),100,1);
        set(handles.edit_y_inc,'String',num2str(y_inc,8))
        set(handles.edit_Nrows,'String','100')
        guidata(hObject, handles);
    end
else                % box is empty, so clear also y_inc and nrows
    set(handles.edit_y_inc,'String','');     set(handles.edit_Nrows,'String','');
    set(hObject,'String','');   guidata(hObject, handles);
end

% -------------------------------------------------------------------------------------
function edit_y_max_Callback(hObject, eventdata, handles)
y_max_or = handles.y_max_or;
xx = get(hObject,'String');     val = test_dms(xx);
if ~isempty(val)
    y_max = 0;
    if str2double(val{1}) > 0
        for i = 1:length(val)   y_max = y_max + str2double(val{i}) / (60^(i-1));    end
    else
        for i = 1:length(val)   y_max = y_max - abs(str2double(val{i})) / (60^(i-1));   end
    end
    if y_max > y_max_or; set(hObject,'String',y_max_or);    return;     end; 
    handles.y_max = y_max;
    if ~isempty(handles.y_min) & y_max <= handles.y_min 
        errordlg('North Latitude <= South Latitude','Error in Latitude limits')
        set(hObject,'String','');   guidata(hObject, handles);  return
    else
        guidata(hObject, handles);
    end
    nr = get(handles.edit_Nrows,'String');
    if ~isempty(handles.y_min) & ~isempty(nr)       % y_min and nrows boxes are filled
        % Compute Nrows, but first must recompute y_inc
        y_inc = ivan_the_terrible((y_max - handles.y_min),round(abs(str2double(nr))),1);
        xx = floor((handles.y_min - str2double(xx)) / (str2double(get(handles.edit_y_inc,'String')))+0.5) + handles.one_or_zero;
        set(handles.edit_y_inc,'String',num2str(y_inc,8))
        guidata(hObject, handles);
    elseif ~isempty(handles.y_min)      % y_min box is filled but nrows is not, so put to the default (100)
        y_inc = ivan_the_terrible((y_max - handles.y_min),100,1);
        set(handles.edit_y_inc,'String',num2str(y_inc,8))
        set(handles.edit_Nrows,'String','100')
        guidata(hObject, handles);
    end
else                % This box is empty, so clear also y_inc and nrows
    set(handles.edit_y_inc,'String','');     set(handles.edit_Nrows,'String','');
    set(hObject,'String','');   guidata(hObject, handles);
end

% -------------------------------------------------------------------------------------
function edit_x_inc_Callback(hObject, eventdata, handles)
dms = 0;
xx = get(hObject,'String');     val = test_dms(xx);
if isempty(val)
    set(hObject, 'String', '');    return
end
% If it survived then ...
if length(val) > 1    dms = 1;      end         % inc given in dd:mm or dd:mm:ss format
x_inc = 0;
for i = 1:length(val)   x_inc = x_inc + str2double(val{i}) / (60^(i-1));    end
if ~isempty(handles.x_min) & ~isempty(handles.x_max)
    % Make whatever x_inc given compatible with GMT_grd_RI_verify
    x_inc = ivan_the_terrible((handles.x_max - handles.x_min), x_inc,2);
    if ~dms         % case of decimal unities
        set(hObject,'String',num2str(x_inc,8))
        ncol = floor((handles.x_max - handles.x_min) / x_inc + 0.5) + handles.one_or_zero;
    else            % inc was in dd:mm or dd:mm:ss format
        ncol = floor((handles.x_max - handles.x_min) / x_inc + 0.5) + handles.one_or_zero;
        ddmm = dec2deg(x_inc);
        set(hObject,'String',ddmm)
    end
    set(handles.edit_Ncols,'String',num2str(ncol))
end
handles.x_inc = x_inc;      handles.dms_xinc = dms;
guidata(hObject, handles);
if isempty(get(handles.edit_y_inc,'String'))     set(handles.edit_y_inc,'String',xx);    end

% -------------------------------------------------------------------------------------
function edit_Ncols_Callback(hObject, eventdata, handles)
xx = get(hObject,'String');
if ~isempty(get(handles.edit_x_min,'String')) & ~isempty(get(handles.edit_x_max,'String')) & ...
        ~isempty(get(handles.edit_x_inc,'String')) & ~isempty(xx)
    x_inc = ivan_the_terrible((handles.x_max - handles.x_min),round(abs(str2double(xx))),1);
    if handles.dms_xinc        % x_inc was given in dd:mm:ss format
        ddmm = dec2deg(x_inc);
        set(handles.edit_x_inc,'String',ddmm)
    else                    % x_inc was given in decimal format
        set(handles.edit_x_inc,'String',num2str(x_inc,8));
    end
    handles.x_inc = x_inc;
    guidata(hObject, handles);
end

% -------------------------------------------------------------------------------------
function edit_y_inc_Callback(hObject, eventdata, handles)
dms = 0;
xx = get(hObject,'String');     val = test_dms(xx);
if isempty(val)
    set(hObject, 'String', '');    return
end
% If it survived then ...
if length(val) > 1    dms = 1;      end         % inc given in dd:mm or dd:mm:ss format
y_inc = 0;
for i = 1:length(val)   y_inc = y_inc + str2double(val{i}) / (60^(i-1));    end
if ~isempty(handles.y_min) & ~isempty(handles.y_max)
    % Make whatever y_inc given compatible with GMT_grd_RI_verify
    y_inc = ivan_the_terrible((handles.y_max - handles.y_min), y_inc,2);
    if ~dms         % case of decimal unities
        set(hObject,'String',num2str(y_inc,8))
        nrow = floor((handles.y_max - handles.y_min) / y_inc + 0.5) + handles.one_or_zero;
    else            % inc was in dd:mm or dd:mm:ss format
        nrow = floor((handles.y_max - handles.y_min) / y_inc + 0.5) + handles.one_or_zero;
        ddmm = dec2deg(y_inc);
        set(hObject,'String',ddmm)
    end
    set(handles.edit_Nrows,'String',num2str(nrow))
end
handles.y_inc = y_inc;      handles.dms_yinc = dms;
guidata(hObject, handles);

% -------------------------------------------------------------------------------------
function edit_Nrows_Callback(hObject, eventdata, handles)
xx = get(hObject,'String');
if ~isempty(get(handles.edit_y_min,'String')) & ~isempty(get(handles.edit_y_max,'String')) & ...
        ~isempty(get(handles.edit_y_inc,'String')) & ~isempty(xx)
    y_inc = ivan_the_terrible((handles.y_max - handles.y_min),round(abs(str2double(xx))),1);
    if handles.dms_yinc        % y_inc was given in dd:mm:ss format
        ddmm = dec2deg(y_inc);
        set(handles.edit_y_inc,'String',ddmm)
    else                    % y_inc was given in decimal format
        set(handles.edit_y_inc,'String',num2str(y_inc,8));
    end
    handles.y_inc = y_inc;
    guidata(hObject, handles);
end

% -------------------------------------------------------------------------------------
function pushbutton_Help_R_F_T_Callback(hObject, eventdata, handles)
message = {'Min and Max, of "X Direction" and "Y Direction" specify the Region of'
    'interest. To specify boundaries in degrees and minutes [and seconds],'
    'use the dd:mm[:ss.xx] format.'
    '"Spacing" sets the grid size for grid output. You may choose different'
    'spacings for X and Y. Also here you can use the dd:mm[:ss.xx] format.'
    'In "#of lines" it is offered the easyeast way of controling the grid'
    'dimensions (lines & columns).'
    '"Toggle grid registration" Toggle the node registration for the output grid so'
    'as to become the opposite of the input grid [Default gives the same registration as the input grid].'};
helpdlg(message,'Help on Grid Line Geometry');

% -------------------------------------------------------------------------------------
function popup_FilterType_Callback(hObject, eventdata, handles)
	val = get(hObject,'Value');     str = get(hObject, 'String');
	switch str{val};
        case 'boxcar',              handles.command{1} = ['-Fb'];
        case 'cosine arc',          handles.command{1} = ['-Fc'];
        case 'gaussian',            handles.command{1} = ['-Fg'];
        case 'median',              handles.command{1} = ['-Fm'];
        case 'maximum likelihood',  handles.command{1} = ['-Fp'];
	end
	guidata(hObject, handles);

% -------------------------------------------------------------------------------------
function edit_FilterWidth_Callback(hObject, eventdata, handles)
	xx = get(hObject,'String');
	if ~isempty(xx)    handles.command{2} = xx;
	else    handles.command{2} = []; end
	guidata(hObject,handles)

% -------------------------------------------------------------------------------------
function pushbutton_Help_FilterType_Callback(hObject, eventdata, handles)
	message = ['Choose one only of for boxcar, cosine arch, gaussian, median, or maximum likelihood ' ...
               'probability (a mode estimator) filter and specify full width.'];
	helpdlg(message,'Help -F option');

% -------------------------------------------------------------------------------------
function popup_Option_D_Callback(hObject, eventdata, handles)
	val = get(hObject,'Value');     str = get(hObject, 'String');
	switch str{val};
        case '0',        handles.command{4} = ['-D0'];
        case '1',        handles.command{4} = ['-D1'];
        case '2',        handles.command{4} = ['-D2'];
        case '3',        handles.command{4} = ['-D3'];
        case '4',        handles.command{4} = ['-D4']; 
	end
	guidata(hObject, handles);

% -------------------------------------------------------------------------------------
function pushbutton_Help_Option_D_Callback(hObject, eventdata, handles)
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

% --------------------------------------------------------------------------------
function pushbutton_Compute_Callback(hObject, eventdata, handles)
	opt_R = ' ';    opt_I = ' ';
	x_min = get(handles.edit_x_min,'String');   x_max = get(handles.edit_x_max,'String');
	y_min = get(handles.edit_y_min,'String');   y_max = get(handles.edit_y_max,'String');
	if isempty(x_min) || isempty(x_max) || isempty(y_min) || isempty(y_max)
        errordlg('One or more grid limits are empty. Open your yes.','Error');    return
	end

	x_inc = get(handles.edit_x_inc,'String');   y_inc = get(handles.edit_y_inc,'String');
	if (isempty(x_inc) || isempty(y_inc))
        errordlg('One or two grid increments are empty. Open your yes.','Error');    return
	end

	nx = str2double(get(handles.edit_Ncols,'String'));
	ny = str2double(get(handles.edit_Nrows,'String'));
	if (isnan(nx) || isnan(ny))      % I think this was already tested, but ...
		errordlg('One (or two) of the grid dimensions are not valid. Do your best.','Error');   return
	end

	if isempty(handles.command{2})
		errordlg('Must specify a Filter width','Error');    return
	end

	% See if grid limits were changed
	if ( (abs(handles.x_min-handles.x_min_or) > 1e-5) || (abs(handles.x_max-handles.x_max_or) > 1e-5) || ...
            (abs(handles.y_min-handles.y_min_or) > 1e-5) || (abs(handles.y_max-handles.y_max_or) > 1e-5))
		opt_R = ['-R' sprintf('%.8g',handles.x_min) '/' sprintf('%.8g',handles.x_max) '/' ...
            sprintf('%.8g',handles.y_min) '/' sprintf('%.8g',handles.y_max)];
		out.x_min = handles.x_min;    out.x_max = handles.x_max;
		out.y_min = handles.y_min;    out.y_max = handles.y_max;
	end

	% See if grid increments were changed
	if ( (abs(handles.x_inc-handles.x_inc_or) > 1e-6) || (abs(handles.y_inc-handles.y_inc_or) > 1e-6) )
		opt_I = ['-I' sprintf('%.12g',handles.x_inc) '/' sprintf('%.12g',handles.y_inc)];
	end

	opt_F = [handles.command{1} handles.command{2}];
	opt_D = handles.command{4};

	set(handles.figure1,'Name','COMPUTING')
	head = handles.head;
	if ( strcmp(opt_R,' ') && strcmp(opt_I,' '))
		Z = grdfilter_m(handles.Z,head,opt_F,opt_D);
	elseif ( strcmp(opt_R,' ') && ~strcmp(opt_I,' '))
		[Z,head] = grdfilter_m(handles.Z,head,opt_I,opt_F,opt_D);
	elseif ( ~strcmp(opt_R,' ') && strcmp(opt_I,' '))
		[Z,head] = grdfilter_m(handles.Z,head,opt_R,opt_F,opt_D);
	else
		errordlg('Bad number of options (My fault also that didn''t guess all possible nonsenses).','Error');
		set(handles.figure1,'Name','Grdfilter')
		return
	end
	set(handles.figure1,'Name','Grdfilter')

	[ny,nx] = size(Z);
	X = linspace(head(1),head(2),nx);       Y = linspace(head(3),head(4),ny);
    zMinMax = grdutils(Z,'-L');
	tmp.head = [head(1:4) zMinMax(1:2)' head(7:9)];
	tmp.X = X;    tmp.Y = Y;    tmp.name = 'grdfiltered grid';
	new_window = mirone(Z,tmp);

% --------------------------------------------------------------------------------
function pushbutton_cancel_Callback(hObject, eventdata, handles)
	delete(handles.figure1);

% --------------------------------------------------------------------
function Menu_Help_Callback(hObject, eventdata, handles)
message = ['grdfilter will filter a .grd file in the time domain using ' ...
       'a boxcar, cosine arch, gaussian, median,  or  mode  filter ' ...
       'and  computing  distances  using  Cartesian  or  Spherical ' ...
       'geometries.  The output .grd file can optionally be gener- ' ...
       'ated  as  a  sub-Region  of  the  input  and/or with a new ' ...
       'Increment. In this way, one may have "extra space" in the ' ...
       'input data so that the edges will not be used and the out- ' ...
       'put can be within one-half-width of the input  edges.  If ' ...
       'the  filter  is low-pass, then the output may be less fre- ' ...
       'quently sampled than the input.'];
helpdlg(message,'Help on grdfilter');

% --------------------------------------------------------------------
function about_window_Callback(hObject, eventdata, handles)
About_box(['grdfilter_Last_modified_at_02_Nov_2004']);

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
	delete(handles.figure1);

% --- Executes on key press over figure1 with no controls selected.
function figure1_KeyPressFcn(hObject, eventdata, handles)
	% Check for "escape"
	if isequal(get(hObject,'CurrentKey'),'escape')
		delete(handles.figure1);
	end

% --- Creates and returns a handle to the GUI figure. 
function grdfilter_Mir_LayoutFcn(h1,handles);

set(h1,'PaperUnits',get(0,'defaultfigurePaperUnits'),...
'CloseRequestFcn',{@figure1_CloseRequestFcn,handles},...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'KeyPressFcn',{@figure1_KeyPressFcn,handles},...
'MenuBar','none',...
'Name','Grdfilter',...
'NumberTitle','off',...
'Position',[265.768111202607 370.94409531652 411 186],...
'RendererMode','manual',...
'Resize','off',...
'Tag','figure1');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[10 78 390 97],...
'Style','frame');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[30 168 121 15],...
'String','Griding Line Geometry',...
'Style','text',...
'Tag','text8');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@grdfilter_Mir_uicallback,h1,'edit_x_min_Callback'},...
'HorizontalAlignment','left',...
'Position',[77 136 80 21],...
'Style','edit',...
'Tag','edit_x_min');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@grdfilter_Mir_uicallback,h1,'edit_x_max_Callback'},...
'HorizontalAlignment','left',...
'Position',[163 136 80 21],...
'Style','edit',...
'Tag','edit_x_max');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@grdfilter_Mir_uicallback,h1,'edit_y_min_Callback'},...
'HorizontalAlignment','left',...
'Position',[77 110 80 21],...
'Style','edit',...
'Tag','edit_y_min');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@grdfilter_Mir_uicallback,h1,'edit_y_max_Callback'},...
'HorizontalAlignment','left',...
'Position',[163 110 80 21],...
'Style','edit',...
'Tag','edit_y_max');

uicontrol('Parent',h1,...
'Enable','inactive',...
'HorizontalAlignment','left',...
'Position',[22 140 55 15],...
'String','X Direction',...
'Style','text',...
'Tag','text2');

uicontrol('Parent',h1,...
'Enable','inactive',...
'HorizontalAlignment','left',...
'Position',[21 114 55 15],...
'String','Y Direction',...
'Style','text',...
'Tag','text3');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[183 157 41 13],...
'String','Max',...
'Style','text',...
'Tag','text4');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[99 157 41 13],...
'String','Min',...
'Style','text',...
'Tag','text5');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@grdfilter_Mir_uicallback,h1,'edit_x_inc_Callback'},...
'HorizontalAlignment','left',...
'Position',[248 136 71 21],...
'Style','edit',...
'TooltipString','DX grid spacing',...
'Tag','edit_x_inc');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@grdfilter_Mir_uicallback,h1,'edit_y_inc_Callback'},...
'HorizontalAlignment','left',...
'Position',[248 110 71 21],...
'Style','edit',...
'TooltipString','DY grid spacing',...
'Tag','edit_y_inc');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@grdfilter_Mir_uicallback,h1,'edit_Ncols_Callback'},...
'HorizontalAlignment','left',...
'Position',[324 136 65 21],...
'Style','edit',...
'TooltipString','Number of columns in the grid',...
'Tag','edit_Ncols');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@grdfilter_Mir_uicallback,h1,'edit_Nrows_Callback'},...
'HorizontalAlignment','left',...
'Position',[324 110 65 21],...
'Style','edit',...
'TooltipString','Number of rows in the grid',...
'Tag','edit_Nrows');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[265 159 41 13],...
'String','Spacing',...
'Style','text',...
'Tag','text6');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[332 159 51 13],...
'String','# of lines',...
'Style','text',...
'Tag','text7');

uicontrol('Parent',h1,...
'BackgroundColor',[0.831372559070587 0.815686285495758 0.7843137383461],...
'Callback',{@grdfilter_Mir_uicallback,h1,'pushbutton_Help_R_F_T_Callback'},...
'FontWeight','bold',...
'ForegroundColor',[0 0 1],...
'Position',[327 85 61 18],...
'String','?',...
'Tag','pushbutton_Help_R_F_T');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@grdfilter_Mir_uicallback,h1,'popup_FilterType_Callback'},...
'HorizontalAlignment','right',...
'Position',[20 36 122 22],...
'String',{'boxcar'; 'cosine arch'; 'gaussian'; 'median'; 'maximum likelihood'},...
'Style','popupmenu',...
'Value',1,...
'Tag','popup_FilterType');

uicontrol('Parent',h1,...
'BackgroundColor',[1 0.49 0.49],...
'Callback',{@grdfilter_Mir_uicallback,h1,'edit_FilterWidth_Callback'},...
'HorizontalAlignment','left',...
'Position',[142 37 47 21],...
'Style','edit',...
'TooltipString','specify filter full width',...
'Tag','edit_FilterWidth');

uicontrol('Parent',h1,...
'Callback',{@grdfilter_Mir_uicallback,h1,'pushbutton_Help_FilterType_Callback'},...
'FontWeight','bold',...
'ForegroundColor',[0 0 1],...
'Position',[190 36 21 23],...
'String','?',...
'Tag','pushbutton_Help_FilterType');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@grdfilter_Mir_uicallback,h1,'popup_Option_D_Callback'},...
'Position',[280 38 72 22],...
'String',{  '0'; '1'; '2'; '3'; '4' },...
'Style','popupmenu',...
'TooltipString','You better read the help availabe on box aside',...
'Value',1,...
'Tag','popup_Option_D');

uicontrol('Parent',h1,...
'Callback',{@grdfilter_Mir_uicallback,h1,'pushbutton_Help_Option_D_Callback'},...
'FontWeight','bold',...
'ForegroundColor',[0 0 1],...
'Position',[359 38 21 23],...
'String','?',...
'Tag','pushbutton_Help_Option_D');

uicontrol('Parent',h1,...
'Callback',{@grdfilter_Mir_uicallback,h1,'pushbutton_cancel_Callback'},...
'Position',[220 6 66 23],...
'String','Cancel',...
'Tag','pushbutton_cancel');

uimenu('Parent',h1,...
'Callback',{@grdfilter_Mir_uicallback,h1,'Menu_Help_Callback'},...
'Label','Help',...
'Tag','Menu_Help');

uimenu('Parent',h1,...
'Callback',{@grdfilter_Mir_uicallback,h1,'about_window_Callback'},...
'Label','About',...
'Tag','about_window');

uicontrol('Parent',h1,...
'Callback',{@grdfilter_Mir_uicallback,h1,'pushbutton_Compute_Callback'},...
'FontWeight','bold',...
'Position',[290 6 111 23],...
'String','Compute',...
'TooltipString','Write the command line into the script file',...
'Tag','pushbutton_Compute');

uicontrol('Parent',h1,...
'Enable','inactive',...
'HorizontalAlignment','left',...
'Position',[53 60 47 15],...
'String','Filter type',...
'Style','text',...
'Tag','text10');

uicontrol('Parent',h1,...
'Enable','inactive',...
'HorizontalAlignment','left',...
'Position',[144 59 53 15],...
'String','Filter width',...
'Style','text',...
'Tag','text11');

uicontrol('Parent',h1,...
'Enable','inactive',...
'HorizontalAlignment','left',...
'Position',[283 61 63 15],...
'String','Distance flag',...
'Style','text',...
'Tag','text12');

function grdfilter_Mir_uicallback(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
feval(callback_name,hObject,[],guidata(h1));

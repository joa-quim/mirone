function varargout = tsunamovie(varargin)
% Create .gif or .avi movies from a series of tsunami grids
% M-File changed by desGUIDE 

    hObject = figure('Tag','figure1','Visible','off');
    handles = guihandles(hObject);
    guidata(hObject, handles);
    tsunamovie_LayoutFcn(hObject,handles);
    handles = guihandles(hObject);
    movegui(hObject,'northeast')

	home_dir = cd;
	f_path = [home_dir filesep 'data' filesep];
	handles.path_tmp = [home_dir filesep 'tmp' filesep];
    handles.stop = 0;           % When the running engine detects it has canged to 1 it stops
    handles.dither = 'nodither';% Default
    handles.checkedMM = 0;      % To signal if need or not to get the ensemble water Min/Max
    handles.usrMM = 0;          % To signal if user has changed the ensemble Min|Max
    handles.lambCteComm = '/0.55/0.6/0.4/10';   % These Lamertian params are here const
    handles.waterIllumComm = '-E0/30/0.55/0.6/0.4/10';      % Starting values
    handles.landIllumComm = '-A0';
    handles.landCurrIllumType  = 'grdgradient_A';
    handles.waterCurrIllumType = 'lambertian';
    handles.fps = 5;            % Frames per second
    handles.dt = 0.2;           % 1/fps
    handles.scaleFactor = 1;    % To shrink or increase the movie dimensions
	handles.Z_bat   = [];
	handles.Z_water = [];
	handles.nameList = [];
    
    % Load some icons and put them in the toggles
    load([f_path 'mirone_icons.mat'],'um_ico','dois_ico','Mfopen_ico');
    set(handles.toggle_1,'CData',um_ico)
    set(handles.toggle_2,'CData',dois_ico)
    set(handles.push_batGrid,'CData',Mfopen_ico)
    set(handles.push_singleWater,'CData',Mfopen_ico)
    set(handles.push_namesList,'CData',Mfopen_ico)
    set(handles.push_movieName,'CData',Mfopen_ico)
    
    % Import background image
    astrolabio = imread([f_path 'astrolabio.jpg']);
    image(astrolabio,'parent',handles.axes1);

    pos = get(handles.axes1,'Position');
    set(handles.axes1,'Visible','off')

	% Draw everything that may be needed for all options. Later, depending on the
	% option selected, only the allowed features will be let visible
	x0 = pos(3)/2;      y0 = pos(4)/2;      radius = pos(3)/2 - 3;
	h_line(1) = line('parent',handles.axes1,'XData',[x0 x0],'YData',[y0 0],'Color','r','Tag','red','LineWidth',3,'Userdata',radius);
	% Now draw, on axes2, a quarter of circle and a line
	t = 0:0.02:pi/2;    x = [0 cos(t) 0];     y = [0 sin(t) 0];
	line('parent',handles.axes2,'XData',x,'YData',y,'HitTest','off','Color','k','LineWidth',1);
    xdataElev = [0 cos(30*pi/180)];     ydataElev = [0 sin(30*pi/180)];
	h_line(2) = line('parent',handles.axes2,'XData',xdataElev,'YData',ydataElev,'Color','k','LineWidth',3,'Visible','off');
	set(h_line(2),'Tag','Elev','Userdata',1)        % save radius of circumscribed circle

    % Backup the graphical info
    handles.landLineAzBack = [x0 x0 y0 0];
    handles.waterLineAzBack = [x0 x0 y0 0];
    handles.landLineElevBack = [xdataElev ydataElev];
    handles.waterLineElevBack = [xdataElev ydataElev];
    handles.landAzStrBack  = '0';
    handles.waterAzStrBack = '0';
    handles.landElevStrBack  = '30';
    handles.waterElevStrBack = '30';
    
    handles.ciclePar = [x0 y0 radius];
    handles.h_line = h_line;
    guidata(hObject, handles);
    show_needed(handles,'grdgradient_A')
    set(handles.toggle_1,'Value',1)         % Start it in a pressed state
    set(hObject,'WindowButtonDownFcn',{@ButtonDown,h_line,handles});

% Choose default command line output for tsunamovie
handles.output = hObject;
guidata(hObject, handles);
set(hObject,'Visible','on');
if (nargout),   varargout{1} = hObject;     end

% -----------------------------------------------------------------------------------------
function edit_batGrid_Callback(hObject, eventdata, handles)
    fname = get(hObject,'String');
    push_batGrid_Callback([], [], handles, fname)

% -----------------------------------------------------------------------------------------
function push_batGrid_Callback(hObject, eventdata, handles, opt)
    if (nargin == 3)        % Direct call
    	[FileName,PathName] = uigetfile({'*.grd;*.GRD', 'Grid files (*.grd,*.GRD)';'*.*',...
                'All Files (*.*)'},'Select GMT grid');
	    pause(0.01);
	    if isequal(FileName,0);     return;     end
    else        % File name on input
        [PathName,FNAME,EXT] = fileparts(opt);
        PathName = [PathName filesep];      % To be coherent with the 'if' branch
        FileName = [FNAME EXT];
    end
	fname = [PathName FileName];
	
	[handles,handles.X_bat,handles.Y_bat,handles.Z_bat,handles.head_bat] = read_gmt_type_grids(handles,fname);
	if (isempty(handles.X_bat)),    return;     end
	
	set(handles.edit_batGrid,'String',fname)
    [dump,FNAME,EXT] = fileparts(FileName);
    EXT = '.gif';
    if (~get(handles.radio_gif,'Value')),   EXT = '.avi';   end
    handles.moviePato = PathName;
    handles.movieName = FNAME;
	set(handles.edit_movieName,'String',[PathName handles.movieName EXT])
	guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------------------
function edit_singleWater_Callback(hObject, eventdata, handles)
    fname = get(hObject,'String');
    push_singleWater_Callback([], [], handles, fname)

% -----------------------------------------------------------------------------------------
function push_singleWater_Callback(hObject, eventdata, handles, opt)
    if (nargin == 3)        % Direct call
        [FileName,PathName] = uigetfile({'*.grd;*.GRD', 'Grid files (*.grd,*.GRD)';'*.*', ...
                'All Files (*.*)'},'Select GMT grid');
        pause(0.01);
        if isequal(FileName,0);     return;     end
    else        % File name on input
        [PathName,FNAME,EXT] = fileparts(opt);
        PathName = [PathName filesep];      % To be coherent with the 'if' branch
        FileName = [FNAME EXT];
    end
    fname = [PathName FileName];

	[handles,X,Y,handles.Z_water,handles.head_water] = read_gmt_type_grids(handles,fname);
	if (isempty(X)),    return;     end
	set(handles.push_clearTestBat,'Visible','on')
	
	set(handles.edit_singleWater,'String',fname)
	guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------------------
function edit_namesList_Callback(hObject, eventdata, handles)
    fname = get(hObject,'String');
    push_namesList_Callback([], [], handles, fname)
    
% -----------------------------------------------------------------------------------------
function push_namesList_Callback(hObject, eventdata, handles, opt)
    if (nargin == 3)        % Direct call
    	str1 = {'*.dat;*.DAT;*.txt;*.TXT', 'Data files (*.dat,*.DAT,*.txt,*.TXT)';'*.*', 'All Files (*.*)'};
        [FileName,PathName] = uigetfile(str1,'File with grids list');
	    if isequal(FileName,0);     return;     end
    else        % File name on input
        [PathName,FNAME,EXT] = fileparts(opt);
        PathName = [PathName filesep];      % To be coherent with the 'if' branch
        FileName = [FNAME EXT];
    end
	fid = fopen([PathName FileName]);
	c = char(fread(fid))';      fclose(fid);
	names = strread(c,'%s','delimiter','\n');   clear c fid;
	m = length(names);
    handles.shortNameList = cell(m,1);      % To hold grid names with path striped
	for (k=1:m)
        [PATH,FNAME,EXT] = fileparts(names{k});
        if (isempty(PATH))
            handles.shortNameList{k} = names{k};
            names{k} = [PathName names{k}];
        else
            handles.shortNameList{k} = [FNAME EXT];
        end
	end

    handles.nameList = names;
    set(handles.listbox1,'String',handles.shortNameList)
    set(handles.edit_namesList,'String',[PathName FileName])
    set(handles.push_checkGlobalMM,'Visible','on')
    guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------------------
function edit_movieName_Callback(hObject, eventdata, handles)
    fname = get(hObject,'String');
    push_movieName_Callback([], [], handles, fname)

% -----------------------------------------------------------------------------------------
function push_movieName_Callback(hObject, eventdata, handles, opt)
    if (nargin == 3)        % Direct call
        [FileName,PathName] = uiputfile({'*.gif;*.avi', 'Grid files (*.gif,*.avi)'},'Select Movie name');
        pause(0.01);
        if isequal(FileName,0);     return;     end
        [dumb,FNAME,EXT]= fileparts(FileName);
    else        % File name on input
        [PathName,FNAME,EXT] = fileparts(opt);
        PathName = [PathName filesep];      % To be coherent with the 'if' branch
        FileName = [FNAME EXT];
    end
    if (~strmatch(lower(EXT),{'.gif' '.avi'}))
        errordlg('Ghrrrrrrrr! Don''t be smart. Only ''.gif'' or ''.avi'' extensions are acepted.','Chico Clever');
        return
    end
    
    handles.moviePato = PathName;
    handles.movieName = FileName;
    if (strcmpi(EXT,'.gif'))
        set(handles.radio_gif,'Value',1)
        radio_gif_Callback(handles.radio_gif, [], handles)
    else
        set(handles.radio_avi,'Value',1)
        radio_avi_Callback(handles.radio_avi, [], handles)
    end
    guidata(handles.figure1,handles)
    
% -----------------------------------------------------------------------------------------
function listbox1_Callback(hObject, eventdata, handles)
    % if this is a doubleclick,
    if ( strcmp(get(gcbf,'SelectionType'),'open') && ~isempty(handles.nameList) )
        val = get(hObject,'Value');
        push_singleWater_Callback([], [], handles, handles.nameList{val})
    end

% -----------------------------------------------------------------------------------------
function toggle_1_Callback(hObject, eventdata, handles)
    show_needed(handles,'grdgradient_A')
    if (get(handles.radio_land,'Value')),   handles.landCurrIllumType = 'grdgradient_A';
    else                                    handles.waterCurrIllumType = 'grdgradient_A';
    end
    guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------------------
function toggle_2_Callback(hObject, eventdata, handles)
    show_needed(handles,'lambertian')
    if (get(handles.radio_land,'Value')),   handles.landCurrIllumType = 'lambertian';
    else                                    handles.waterCurrIllumType = 'lambertian';
    end
    guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------------------
function radio_land_Callback(hObject, eventdata, handles)
    if (get(hObject,'Value'))
        set(handles.radio_water,'Value',0)
        set(handles.h_line(1),'XData',handles.landLineAzBack(1:2),'YData',handles.landLineAzBack(3:4))
        set(handles.h_line(2),'XData',handles.landLineElevBack(1:2),'YData',handles.landLineElevBack(3:4))
        show_needed(handles,handles.landCurrIllumType)
        set(handles.edit_azim,'String', handles.landAzStrBack)
        set(handles.edit_elev,'String', handles.landElevStrBack)
        str = sprintf(['Set parametrs for Land illumination\nCurrent selection (in GMT parlance) is:\n',...
            handles.landIllumComm]);
        set(hObject,'TooltipString',str)
    else
        set(hObject,'Value',1)
    end

% -----------------------------------------------------------------------------------------
function radio_water_Callback(hObject, eventdata, handles)
    if (get(hObject,'Value'))
        set(handles.radio_land,'Value',0)
        set(handles.h_line(1),'XData',handles.waterLineAzBack(1:2),'YData',handles.waterLineAzBack(3:4))
        set(handles.h_line(2),'XData',handles.waterLineElevBack(1:2),'YData',handles.waterLineElevBack(3:4))
        show_needed(handles,handles.waterCurrIllumType)
        set(handles.edit_azim,'String', handles.waterAzStrBack)
        set(handles.edit_elev,'String', handles.waterElevStrBack)
        str = sprintf(['Set parametrs for Water illumination\nCurrent selection (in GMT parlance) is:\n',...
            handles.waterIllumComm]);
        set(hObject,'TooltipString',str)
    else
        set(hObject,'Value',1)
    end

% -----------------------------------------------------------------------------------------
function show_needed(handles,opt)
	h_all = handles.h_line;
	if (strncmp(opt,'grdgradient',11))
        set(handles.edit_elev,'Enable','off');          set(handles.edit_azim,'Visible','on')
        set(handles.text_elev,'Enable','on');           set(handles.edit_azim,'Enable','on');
        set(handles.text_azim,'Enable','on');
        set(h_all(1),'Visible','on');                   set(h_all(2),'Visible','off')
        set(handles.toggle_1,'Value',1);                set(handles.toggle_2,'Value',0)
	else                            % Lambertian
        set(handles.edit_elev,'Enable','on');           set(handles.edit_azim,'Visible','on')
        set(handles.text_elev,'Enable','on');           set(handles.edit_azim,'Enable','on');
        set(handles.text_azim,'Enable','on');           set(h_all,'Visible','on');
        set(handles.toggle_1,'Value',0);                set(handles.toggle_2,'Value',1)
	end
	guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------------------
function ButtonDown(obj,eventdata,h_all,handles)
	% It could be cleverer.
	pt = get(gca, 'CurrentPoint');
	x_lim = get(gca,'xlim');      y_lim = get(gca,'ylim');
	% check if x,y is inside of axis
	if ~((pt(1,1)>=x_lim(1)) & (pt(1,1)<=x_lim(2)) & (pt(1,2)>=y_lim(1)) & (pt(1,2)<=y_lim(2)))    % outside axis limits
        return
	end
	if any(h_all == gco)
        h = h_all(h_all == gco);    % When more than one line handle exists, find only the selected one
        set(gcf,'WindowButtonMotionFcn',{@ButtonMotion,h,handles},'WindowButtonUpFcn',{@ButtonUp,h_all,handles},...
            'Pointer', 'crosshair');
	else
        return
	end

% -----------------------------------------------------------------------------------------
function ButtonMotion(obj,eventdata,h,handles)
	selectionType = get(gcf, 'SelectionType');
	pt = get(gca, 'CurrentPoint');
	if strcmp(selectionType, 'normal')      % right-cick
        xx = get(h,'XData');    yy = get(h,'YData');
        theta = cart2pol(pt(1,1)-xx(1),pt(1,2)-yy(1));
        radius = get(h,'Userdata');
        x2 = xx(1) + radius * cos(theta);      y2 = yy(1) + radius * sin(theta);
        if strcmp(get(h,'Tag'),'Elev') & (theta >= 0 & theta <= pi/2)   % Elevation line
            set(h,'XData',[xx(1) x2],'YData',[yy(1) y2]);
            set(handles.edit_elev,'String',num2str(fix(theta *180/pi)) )
        elseif ~strcmp(get(h,'Tag'),'Elev')     % Azimuth line(s)
            set(h,'XData',[xx(1) x2],'YData',[yy(1) y2]);
            % NOTE to if I ever want to reuse this code. Normally ang_2pi should be = pi/2 - (pi*.....)
            % for the normal y origin at bottm left corner. However, due to the stupid habit of using y=0
            % at top left corner when dealing with images, to get an azimuth angle we have to do like following. 
	
            % truncate angles into [-pi pi] range
            ang_2pi = pi/2 + ( pi*((abs(theta)/pi) - 2*ceil(((abs(theta)/pi)-1)/2)) * sign(theta) );
            epsilon = -1e-7;        %  Allow points near zero to remain there
            indx = find(ang_2pi < epsilon);
            %  Shift the points in the [-pi 0] range to [pi 2pi] range
            if ~isempty(indx);  ang_2pi(indx) = ang_2pi(indx) + 2*pi;  end;
            set(handles.edit_azim,'String',sprintf('%.0f',ang_2pi *180/pi) )
        end
	end

% -----------------------------------------------------------------------------------------
function ButtonUp(obj,eventdata,h,handles)
    handles = guidata(handles.figure1);     % We need an updated version
    set(handles.figure1,'WindowButtonMotionFcn','','WindowButtonDownFcn', ...
        {@ButtonDown,h,handles},'WindowButtonUpFcn','');
    set(handles.figure1,'Pointer', 'arrow')
    azim = get(handles.edit_azim,'String');
    xdataAz = get(h(1),'XData');        ydataAz = get(h(1),'YData');
    xdataElev = get(h(2),'XData');      ydataElev = get(h(2),'YData');
    if (get(handles.radio_land,'Value'))    % We are setting the LAND illum params
        if (get(handles.toggle_1,'Value'))  % grdgradient classic
            handles.landIllumComm = ['-A' azim];
        else                                % Lambertian
            elev = get(handles.edit_elev,'String');
            handles.landIllumComm = ['-E' azim '/' elev handles.lambCteComm];
        end
        handles.landLineAzBack = [xdataAz ydataAz];         % Backup the graphical info
        handles.landLineElevBack = [xdataElev ydataElev];
        handles.landAzStrBack = get(handles.edit_azim,'String'); % And string-numeric
        handles.landElevStrBack = get(handles.edit_elev,'String');
    else                                    % WATER
        if (get(handles.toggle_1,'Value'))  % grdgradient classic
            handles.waterIllumComm = ['-A' azim];
        else                                % Lambertian
            elev = get(handles.edit_elev,'String');
            handles.waterIllumComm = ['-E' azim '/' elev handles.lambCteComm];
        end
        handles.waterLineAzBack = [xdataAz ydataAz];         % Backup the graphical info
        handles.waterLineElevBack = [xdataElev ydataElev];
        handles.waterAzStrBack = get(handles.edit_azim,'String'); % And string-numeric
        handles.waterElevStrBack = get(handles.edit_elev,'String');
    end
    guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------------------
function edit_elev_Callback(hObject, eventdata, handles)
    % No need
    
% -----------------------------------------------------------------------------------------
function edit_azim_Callback(hObject, eventdata, handles)
    % No need

% -----------------------------------------------------------------------------------------
function radio_gif_Callback(hObject, eventdata, handles)
    if (get(hObject,'Value')),      set(handles.radio_avi,'Value',0)
    else                            set(hObject,'Value',1)
    end
    mname = get(handles.edit_movieName,'String');
    if (~isempty(mname))
        mname = [handles.moviePato handles.movieName '.gif'];
        set(handles.edit_movieName,'String',mname)
    end

% -----------------------------------------------------------------------------------------
function radio_avi_Callback(hObject, eventdata, handles)
    if (get(hObject,'Value')),      set(handles.radio_gif,'Value',0)
    else                            set(hObject,'Value',1)
    end
    mname = get(handles.edit_movieName,'String');
    if (~isempty(mname))
        mname = [handles.moviePato handles.movieName '.avi'];
        set(handles.edit_movieName,'String',mname)
    end

% -----------------------------------------------------------------------------------------
function checkbox_dither_Callback(hObject, eventdata, handles)
    if (get(hObject,'Value')),      handles.dither = 'dither';
    else                            handles.dither = 'nodither';
    end
    guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------------------
function push_stop_Callback(hObject, eventdata, handles)
    handles.stop = 1;
    guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------------------
function edit_fps_Callback(hObject, eventdata, handles)
    % Frames per second
    fps = round(str2double(get(hObject,'String')));
    if (isnan(fps))
        set(hObject,'String',num2str(handles.fps))
        return
    end
    set(hObject,'String',num2str(fps))      % In case there were decimals
    handles.fps = fps;
    handles.dt = 1/fps;
    guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------------------
function push_clearTestBat_Callback(hObject, eventdata, handles)
    set(handles.edit_singleWater,'String','')
    set(hObject,'Visible','off')
    handles.Z_water = [];
    guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------------------
function slider_transparency_Callback(hObject, eventdata, handles)
    val = get(hObject,'Value');
    set(handles.text_WaterTrans,'String',['Water transparency  ' sprintf('%.0f',val*100) '%'])

% -----------------------------------------------------------------------------------------
function edit_globalWaterMax_Callback(hObject, eventdata, handles)
    % User can change the ensemble limits, but we must keep track of that
    xx = str2double(get(hObject,'String'));
    if (~isnan(xx))
        handles.maxWater = xx;
        handles.minWater = str2double(get(handles.edit_globalWaterMin,'String'));
        handles.usrMM = 1;          % Flag that user has changed the ensemble Min|Max
        guidata(handles.figure1,handles)
    end

% -----------------------------------------------------------------------------------------
function edit_globalWaterMin_Callback(hObject, eventdata, handles)
    % User can change the ensemble limits, but we must keep track of that
    xx = str2double(get(hObject,'String'));
    if (~isnan(xx))
        handles.minWater = xx;
        handles.maxWater = str2double(get(handles.edit_globalWaterMax,'String'));
        handles.usrMM = 1;          % Flag that user has changed the ensemble Min|Max
        guidata(handles.figure1,handles)
    end

% -----------------------------------------------------------------------------------------
function [minWater, maxWater, heads] = push_checkGlobalMM_Callback(hObject, eventdata, handles)
    [minWater, maxWater, heads] = get_globalMinMax(handles);
    set(handles.edit_globalWaterMin,'String',sprintf('%.2f',minWater))
    set(handles.edit_globalWaterMax,'String',sprintf('%.2f',maxWater))
    set(hObject,'Visible','off')
    handles.checkedMM = 1;          % Signal that ensemble Min/Max is now known
    handles.minWater = minWater;    % Make a copy that can be changed by a user
    handles.maxWater = maxWater;    % direct intervenction on the edit boxes
    handles.gridHeaders = heads;    % a cell array
    guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------------------
function [minWater, maxWater, heads] = get_globalMinMax(handles)
    % Run trough all water level grids and find the ensemble Min/Max
    minWater = 1000;        maxWater = -1000;
    nGrids = numel(handles.nameList);
    heads = cell(nGrids,1);
    for (i=1:nGrids)
        [handles,heads{i}] = read_gmt_type_grids(handles,handles.nameList{i},'hdr');
        minWater = min(minWater,heads{i}(5));
        maxWater = max(maxWater,heads{i}(6));
    end

% -----------------------------------------------------------------------------------------
function pushbutton_OK_Callback(hObject, eventdata, handles)

    if (isempty(handles.Z_bat))
        errordlg('Noooo! Where is the bathymetry file? Do you think I''m bruxo?','ERROR')
        return
    end
    if (isempty(handles.movieName) && isempty(handles.Z_water))
        errordlg('Hei! what shoult it be the movie name?','ERROR');        return
    end
    
    is_gif = 1;     % Veracity of this is checked bellow
    alfa = get(handles.slider_transparency,'Value');
    ind = (handles.Z_water == 0);
    tmp.head = [handles.head_bat(1:4) 0 255 handles.head_bat(7:9)];
    tmp.X = handles.X_bat;      tmp.Y = handles.Y_bat;
    tmp.geog = 1;               tmp.name = 'Tsu transparente';
	
	% Do ground illum
	imgBat = scaleto8(handles.Z_bat);
    imgBat = ind2rgb8(imgBat,jet(256));
	R = grdgradient_m(handles.Z_bat,handles.head_bat,'-M',handles.landIllumComm,'-Nt');
	imgBat = shading_mat(imgBat,R,'no_scale');

    % Initialize the circular waitbar
    x0 = handles.ciclePar(1);
    y0 = handles.ciclePar(2);
    radius = handles.ciclePar(3);
    x = [x0 x0 x0];     y = [y0 0 y0];
    hWC = patch('parent',handles.axes1,'XData',x,'YData',y,'FaceColor','b', 'EdgeColor','b');

    % Do the water illum
    if (~isempty(handles.Z_water))
    	imgWater = scaleto8(handles.Z_water);
        imgWater = ind2rgb8(imgWater,jet(256));
    	R = grdgradient_m(handles.Z_water,handles.head_water,handles.waterIllumComm);
    	imgWater = shading_mat(imgWater,R,'no_scale');
    	clear R;
        imgWater = mixe_images(handles, imgBat, imgWater, ind, alfa);
        %imgWater = cvlib_mex('text',imgWater,'1:32:45',[10 30]);
        mirone(imgWater,tmp);
    elseif (~isempty(handles.nameList))     % If we have a list of names
        nGrids = numel(handles.nameList);
        if (~handles.checkedMM)         % We don't know yet the water ensemble Min/Max
            [minWater, maxWater,heads] = push_checkGlobalMM_Callback([], [], handles);
        else                            % Get what we know
            minWater = handles.minWater;
            maxWater = handles.maxWater;
            heads = handles.gridHeaders;
        end
        set(handles.push_stop,'Visible','on')       % To allow interruptions
        if (~get(handles.radio_gif,'Value')),    is_gif = 0;     end
        str_nGrids = sprintf('%d',nGrids);
        for (i=1:nGrids)
            % Check if meanwhile the stop button has been pressed 
            pause(0.01)     % Little pause so that it can listen a eventual STOP request
            handles = guidata(handles.figure1); % We need an updated version
            if (handles.stop)       % Yes. Stop and reset things
                handles.stop = 0;   % Reset
                set(handles.figure1,'Name','Tsunamovie')
                set(handles.push_stop,'Visible','off')
                guidata(handles.figure1,handles)
                break
            end
            
            % Show visualy the processing advance
            set(handles.figure1,'Name',['Processing grid ' sprintf('%d',i) ' of ' str_nGrids]);
            theta = 2*pi * i / nGrids;
            t = (0:0.02:theta) -pi/2;
            set(hWC,'XData',[x0 x0+radius*cos(t) x0],'YData',[y0 y0+radius*sin(t) y0])

            [handles,X,Y,Z] = read_gmt_type_grids(handles,handles.nameList{i});
        	imgWater = scaleto8(Z,8,[minWater maxWater]);
            imgWater = ind2rgb8(imgWater,jet(256));
        	R = grdgradient_m(Z,heads{i},handles.waterIllumComm);
        	imgWater = shading_mat(imgWater,R,'no_scale');
            ind = (Z == 0);
%             imgWater = mixe_images(handles, imgBat, imgWater, ind, alfa);
%             imgWater = flipdim(imgWater,1);     % The stupid UL origin
%             imgWater = cvlib_mex('text',imgWater,'1:32:45',[10 30]);
            if (is_gif)
                [imgWater,map] = img_fun('rgb2ind',imgWater,256,handles.dither);
            end
            imgWater = flipdim(imgWater,1);     % The stupid UL origin
            if (is_gif)
                mname = [handles.moviePato handles.movieName '.gif'];
                if (i == 1)
                    writegif(imgWater,map,mname,'loopcount',Inf)
                else
                    writegif(imgWater,map,mname,'WriteMode','append','DelayTime',handles.dt)
                end
            else        % AVI
                M(i) = im2frame(imgWater);
            end
        end
        if (~is_gif)    % AVI
            mname = [handles.moviePato handles.movieName '.avi'];
      	    movie2avi_j(M,mname,'compression','none','fps',handles.fps)
        end
        set(handles.figure1,'Name','Tsunamovie')
        set(handles.push_stop,'Visible','off')
    else
        errordlg('Where are the actors? There is no water to make the WaterWorld movie!','ERROR')
    end
    try,    delete(hWC);    end

% -----------------------------------------------------------------------------------------
function imgWater = mixe_images(handles, imgBat, imgWater, ind, alfa)
    % Mixes land and water images simulating transparency.
    % It also resizes the image if the scale factor is ~= 1
    
    try
        if (alfa > 0.01)    % Only if transparency is greater than 1%
            cvlib_mex('addweighted',imgWater,(1 - alfa),imgBat,alfa);     % In-place
        end
        
        tmpW = imgWater(:,:,1);     tmpB = imgBat(:,:,1);           % R
        tmpW(ind) = tmpB(ind);      imgWater(:,:,1) = tmpW;
        
        tmpW = imgWater(:,:,2);     tmpB = imgBat(:,:,2);           % G
        tmpW(ind) = tmpB(ind);      imgWater(:,:,2) = tmpW;
        
        tmpW = imgWater(:,:,3);     tmpB = imgBat(:,:,3);           % B
        tmpW(ind) = tmpB(ind);      imgWater(:,:,3) = tmpW;
        
        if (handles.scaleFactor ~= 1)
            imgWater = cvlib_mex('resize',imgWater,handles.scaleFactor);
        end
    catch
        errordlg(lasterr,'Error')
    end
    
    %imgWater = cvlib_mex('text',imgWater,'1:32:45',[10 30]);    

% -----------------------------------------------------------------------------------------
function popup_resize_Callback(hObject, eventdata, handles)
    contents = get(hObject,'String');
    handles.scaleFactor = str2double(contents{get(hObject,'Value')});
    guidata(handles.figure1,handles)


% --- Creates and returns a handle to the GUI figure. 
function tsunamovie_LayoutFcn(h1,handles);

set(h1,...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'DoubleBuffer','on',...
'MenuBar','none',...
'Name','TsunaMovie',...
'NumberTitle','off',...
'PaperSize',[20.98404194812 29.67743169791],...
'Position',[520 281 451 409],...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1');

uicontrol('Parent',h1,'Position',[370 131 71 61],'Style','frame','Tag','frame5');
uicontrol('Parent',h1,'Position',[10 8 221 41],'Style','frame','Tag','frame3');
uicontrol('Parent',h1,'Position',[260 131 101 61],'Style','frame','Tag','frame2');
uicontrol('Parent',h1,'Position',[260 57 181 61],'Style','frame','Tag','frame1');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@tsunamovie_uicallback,h1,'edit_batGrid_Callback'},...
'HorizontalAlignment','left',...
'Position',[10 369 200 21],...
'Style','edit',...
'TooltipString','Name of bathymetry grid',...
'Tag','edit_batGrid');

uicontrol('Parent',h1,...
'Callback',{@tsunamovie_uicallback,h1,'push_batGrid_Callback'},...
'Position',[210 369 21 21],...
'TooltipString','Browse for a bathymetry grid',...
'Tag','push_batGrid');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@tsunamovie_uicallback,h1,'edit_singleWater_Callback'},...
'HorizontalAlignment','left',...
'Position',[20 13 180 21],...
'Style','edit',...
'TooltipString','Name of one water level grid for test illumination purposes',...
'Tag','edit_singleWater');

uicontrol('Parent',h1,...
'Callback',{@tsunamovie_uicallback,h1,'push_singleWater_Callback'},...
'Position',[200 13 21 21],...
'Tag','push_singleWater');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@tsunamovie_uicallback,h1,'edit_namesList_Callback'},...
'HorizontalAlignment','left',...
'Position',[10 319 200 21],...
'Style','edit',...
'TooltipString','Name of a file with the water level grids list',...
'Tag','edit_namesList');

uicontrol('Parent',h1,...
'Callback',{@tsunamovie_uicallback,h1,'push_namesList_Callback'},...
'Position',[210 319 21 21],...
'TooltipString','Browse for a water level grids list file',...
'Tag','push_namesList');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@tsunamovie_uicallback,h1,'listbox1_Callback'},...
'Position',[10 127 221 191],...
'Style','listbox',...
'Value',1,...
'Tag','listbox1');

axes('Parent',h1,...
'Units','pixels',...
'CameraPosition',[0.5 0.5 9.16025403784439],...
'Position',[260 287 91 91],...
'Tag','axes1');

axes('Parent',h1,...
'Units','pixels',...
'CameraPosition',[0.5 0.5 9.16025403784439],...
'Position',[370 288 51 51],...
'Tag','axes2',...
'Visible','off');

uicontrol('Parent',h1,...
'Callback',{@tsunamovie_uicallback,h1,'toggle_1_Callback'},...
'Position',[260 388 20 20],...
'Style','togglebutton',...
'TooltipString','GMT grdgradient classic Illumination',...
'Tag','toggle_1');

uicontrol('Parent',h1,...
'Callback',{@tsunamovie_uicallback,h1,'toggle_2_Callback'},...
'Position',[281 388 20 20],...
'Style','togglebutton',...
'TooltipString','Lambertian with lighting Illumination',...
'Tag','toggle_2');


uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@tsunamovie_uicallback,h1,'edit_elev_Callback'},...
'Position',[380 264 30 18],...
'String','30',...
'Style','edit',...
'TooltipString','Elevation light direction',...
'Tag','edit_elev');

uicontrol('Parent',h1,...
'Callback',{@tsunamovie_uicallback,h1,'pushbutton_OK_Callback'},...
'FontName','Helvetica',...
'FontSize',10,...
'Position',[375 9 66 23],...
'String','OK',...
'Tag','pushbutton_OK');

uicontrol('Parent',h1,...
'FontName','Elvetica',...
'FontSize',9,...
'Position',[367 343 50 16],...
'String','Elevation',...
'Style','text',...
'Tag','text_elev');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@tsunamovie_uicallback,h1,'edit_azim_Callback'},...
'Position',[305 264 34 18],...
'String','0',...
'Style','edit',...
'TooltipString','Azimuth direction',...
'Tag','edit_azim');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'FontSize',9,...
'HorizontalAlignment','left',...
'Position',[261 264 42 16],...
'String','Azimuth',...
'Style','text',...
'Tag','text_azim');

uicontrol('Parent',h1,...
'Callback',{@tsunamovie_uicallback,h1,'radio_avi_Callback'},...
'FontName','Helvetica',...
'Position',[271 65 41 15],...
'String','AVI',...
'Style','radiobutton',...
'TooltipString','Write movie file in RGB AVI format',...
'Tag','radio_avi');

uicontrol('Parent',h1,...
'Callback',{@tsunamovie_uicallback,h1,'radio_land_Callback'},...
'FontName','Helvetica',...
'FontSize',9,...
'Position',[317 388 50 17],...
'String','Land',...
'Style','radiobutton',...
'TooltipString','Set parametrs for Land illumination',...
'Value',1,...
'Tag','radio_land');

uicontrol('Parent',h1,...
'Callback',{@tsunamovie_uicallback,h1,'radio_water_Callback'},...
'FontName','Helvetica',...
'FontSize',9,...
'Position',[379 388 55 17],...
'String','Water',...
'Style','radiobutton',...
'TooltipString','Set parametrs for Water illumination',...
'Tag','radio_water');

uicontrol('Parent',h1,...
'Callback',{@tsunamovie_uicallback,h1,'radio_gif_Callback'},...
'FontName','Helvetica',...
'Position',[271 86 41 15],...
'String','GIF',...
'Style','radiobutton',...
'TooltipString','Write movie file in animated GIF format',...
'Value',1,...
'Tag','radio_gif');

uicontrol('Parent',h1,...
'Callback',{@tsunamovie_uicallback,h1,'checkbox_dither_Callback'},...
'FontName','Helvetica',...
'Position',[390 89 48 15],...
'String','Dither',...
'Style','checkbox',...
'TooltipString','If you don''t know what this is, ask google',...
'Tag','checkbox_dither');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'FontSize',9,...
'Position',[272 110 70 15],...
'String','Movie type',...
'Style','text',...
'Tag','text_MovType');

uicontrol('Parent',h1,...
'Callback',{@tsunamovie_uicallback,h1,'push_stop_Callback'},...
'FontName','Helvetica',...
'FontSize',10,...
'ForegroundColor',[1 0 0],...
'Position',[260 9 66 23],...
'String','STOP',...
'Tag','push_stop',...
'Visible','off');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@tsunamovie_uicallback,h1,'edit_fps_Callback'},...
'Position',[390 66 30 18],...
'String','5',...
'Style','edit',...
'TooltipString','Frames per second',...
'Tag','edit_fps');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'HorizontalAlignment','left',...
'Position',[333 68 55 15],...
'String','Frames p/s',...
'Style','text',...
'Tag','text9');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'FontSize',9,...
'HorizontalAlignment','left',...
'Position',[10 392 82 17],...
'String','Bathymetry file',...
'Style','text',...
'Tag','text10');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'FontSize',9,...
'HorizontalAlignment','left',...
'Position',[11 341 82 17],...
'String','Water files list',...
'Style','text',...
'Tag','text11');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'FontSize',9,...
'Position',[20 40 100 17],...
'String','Test with this file',...
'Style','text',...
'Tag','text_testFile');

uicontrol('Parent',h1,...
'Callback',{@tsunamovie_uicallback,h1,'push_clearTestBat_Callback'},...
'FontName','Helvetica',...
'FontSize',7,...
'Position',[181 33 40 16],...
'String','Clear',...
'TooltipString','Remove the test grid so that you can compute a movie',...
'Tag','push_clearTestBat',...
'Visible','off');

 uicontrol('Parent',h1,...
'BackgroundColor',[0.9 0.9 0.9],...
'Callback',{@tsunamovie_uicallback,h1,'slider_transparency_Callback'},...
'Position',[260 217 180 16],...
'Style','slider',...
'Value',0.2,...
'Tag','slider_transparency');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'FontSize',9,...
'HorizontalAlignment','left',...
'Position',[260 234 150 17],...
'String','Water transparency  20%',...
'Style','text',...
'Tag','text_WaterTrans');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@tsunamovie_uicallback,h1,'edit_globalWaterMax_Callback'},...
'Position',[300 159 50 19],...
'Style','edit',...
'TooltipString','Global maximum water level',...
'Tag','edit_globalWaterMax');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'FontSize',9,...
'HorizontalAlignment','left',...
'Position',[272 160 25 16],...
'String','Max',...
'Style','text',...
'Tag','text14');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@tsunamovie_uicallback,h1,'edit_globalWaterMin_Callback'},...
'Position',[300 138 50 19],...
'Style','edit',...
'TooltipString','Global minimum water level',...
'Tag','edit_globalWaterMin');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'FontSize',9,...
'HorizontalAlignment','left',...
'Position',[272 140 25 16],...
'String','Min',...
'Style','text',...
'Tag','text15');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'FontSize',9,...
'Position',[265 184 90 17],...
'String','Global min/max',...
'Style','text',...
'Tag','text_globalMM');

uicontrol('Parent',h1,...
'Callback',{@tsunamovie_uicallback,h1,'push_checkGlobalMM_Callback'},...
'FontName','Helvetica',...
'Position',[111 341 120 18],...
'String','Check global Min/Max',...
'TooltipString','Run trough all grids to find ensemble Min/Max',...
'Tag','push_checkGlobalMM',...
'Visible','off');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@tsunamovie_uicallback,h1,'popup_resize_Callback'},...
'Position',[380 147 50 22],...
'String',{'0.25'; '0.33'; '0.5'; '1.0'; '2.0' },...
'Style','popupmenu',...
'TooltipString','Resize output images by this value',...
'Value',4,...
'Tag','popup_resize');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@tsunamovie_uicallback,h1,'edit_movieName_Callback'},...
'HorizontalAlignment','left',...
'Position',[10 80 200 21],...
'Style','edit',...
'TooltipString','Name of movie file',...
'Tag','edit_movieName');

uicontrol('Parent',h1,...
'Callback',{@tsunamovie_uicallback,h1,'push_movieName_Callback'},...
'Position',[210 80 21 21],...
'TooltipString','Browse for a movie file name (extention is ignored)',...
'Tag','push_movieName');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'FontSize',9,...
'FontWeight','bold',...
'HorizontalAlignment','left',...
'Position',[10 103 115 17],...
'String','Output movie name',...
'Style','text',...
'Tag','text17');

uicontrol('Parent',h1,'ForegroundColor',[0 0.5 0],'Position',[240 8 2 401],'Style','frame','Tag','frame4');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'FontSize',9,...
'Position',[377 183 50 17],...
'String','Scale it?',...
'Style','text',...
'Tag','text18');

% ---------------------------------------------------------------------------------
function tsunamovie_uicallback(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
feval(callback_name,hObject,[],guidata(h1));

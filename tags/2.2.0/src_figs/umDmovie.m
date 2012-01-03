function varargout = umDmovie(varargin)
% ...

%	Copyright (c) 2004-2012 by J. Luis
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

    hObject = figure('Tag','figure1','Visible','off');
    umDmovie_LayoutFcn(hObject);
    handles = guihandles(hObject);
    move2side(hObject,'right')
 
    if (numel(varargin) > 0 && isstruct(varargin{1}))
        handMir = varargin{1};
        handles.work_dir = handMir.work_dir;
        handles.last_dir = handMir.last_dir;
    	handles.home_dir = handMir.home_dir;
    else
    	handles.home_dir = cd;
        handles.last_dir = handles.home_dir;
        handles.work_dir = handles.home_dir;
    end
    
	f_path = [handles.home_dir filesep 'data' filesep];
    handles.stop = 0;           % When the running engine detects it has canged to 1 it stops
    handles.fps = 5;            % Frames per second
    handles.dt = 0.2;           % 1/fps
    handles.imgW = 600;         % Default figure size (will also be image size)
    handles.imgH = 300;
	handles.profile  = [];
	handles.Z_bat    = [];
	handles.Z_water  = [];
	handles.nameList = [];
    handles.testTime = [];
    
    % Load some icons and put them in the toggles
    load([f_path 'mirone_icons.mat'],'Mfopen_ico');
    set(handles.push_batGrid,'CData',Mfopen_ico)
    set(handles.push_singleWater,'CData',Mfopen_ico)
    set(handles.push_namesList,'CData',Mfopen_ico)
    set(handles.push_movieName,'CData',Mfopen_ico)
    set(handles.push_profile,'CData',Mfopen_ico)

    guidata(hObject, handles);
    set(hObject,'Visible','on');
	if (nargout),	varargout{1} = hObject;		end

% -----------------------------------------------------------------------------------------
function edit_batGrid_CB(hObject, handles)
    fname = get(hObject,'String');
    push_batGrid_CB([], [], handles, fname)

% -----------------------------------------------------------------------------------------
function push_batGrid_CB(hObject, handles, opt)
    if (nargin == 3)        % Direct call
        cd(handles.last_dir)
    	[FileName,PathName] = uigetfile({'*.grd;*.GRD', 'Grid files (*.grd,*.GRD)';'*.*',...
                'All Files (*.*)'},'Select GMT grid');
	    pause(0.01);        cd(handles.home_dir);
	    if isequal(FileName,0);     return;     end
        if (PathName ~= 0),         handles.last_dir = PathName;    end
    else        % File name on input
        [PathName,FNAME,EXT] = fileparts(opt);
        PathName = [PathName filesep];      % To be coherent with the 'if' branch
        FileName = [FNAME EXT];
    end
	fname = [PathName FileName];
	
	[handles,handles.X_bat,handles.Y_bat,handles.Z_bat,handles.head_bat] = read_gmt_type_grids(handles,fname);
	if (isempty(handles.X_bat)),    return;     end
	
	set(handles.edit_batGrid,'String',fname)
    [dump,FNAME] = fileparts(FileName);
    EXT = '.gif';
    if (~get(handles.radio_gif,'Value')),   EXT = '.avi';   end
    handles.moviePato = PathName;
    handles.movieName = FNAME;
	set(handles.edit_movieName,'String',[PathName handles.movieName EXT])
	guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------------------
function edit_singleWater_CB(hObject, handles)
    fname = get(hObject,'String');
    push_singleWater_CB([], [], handles, fname)

% -----------------------------------------------------------------------------------------
function push_singleWater_CB(hObject, handles, opt)
    if (nargin == 3)        % Direct call
        cd(handles.last_dir)
        [FileName,PathName] = uigetfile({'*.grd;*.GRD', 'Grid files (*.grd,*.GRD)';'*.*', ...
                'All Files (*.*)'},'Select GMT grid');
	    pause(0.01);        cd(handles.home_dir);
        if isequal(FileName,0);     return;     end
        if (PathName ~= 0),         handles.last_dir = PathName;    end
    else        % File name on input
        [PathName,FNAME,EXT] = fileparts(opt);
        PathName = [PathName filesep];      % To be coherent with the 'if' branch
        FileName = [FNAME EXT];
    end
    fname = [PathName FileName];

	[handles,X,Y,handles.Z_water,handles.head_water] = read_gmt_type_grids(handles,fname);
	if (isempty(X)),    return;     end
    
	set(handles.push_clearTestBat,'Visible','on')
	set(handles.edit_singleWater,'String',FileName)
	guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------------------
function edit_namesList_CB(hObject, handles)
    fname = get(hObject,'String');
    push_namesList_CB([], [], handles, fname)

% -----------------------------------------------------------------------------------------
function push_namesList_CB(hObject, handles, opt)
    if (nargin == 3)        % Direct call
        cd(handles.last_dir)
    	str1 = {'*.dat;*.DAT;*.txt;*.TXT', 'Data files (*.dat,*.DAT,*.txt,*.TXT)';'*.*', 'All Files (*.*)'};
        [FileName,PathName] = uigetfile(str1,'File with grids list');
        cd(handles.home_dir);
	    if isequal(FileName,0);     return;     end
        if (PathName ~= 0),         handles.last_dir = PathName;    end
    else        % File name on input
        [PathName,FNAME,EXT] = fileparts(opt);
        PathName = [PathName filesep];      % To be coherent with the 'if' branch
        FileName = [FNAME EXT];
    end
	fname = [PathName FileName];

    [bin,n_column,multi_seg,n_headers] = guess_file(fname);
    % If error in reading file
    if isempty(bin)
        errordlg(['Error reading file ' fname],'Error');    return
    end
    
    fid = fopen([PathName FileName]);
	c = char(fread(fid))';      fclose(fid);
	names = strread(c,'%s','delimiter','\n');   clear c fid;
	m = length(names);
    
    handles.strTimes = [];          % To hold time steps as strings
    if (n_column > 1)
        handles.strTimes = cell(m,1);
        c = false(m,1);
	    for (k=1:m)
            [t,r] = strtok(names{k});
            if (t(1) == '#'),  c(k) = true;  continue;   end
            names{k} = t;
            handles.strTimes{k} = r;
        end
        % Remove eventual commented lines
        if (any(c))
            names(c) = [];          handles.strTimes(c) = [];
            m = length(names);      % Count remaining ones
        end
    end
    
    handles.shortNameList = cell(m,1);      % To hold grid names with path striped
    c = false(m,1);
	for (k=1:m)
        if (n_column == 1 && names{k}(1) == '#')    % If n_column > 1, this test was already done above
            c(k) = true;    continue;
        end
        [PATH,FNAME,EXT] = fileparts(names{k});
        if (isempty(PATH))
            handles.shortNameList{k} = names{k};
            names{k} = [PathName names{k}];
        else
            handles.shortNameList{k} = [FNAME EXT];
        end
        if (any(c))
            names(c) = [];          handles.shortNameList(c) = [];
        end
	end
    
    % Check that at least the files in provided list do exist
    c = false(m,1);
    for (k=1:m)
        c(k) = (exist(names{k},'file') ~= 2);
    end
    names(c) = [];      handles.shortNameList(c) = [];

    handles.nameList = names;
    set(handles.listbox1,'String',handles.shortNameList)
    set(handles.edit_namesList,'String',[PathName FileName])
    guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------------------
function listbox1_CB(hObject, handles)
    % if this is a doubleclick,
    if ( strcmp(get(gcbf,'SelectionType'),'open') && ~isempty(handles.nameList) )
        val = get(hObject,'Value');
        if (~isempty(handles.strTimes))
            handles.testTime = handles.strTimes{val};
        end
        push_singleWater_CB([], [], handles, handles.nameList{val})
    end

% -----------------------------------------------------------------------------------------
function push_clearTestBat_CB(hObject, handles)
    set(handles.edit_singleWater,'String','')
    set(hObject,'Visible','off')
    handles.Z_water = [];
    guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------------------
function edit_movieName_CB(hObject, handles)
    fname = get(hObject,'String');
    push_movieName_CB([], [], handles, fname)

% -----------------------------------------------------------------------------------------
function push_movieName_CB(hObject, handles, opt)
    if (nargin == 3)        % Direct call
        cd(handles.work_dir)
        [FileName,PathName] = uiputfile({'*.gif;*.avi,*.mpg,*.mpeg', ...
                'Grid files (*.gif,*.avi,*.mpg,*.mpeg)'},'Select Movie name');
        pause(0.01);        cd(handles.home_dir);
        if isequal(FileName,0);     return;     end
        if (PathName ~= 0),         handles.last_dir = PathName;    end
        [dumb,FNAME,EXT]= fileparts(FileName);
    	fname = [PathName FileName];
    else        % File name on input
        [PathName,FNAME,EXT] = fileparts(opt);
        PathName = [PathName filesep];      % To be coherent with the 'if' branch
    	fname = opt;
    end
    if (~strmatch(lower(EXT),{'.gif' '.avi' '.mpg' '.mpeg'}))
        errordlg('Ghrrrrrrrr! Don''t be smart. Only ''.gif'', ''.avi'', ''.mpg'' or ''mpeg'' extensions are acepted.', ...
            'Chico Clever');
        return
    end
    
    handles.moviePato = PathName;
    handles.movieName = FNAME;
    if (strcmpi(EXT,'.gif'))
        set(handles.radio_gif,'Value',1)
        radio_gif_CB(handles.radio_gif, [], handles)
    elseif (strcmpi(EXT,'.avi'))
        set(handles.radio_avi,'Value',1)
        radio_avi_CB(handles.radio_avi, [], handles)
    else
        set(handles.radio_avi,'Value',1)
        radio_mpg_CB(handles.radio_mpg, [], handles)
    end
	set(handles.edit_movieName,'String',fname)
    guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------------------
function radio_avi_CB(hObject, handles)
    if (get(hObject,'Value')),      set([handles.radio_gif handles.radio_mpg],'Value',0)
    else                            set(hObject,'Value',1)
    end
    mname = get(handles.edit_movieName,'String');
    if (~isempty(mname))
        mname = [handles.moviePato handles.movieName '.avi'];
        set(handles.edit_movieName,'String',mname)
    end

% -----------------------------------------------------------------------------------------
function radio_gif_CB(hObject, handles)
    if (get(hObject,'Value')),      set([handles.radio_avi handles.radio_mpg],'Value',0)
    else                            set(hObject,'Value',1)
    end
    mname = get(handles.edit_movieName,'String');
    if (~isempty(mname))
        mname = [handles.moviePato handles.movieName '.gif'];
        set(handles.edit_movieName,'String',mname)
    end

% -----------------------------------------------------------------------------------------
function radio_mpg_CB(hObject, handles)
    if (get(hObject,'Value')),      set([handles.radio_avi handles.radio_gif],'Value',0)
    else                            set(hObject,'Value',1)
    end
    mname = get(handles.edit_movieName,'String');
    if (~isempty(mname))
        mname = [handles.moviePato handles.movieName '.mpg'];
        set(handles.edit_movieName,'String',mname)
    end

% -----------------------------------------------------------------------------------------
function edit_fps_CB(hObject, handles)
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
function edit_profile_CB(hObject, handles)
    fname = get(hObject,'String');
    push_profile_CB([], [], handles, fname)

% -----------------------------------------------------------------------------------------
function push_profile_CB(hObject, handles, opt)
    % Read a file with the track coordinates
    if (nargin == 3)        % Direct call
        cd(handles.last_dir)        
        [FileName,PathName] = uigetfile({'*.dat;*.xy', 'Profile file (*.dat,*.xy)';'*.*',...
                'All Files (*.*)'},'Select Profile');
	    pause(0.01);        cd(handles.home_dir);
	    if isequal(FileName,0);     return;     end
        if (PathName ~= 0),         handles.last_dir = PathName;    end
    else        % File name on input
        [PathName,FNAME,EXT] = fileparts(opt);
        PathName = [PathName filesep];      % To be coherent with the 'if' branch
        FileName = [FNAME EXT];
    end
	fname = [PathName FileName];
    
    [bin,n_column,multi_seg,n_headers] = guess_file(fname);
    % If error in reading file
    if isempty(bin) && isempty(n_column) && isempty(multi_seg) && isempty(n_headers)
        errordlg(['Error reading file ' fname],'Error');    return
    end
    if (isa(bin,'struct') || bin ~= 0)   % NOT ASCII
        errordlg('Sorry, reading binary files is not allowed','Error');   return
    end
    if (n_column < 2)
        errordlg('File error. Your file doesn''t have at least 2 columns','Error'); return
    end
    if (multi_seg)
        errordlg('File error. Multisegments file make no sense here.','Error'); return
    end
    if (isempty(n_headers)),    n_headers = NaN;    end
    numeric_data = text_read(fname,NaN,n_headers);
    handles.profile = numeric_data(:,1:2);
	set(handles.edit_profile,'String',fname)
    guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------------------
function edit_imgHeight_CB(hObject, handles)
    h = round(abs(str2double(get(hObject,'String'))));
    if (isnan(h) || h < 50),  return;     end
    handles.imgH = h;
    guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------------------
function edit_imgWidth_CB(hObject, handles)
    w = round(abs(str2double(get(hObject,'String'))));
    if (isnan(w) || w < 50),  return;     end
    handles.imgW = w;
    guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------------------
function pushbutton_OK_CB(hObject, handles)
    % 
    yMax = 0.05;						% Maximum height displayed
    depth_clip = -.07;					% Nealy the maximum depth shown in the plot
    depth_min = depth_clip - 1;			% Now this is the max depth ploted
    depth_min = -0.15;
    if (isempty(handles.profile))
        errordlg('Noooo! Where is the profile file? Do you think I''m bruxo?','ERROR')
        return
    end

    xd = diff(handles.profile(:,1));    yd = diff(handles.profile(:,2));
    tmp = sqrt(xd.*xd + yd.*yd);        rd = zeros(numel(xd)+1,1);
    for i=2:numel(rd),                  rd(i) = rd(i-1) + tmp(i-1);        end
    
    hFig = figure('Pos',[100 100 handles.imgW handles.imgH],'DoubleBuffer','on');
    hAxes = axes('Parent',hFig,'xlim',[rd(1) rd(end)],'ylim',[depth_min yMax]);
    txtPos = [(rd(1) + rd(end)) / 2  yMax * 0.9];

    if (~isempty(handles.Z_bat))
        bot = grdtrack_m(handles.Z_bat,handles.head_bat,handles.profile,'-Z');
        bot(bot < depth_clip) = depth_clip;
        xPB = [rd(1); rd; flipud(rd)];
        yPB = [depth_min; bot; depth_min(ones(numel(bot),1))];
    end
    
    % ------------- If we are testing with one grid only ------------------------
    if (~isempty(handles.Z_water))
        zz = grdtrack_m(handles.Z_water,handles.head_water,handles.profile,'-Z');
        if (isempty(handles.Z_bat))         % We still need to compute X
            xPB = [rd(1); rd; flipud(rd)];
        end
        yy = [depth_min; zz; depth_min(ones(numel(zz),1))];
        patch('xdata',xPB,'ydata',yy,'FaceColor',[0.3098 0.3098 0.8510],'EdgeColor','y','LineWidth',1);
        if (~isempty(handles.Z_bat))
            hold on
            patch('xdata',xPB,'ydata',yPB,'FaceColor',[.4863 .1961 .1961],'EdgeColor','k','LineWidth',1);
        end
        if (~isempty(handles.strTimes))
            text(txtPos(1),txtPos(2),handles.testTime)
        end
        return
    end
    % --------------------------------------------------------------------------

    if (isempty(handles.movieName) && isempty(handles.Z_water))
        errordlg('Hei! what shoult it be the movie name?','ERROR');     return
    end
    
    is_gif = get(handles.radio_gif,'Value');
    is_avi = get(handles.radio_avi,'Value');
    is_mpg = get(handles.radio_mpg,'Value');
    nGrids = numel(handles.nameList);
    if (~isempty(handles.strTimes)),    tempos = handles.strTimes;
    else                                tempos = cell(nGrids,1);     % Just make an empty one
    end
    set(handles.push_stop,'Visible','on')       % To allow interruptions
    for (i=1:nGrids)
        % Check if meanwhile the stop button has been pressed 
        pause(0.01)     % Little pause so that it can listen a eventual STOP request
        handles = guidata(handles.figure1); % We need an updated version
        if (handles.stop)       % Yes. Stop and reset things
            handles.stop = 0;   % Reset
            set(handles.figure1,'Name','1D-Tsunamovie')
            set(handles.push_stop,'Visible','off')
            guidata(handles.figure1,handles)
            break
        end
        
        % Read and interpolate water level grid
        [handles,X,Y,Z,head] = read_gmt_type_grids(handles,handles.nameList{i});
        zz = grdtrack_m(Z,head,handles.profile,'-Z');
        yy = [depth_min; zz; depth_min(ones(numel(zz),1))];
        
        if (i == 1)
            stripgray(hFig, hAxes)
            if (isempty(handles.Z_bat))         % We still need to compute X
                xPB = [rd(1); rd; flipud(rd)];
            end
            h = patch('xdata',xPB,'ydata',yy,'Parent',hAxes,'FaceColor',[0.3098 0.3098 0.8510],'EdgeColor','y','LineWidth',1);
            if (~isempty(tempos{1}))
                hTxt = text(txtPos(1),txtPos(2),tempos{i},'Parent',hAxes);
            end
        else
            set(h,'xdata',xPB,'ydata',yy,'Parent',hAxes)
            if (~isempty(tempos{i}))
                set(hTxt,'String',tempos{i});
            end
        end
        
        if (~isempty(handles.Z_bat))
            if (i > 1),     delete(hPatchBat);      end     % Easier to delete an recreate than uistack
            hPatchBat = patch('xdata',xPB,'ydata',yPB,'Parent',hAxes,'FaceColor',[.4863 .1961 .1961],'EdgeColor','k','LineWidth',1);
        end
        
        %imgWater = imcapture(hFig,'imgAx',[200 400]);
        %imgWater = imcapture(hFig);
        F = getframe(hFig);
        imgWater = F.cdata;
        
        if (is_gif || is_mpg)
            [imgWater,map] = img_fun('rgb2ind',imgWater,8,'nodither');
        end
        
        if (is_gif)
            mname = [handles.moviePato handles.movieName '.gif'];
            if (i == 1)
                writegif(imgWater,map,mname,'loopcount',Inf)
            else
                writegif(imgWater,map,mname,'WriteMode','append','DelayTime',handles.dt)
            end
        elseif (is_avi)        % AVI
            M(i) = im2frame(imgWater);
        else
            M(i) = im2frame(imgWater,map);
        end
    end

    if (is_avi)    % AVI
        mname = [handles.moviePato handles.movieName '.avi'];
  	    movie2avi_j(M,mname,'compression','none','fps',handles.fps)
    elseif (is_mpg)
        mname = [handles.moviePato handles.movieName '.mpg'];
        opt = [1, 0, 1, 0, 10, 5, 5, 5];
  	    mpgwrite(M,map,mname,opt)
    end
    
    set(handles.push_stop,'Visible','off')
    
% -----------------------------------------------------------------------------------------
function push_stop_CB(hObject, handles)
    handles.stop = 1;
    guidata(handles.figure1,handles)
    
% -----------------------------------------------------------------------------------------
function stripgray(hFig, hAxes)

    tenSize = 0;             % Wen axes labels have 10^n this will hold its ~ text height
    axUnit = get(hAxes,'Units');
    axPos = get(hAxes,'pos');           % Save this because we will have to restore it later
    set(hAxes,'Units','Normalized')     % This is normally the default, but be sure
    fig_c = get(hFig,'Color');       set(hFig,'Color','w')

	h_Xlabel = get(hAxes,'Xlabel');         h_Ylabel = get(hAxes,'Ylabel');
	units_save = get(h_Xlabel,'units');
	set(h_Xlabel,'units','pixels');         set(h_Ylabel,'units','pixels');
	Xlabel_pos = get(h_Xlabel,'pos');
	Ylabel_pos = get(h_Ylabel,'Extent');
    
    XTickLabel = get(hAxes,'XTickLabel');    XTick = get(hAxes,'XTick');
    YTickLabel = get(hAxes,'YTickLabel');    YTick = get(hAxes,'YTick');
    if ( str2double(XTickLabel(end,:)) / XTick(end) < 0.1 )
        % We have a 10 power. Note that that's the only way I found to detect
        % the presence of this otherwise completely ghost text.
        tenSize = 20;       % Take into account the 10 power text size
    end
	
	if (abs(Ylabel_pos(1)) < 30)    % Stupid hack, but there is a bug somewhere
        Ylabel_pos(1) = 30;
	end
	
	y_margin = abs(Xlabel_pos(2))+get(h_Xlabel,'Margin') + tenSize;  % To hold the Xlabel height
	x_margin = abs(Ylabel_pos(1))+get(h_Ylabel,'Margin');  % To hold the Ylabel width
	%y_margin = min(max(y_margin,20),30);            % Another hack due to the LabelPos non-sense
    
    figUnit = get(hFig,'Units');        set(hFig,'Units','pixels')
    figPos = get(hFig,'pos');           set(hFig,'Units',figUnit)
    x0 = x_margin / figPos(3);
    y0 = y_margin / figPos(4);
    set(hAxes,'pos',[x0 y0 1-[x0 y0]-1e-2])
    set(h_Xlabel,'units',units_save);     set(h_Ylabel,'units',units_save);
    
    %set(hAxes,'Units',axUnit,'pos',axPos,'Visible',axVis)



% --- Creates and returns a handle to the GUI figure. 
function umDmovie_LayoutFcn(h1)

set(h1,...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'MenuBar','none',...
'Name','umDmovie',...
'NumberTitle','off',...
'PaperSize',[20.98404194812 29.67743169791],...
'Position',[520 471 470 329],...
'Resize','off',...
'HandleVisibility','Call',...
'Tag','figure1');

uicontrol('Parent',h1,'Position',[260 157 181 41],'Style','frame','Tag','frame4');
uicontrol('Parent',h1,'Position',[10 10 221 41],'Style','frame','Tag','frame1');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@main_uiCB,h1,'edit_profile_CB'},...
'HorizontalAlignment','left',...
'Position',[10 289 200 21],...
'Style','edit',...
'Tooltip','Name of a x,y file where to build the transect movie ',...
'Tag','edit_profile');

uicontrol('Parent',h1,...
'Call',{@main_uiCB,h1,'push_profile_CB'},...
'Position',[210 289 21 21],...
'Tooltip','Browse for a profile file name',...
'Tag','push_profile');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@main_uiCB,h1,'edit_namesList_CB'},...
'HorizontalAlignment','left',...
'Position',[10 239 200 21],...
'Style','edit',...
'Tooltip','Name of a file with the water level grids list',...
'Tag','edit_namesList');

uicontrol('Parent',h1,...
'Call',{@main_uiCB,h1,'push_namesList_CB'},...
'Position',[210 239 21 21],...
'Tooltip','Browse for a water level grids list file',...
'Tag','push_namesList');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@main_uiCB,h1,'listbox1_CB'},...
'Position',[10 67 221 171],...
'Style','listbox',...
'Value',1,...
'Tag','listbox1');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@main_uiCB,h1,'edit_singleWater_CB'},...
'HorizontalAlignment','left',...
'Position',[20 15 180 21],...
'Style','edit',...
'Tooltip','Name of one water level grid for test testing purposes',...
'Tag','edit_singleWater');

uicontrol('Parent',h1,...
'Call',{@main_uiCB,h1,'push_singleWater_CB'},...
'Position',[200 15 21 21],...
'Tag','push_singleWater');

uicontrol('Parent',h1,...
'Call',{@main_uiCB,h1,'push_clearTestBat_CB'},...
'FontName','Helvetica',...
'FontSize',7,...
'Position',[181 35 40 16],...
'String','Clear',...
'Tooltip','Remove the test grid so that you can compute a movie',...
'Tag','push_clearTestBat',...
'Visible','off');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@main_uiCB,h1,'edit_movieName_CB'},...
'HorizontalAlignment','left',...
'Position',[251 239 190 21],...
'Style','edit',...
'Tooltip','Name of movie file',...
'Tag','edit_movieName');

uicontrol('Parent',h1,...
'Call',{@main_uiCB,h1,'push_movieName_CB'},...
'Position',[441 239 21 21],...
'Tooltip','Browse for a movie file name (extention is ignored)',...
'Tag','push_movieName');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@main_uiCB,h1,'edit_batGrid_CB'},...
'HorizontalAlignment','left',...
'Position',[251 289 190 21],...
'Style','edit',...
'Tooltip','Name of bathymetry grid',...
'Tag','edit_batGrid');

uicontrol('Parent',h1,...
'Call',{@main_uiCB,h1,'push_batGrid_CB'},...
'Position',[441 289 21 21],...
'Tooltip','Browse for a bathymetry grid',...
'Tag','push_batGrid');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'FontSize',9,...
'HorizontalAlignment','left',...
'Position',[251 311 140 17],...
'String','Bathymetry file (optional)',...
'Style','text',...
'Tag','text1');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'FontSize',9,...
'HorizontalAlignment','left',...
'Position',[11 261 82 17],...
'String','Water files list',...
'Style','text',...
'Tag','text2');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'FontSize',9,...
'Position',[20 42 100 17],...
'String','Test with this file',...
'Style','text',...
'Tag','text_testFile');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'FontSize',9,...
'FontWeight','bold',...
'HorizontalAlignment','left',...
'Position',[251 261 115 17],...
'String','Output movie name',...
'Style','text',...
'Tag','text4');

uicontrol('Parent',h1,'Position',[260 68 181 61],'Style','frame','Tag','frame3');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@main_uiCB,h1,'edit_imgHeight_CB'},...
'Position',[302 165 40 20],...
'String','300',...
'Style','edit',...
'Tooltip','Image height in pixels',...
'Tag','edit_imgHeight');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@main_uiCB,h1,'edit_imgWidth_CB'},...
'Position',[392 165 40 20],...
'String','600',...
'Style','edit',...
'Tooltip','Image width in pixels',...
'Tag','edit_imgWidth');

uicontrol('Parent',h1,...
'Call',{@main_uiCB,h1,'radio_gif_CB'},...
'FontName','Helvetica','Position',[271 97 41 15],...
'String','GIF','Style','radiobutton',...
'Tooltip','Write movie file in animated GIF format',...
'Value',1,'Tag','radio_gif');

uicontrol('Parent',h1,...
'Call',{@main_uiCB,h1,'radio_avi_CB'},...
'FontName','Helvetica','Position',[271 76 41 15],...
'String','AVI','Style','radiobutton',...
'Tooltip','Write movie file in RGB AVI format',...
'Tag','radio_avi');

uicontrol('Parent',h1,...
'FontName','Helvetica','FontSize',9,...
'Position',[272 121 70 15],...
'String','Movie type',...
'Style','text',...
'Tag','text_MovType');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@main_uiCB,h1,'edit_fps_CB'},...
'Position',[400 77 30 18],...
'String','5',...
'Style','edit',...
'Tooltip','Frames per second (ignored in MPEG)',...
'Tag','edit_fps');

uicontrol('Parent',h1,...
'Call',{@main_uiCB,h1,'push_stop_CB'},...
'FontName','Helvetica',...
'FontSize',10,...
'ForegroundColor',[1 0 0],...
'Position',[260 14 66 23],...
'String','STOP',...
'Tag','push_stop',...
'Visible','off');

uicontrol('Parent',h1,...
'Call',{@main_uiCB,h1,'pushbutton_OK_CB'},...
'FontName','Helvetica',...
'FontSize',10,...
'Position',[375 14 66 23],...
'String','OK',...
'Tag','pushbutton_OK');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'HorizontalAlignment','left',...
'Position',[343 79 55 15],...
'String','Frames p/s',...
'Style','text',...
'Tag','text10');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'FontSize',9,...
'HorizontalAlignment','left',...
'Position',[10 311 82 17],...
'String','Profile file',...
'Style','text',...
'Tag','text11');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'FontSize',9,...
'Position',[270 188 100 17],...
'String','Movie frame size',...
'Style','text',...
'Tag','text12');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'HorizontalAlignment','left',...
'Position',[269 167 31 15],...
'String','Height',...
'Style','text',...
'Tag','text13');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'HorizontalAlignment','left',...
'Position',[358 167 30 15],...
'String','Width',...
'Style','text',...
'Tag','text14');

uicontrol('Parent',h1,...
'Call',{@main_uiCB,h1,'radio_mpg_CB'},...
'FontName','Helvetica',...
'Position',[342 97 50 16],...
'String','MPEG',...
'Style','radiobutton',...
'Tooltip','Write movie file in MPEG format',...
'Tag','radio_mpg');


function main_uiCB(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
	feval(callback_name,hObject,guidata(h1));

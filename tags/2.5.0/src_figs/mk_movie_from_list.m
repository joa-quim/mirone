function varargout = mk_movie_from_list(varargin)
% Helper Window to create movies from a list of files

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
	mk_movie_from_list_LayoutFcn(hObject);
	handles = guihandles(hObject);
 
	if (numel(varargin) > 0 && isstruct(varargin{1}))
		handMir = varargin{1};
		handles.work_dir = handMir.work_dir;
		handles.last_dir = handMir.last_dir;
		handles.home_dir = handMir.home_dir;
		handles.hMirFig = handMir.figure1;
		move2side(handMir.figure1, hObject)
	else
		handles.home_dir = cd;
		handles.last_dir = handles.home_dir;
		handles.work_dir = handles.home_dir;
		handles.hMirFig = [];
	end

	f_path = [handles.home_dir filesep 'data' filesep];
	handles.stop = 0;           % When the running engine detects it has canged to 1 it stops
	handles.dither = 'nodither';% Default
	handles.fps = 5;            % Frames per second
	handles.dt = 0.2;           % 1/fps
	handles.barPos = get(handles.frame2,'Pos');		% Kind of waitbar
	
	% Load some icons and put them in the toggles
	load([f_path 'mirone_icons.mat'],'Mfopen_ico');
	set(handles.push_namesList,'CData',Mfopen_ico)
	set(handles.push_movieName,'CData',Mfopen_ico)

	%------------ Give a Pro look (3D) to the frame boxes  -------------------------------
	new_frame3D(hObject, handles.text1, handles.frame1)
	%------------- END Pro look (3D) -----------------------------------------------------

	% Choose default command line output for mk_movie_from_list_export
	handles.output = hObject;
	guidata(hObject, handles);
	
	set(hObject,'Visible','on');
	if (nargout),	varargout{1} = hObject;		end

% -----------------------------------------------------------------------------------------
function edit_namesList_CB(hObject, handles)
    fname = get(hObject,'String');
    push_namesList_CB([], handles, fname)

% -----------------------------------------------------------------------------------------
function push_namesList_CB(hObject, handles, opt)
	if (nargin == 2)		% Direct call
		if (~isempty(handles.hMirFig) && ishandle(handles.hMirFig))
			hand = guidata(handles.hMirFig);
		else
			hand = handles;
		end
		[FileName,PathName] = put_or_get_file(hand, ...
			{'*.dat;*.DAT;*.txt;*.TXT', 'Data files (*.dat,*.DAT,*.txt,*.TXT)';'*.*', 'All Files (*.*)'},'File with grids list','get');
		if isequal(FileName,0),		return,		end

	else        % File name on input
		[PathName,FNAME,EXT] = fileparts(opt);
		PathName = [PathName filesep];      % To be coherent with the 'if' branch
		FileName = [FNAME EXT];
	end
	fname = [PathName FileName];

	[bin,n_column] = guess_file(fname);
	% If error in reading file
	if isempty(bin)
		errordlg(['Error reading file ' fname],'Error'),	return
	end

	fid = fopen(fname);
	c = char(fread(fid))';
	fclose(fid);
	names = strread(c,'%s','delimiter','\n');   clear c fid;
	m = length(names);
	while (isempty(names{m}) && m > 0)
		names(m) = [];
		m = m - 1;
	end
	if (isempty(names))
		errordlg('The list file is ... EMPTY!!!','Error'),	return
	end

	if (n_column > 1)
		c = false(m,1);
		for (k=1:m)
			[t,r] = strtok(names{k});
			if (t(1) == '#'),  c(k) = true;  continue;   end
			names{k} = t;
		end
		% Remove eventual commented lines
		if (any(c))
			names(c) = [];
			m = length(names);			% Count remaining ones
		end
	end

	c = false(m,1);
	for (k=1:m)
		if (n_column == 1 && names{k}(1) == '#')    % If n_column > 1, this test was already done above
			c(k) = true;	continue
		end
		PATH = fileparts(names{k});
		if (isempty(PATH))
			names{k} = [PathName names{k}];
		end
		if (any(c)),	names(c) = [];		end
	end

	% Check that at least the files in provided list do exist
	c = false(m,1);
	for (k=1:m)
		c(k) = (exist(names{k},'file') ~= 2);
	end
	names(c) = [];

	handles.nameList = names;
	set(handles.edit_namesList,'String',[PathName FileName])
	guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------------------
function radio_mpg_CB(hObject, handles)
	if (get(hObject,'Value')),		set([handles.radio_avi handles.radio_gif],'Value',0)
	else							set(hObject,'Value',1)
	end
	mname = get(handles.edit_movieName,'String');
	if (~isempty(mname))
		mname = [handles.moviePato handles.movieName '.mpg'];
		set(handles.edit_movieName,'String',mname)
	end

% -----------------------------------------------------------------------------------------
function radio_avi_CB(hObject, handles)
	if (get(hObject,'Value')),		set([handles.radio_gif handles.radio_mpg],'Value',0)
	else							set(hObject,'Value',1)
	end
	mname = get(handles.edit_movieName,'String');
	if (~isempty(mname))
		mname = [handles.moviePato handles.movieName '.avi'];
		set(handles.edit_movieName,'String',mname)
	end

% -----------------------------------------------------------------------------------------
function radio_gif_CB(hObject, handles)
	if (get(hObject,'Value')),		set([handles.radio_avi handles.radio_mpg],'Value',0)
	else							set(hObject,'Value',1)
	end
	mname = get(handles.edit_movieName,'String');
	if (~isempty(mname))
		mname = [handles.moviePato handles.movieName '.gif'];
		set(handles.edit_movieName,'String',mname)
	end

% -----------------------------------------------------------------------------------------
function checkbox_dither_CB(hObject, handles)
	if (get(hObject,'Value')),		handles.dither = 'dither';
	else							handles.dither = 'nodither';
	end
	guidata(handles.figure1,handles)

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
function edit_movieName_CB(hObject, handles)
	fname = get(hObject,'String');
	push_movieName_CB([], handles, fname)

% -----------------------------------------------------------------------------------------
function push_movieName_CB(hObject, handles, opt)
	if (nargin == 2)		% Direct call
		if (~isempty(handles.hMirFig) && ishandle(handles.hMirFig))
			hand = guidata(handles.hMirFig);
		else
			hand = handles;
		end
		[FileName,PathName] = put_or_get_file(hand, ...
			{'*.gif;*.avi', 'Grid files (*.gif,*.avi)'},'Select Movie name','put');
		if isequal(FileName,0),		return,		end
		[dumb,FNAME,EXT]= fileparts(FileName);
	else					% File name on input
		[PathName,FNAME,EXT] = fileparts(opt);
		if (isempty(PathName))	% If no pato, get the one from the list file
			s = get(handles.edit_namesList,'Str');
			PathName = fileparts(s);
		end
		if (~isempty(PathName))
			PathName = [PathName filesep];		% To be coherent with the 'if' branch
		end
	end
	if ( ~isempty(EXT) && ~any(strcmpi(EXT,{'.gif' '.avi' '.mpg' '.mpeg'})) )
		errordlg('Ghrrrrrrrr! Don''t be smart. Only ''.gif'', ''.avi'', ''.mpg'' or ''mpeg'' extensions are acepted.', ...
			'Chico Clever');
		return
	end

	handles.moviePato = PathName;
	handles.movieName = FNAME;
	if (strcmpi(EXT,'.gif'))
		set(handles.radio_gif,'Value',1)
		radio_gif_CB(handles.radio_gif, handles)
	elseif (strcmpi(EXT,'.avi'))
		set(handles.radio_avi,'Value',1)
		radio_avi_CB(handles.radio_avi, handles)
	elseif (strcmpi(EXT,'.mpg') || strcmpi(EXT,'.mpeg'))
		set(handles.radio_mpg,'Value',1)
		radio_mpg_CB(handles.radio_mpg, handles)
	else
		if (get(handles.radio_gif,'Val'))
			radio_gif_CB(handles.radio_gif, handles)
		elseif (get(handles.radio_avi,'Val'))
			radio_avi_CB(handles.radio_avi, handles)
		else
			radio_mpg_CB(handles.radio_mpg, handles)
		end
	end
	guidata(handles.figure1,handles)
	
% -----------------------------------------------------------------------------------------
function push_OK_CB(hObject, handles)

	if (isempty(handles.nameList))
		errordlg('Where are the images to make the movie?','ERROR'),	return
	end
	if (isempty(handles.movieName))
		errordlg('Hei! what shoult it be the movie name?','ERROR');     return
	end

	is_gif = get(handles.radio_gif,'Value');
	is_avi = get(handles.radio_avi,'Value');
	is_mpg = get(handles.radio_mpg,'Value');

	nFrames = numel(handles.nameList);
	for (k = 1:nFrames)
		img = imread(handles.nameList{k});

		if (is_gif || is_mpg)
			[img, map] = img_fun('rgb2ind',img,256,handles.dither);
		end

		if (is_gif)
			mname = [handles.moviePato handles.movieName '.gif'];
			if (k == 1)
				writegif(img,map,mname,'loopcount',Inf)
			else
				writegif(img,map,mname,'WriteMode','append','DelayTime',handles.dt)
			end
		elseif (is_avi)			% AVI
			M(k) = im2frame(img);
		else					% MPEG
			M(k) = im2frame(img,map);
		end

		% Show visualy the processing advance
		set(handles.frame2,'BackgroundColor',[1 0 0],'Pos', [handles.barPos(1:2) k*handles.barPos(3)/nFrames handles.barPos(4)]);
		pause(0.05)			% Otherwise it won't update
	end

	if (is_avi)
		mname = [handles.moviePato handles.movieName '.avi'];
		movie2avi_j(M,mname,'compression','none','fps',handles.fps)
	elseif (is_mpg)
		mname = [handles.moviePato handles.movieName '.mpg'];
		opt = [1, 0, 1, 0, 10, 5, 5, 5];
		mpgwrite(M,map,mname,opt)
	end

	set(handles.frame2,'BackgroundColor',[0 1 0],'Pos', handles.barPos);		% Reset it to green
	

% --- Creates and returns a handle to the GUI figure. 
function mk_movie_from_list_LayoutFcn(h1)

set(h1,...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'MenuBar','none',...
'Name','Movie from list',...
'NumberTitle','off',...
'Position',[520 562 241 238],...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1');

uicontrol('Parent',h1,...
'Position',[10 37 221 61],...
'Style','frame',...
'Tag','frame1');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@mk_movie_from_list_uiCB,...
'CData',[],...
'HorizontalAlignment','left',...
'Position',[10 169 200 21],...
'Style','edit',...
'Tooltip','Name of a file with the water level grids list',...
'Tag','edit_namesList');

uicontrol('Parent',h1,...
'Call',@mk_movie_from_list_uiCB,...
'Position',[210 169 21 21],...
'Tooltip','Browse for a water level grids list file',...
'Tag','push_namesList');

uicontrol('Parent',h1,...
'Call',@mk_movie_from_list_uiCB,...
'FontName','Helvetica',...
'FontSize',10,...
'Position',[165 8 66 21],...
'String','OK',...
'Tag','push_OK');

uicontrol('Parent',h1,...
'Call',@mk_movie_from_list_uiCB,...
'FontName','Helvetica',...
'Position',[19 46 45 15],...
'String','AVI',...
'Style','radiobutton',...
'Tooltip','Write movie file in RGB AVI format',...
'Tag','radio_avi');

uicontrol('Parent',h1,...
'Call',@mk_movie_from_list_uiCB,...
'FontName','Helvetica',...
'Position',[19 67 45 15],...
'String','GIF',...
'Style','radiobutton',...
'Tooltip','Write movie file in animated GIF format',...
'Value',1,...
'Tag','radio_gif');

uicontrol('Parent',h1,...
'Call',@mk_movie_from_list_uiCB,...
'FontName','Helvetica',...
'Position',[135 70 60 15],...
'String','Dither',...
'Style','checkbox',...
'Tooltip','If you don''t know what this is, ask google',...
'Tag','checkbox_dither');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'FontSize',9,...
'Position',[22 89 70 15],...
'String','Movie type',...
'Style','text',...
'Tag','text1');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@mk_movie_from_list_uiCB,...
'Position',[135 47 30 20],...
'String','5',...
'Style','edit',...
'Tooltip','Frames per second',...
'Tag','edit_fps');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'HorizontalAlignment','left',...
'Position',[168 49 60 15],...
'String','Frames p/s',...
'Style','text');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'FontSize',9,...
'FontWeight','bold',...
'HorizontalAlignment','left',...
'Position',[11 191 110 17],...
'String','Images file list',...
'Style','text');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@mk_movie_from_list_uiCB,...
'CData',[],...
'HorizontalAlignment','left',...
'Position',[10 118 200 21],...
'Style','edit',...
'Tooltip','Name of movie file',...
'Tag','edit_movieName');

uicontrol('Parent',h1,...
'Call',@mk_movie_from_list_uiCB,...
'Position',[210 118 21 21],...
'Tooltip','Browse for a movie file name (extention is ignored)',...
'Tag','push_movieName');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'FontSize',9,...
'FontWeight','bold',...
'HorizontalAlignment','left',...
'Position',[10 141 121 17],...
'String','Output movie name',...
'Style','text');

uicontrol('Parent',h1,...
'Call',@mk_movie_from_list_uiCB,...
'FontName','Helvetica',...
'Position',[73 67 60 16],...
'String','MPEG',...
'Style','radiobutton',...
'Tooltip','Write movie file in mpeg format',...
'Tag','radio_mpg');

uicontrol('Parent',h1,...
'BackgroundColor',[0 1 0],...
'Position',[10 218 221 10],...
'Style','frame',...
'Tag','frame2');

function mk_movie_from_list_uiCB(hObject, eventdata)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));

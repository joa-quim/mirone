function varargout = escadeirar(varargin)
% Interface figure to create grids that have constant values in the nodes between contours
% That is as if they were rice fields on mountain areas.

	hObject = figure('Tag','figure1','Visible','off');
	escadeirar_LayoutFcn(hObject);
	handles = guihandles(hObject);
	movegui(hObject,'east')

	if (numel(varargin) > 0)
		handMir = varargin{1};
		handles.home_dir = handMir.home_dir;
		handles.last_dir = handMir.last_dir;
		handles.work_dir = handMir.work_dir;
		handles.head = handMir.head;
	end

	% copy the next fields into this figure handles to be used in
	% load_grd() as if it was being called with a Mirone handles.
	handles.image_type = handMir.image_type;
	handles.grdname = handMir.grdname;
	handles.computed_grid = handMir.grdname;
	handles.hMirFig1 = handMir.figure1;
	handles.IamCompiled = handMir.IamCompiled;
	
	% -------------- Import/set icons --------------------------------------------
	load([handMir.path_data 'mirone_icons.mat'],'Mfopen_ico');
	set(handles.push_file, 'CData',Mfopen_ico)

	Z = zeros(2);
	Z(1,1) = handles.head(5);
	Z(2,2) = handles.head(6);
	Z([2 3]) = sum(handles.head(5:6)) / 2;
	c = contourc(Z);
	handles.stairs = [handles.head(5) c(1,1:3:end) handles.head(6)];

	handles.stairMin = handles.stairs(1);
	handles.stairMax = handles.stairs(end);
	handles.stairDz = handles.stairs(3)-handles.stairs(2);
	set(handles.edit_min,'String', handles.stairMin)
	set(handles.edit_max,'String', handles.stairMax)
	set(handles.edit_dz,'String', handles.stairDz)

	% -------------- Give a Pro look (3D) to the frame boxes ----------------------- 
	new_frame3D(hObject,NaN, [handles.frame1 handles.frame2 handles.frame3])

	set(hObject,'Visible','on');
	guidata(hObject, handles);
	if (nargout),   varargout{1} = hObject;     end

% -----------------------------------------------------------------------------------------
function radio_MinMax_Callback(hObject, eventdata, handles)
	if (~get(hObject,'Value')),		set(hObject,'Value',1),		return,		end
	set([handles.radio_ML handles.radio_file],'Val',0)
	set([handles.edit_ML handles.edit_file handles.push_file handles.text_MLcmd handles.text_import],'Enable','off')
	set([handles.edit_min handles.edit_dz handles.edit_max],'Enable','on')
	set([handles.text_min handles.text_dz handles.text_max],'Enable','on')

% -----------------------------------------------------------------------------------------
function radio_ML_Callback(hObject, eventdata, handles)
	if (~get(hObject,'Value')),		set(hObject,'Value',1),		return,		end
	set([handles.radio_MinMax handles.radio_file],'Val',0)
	set([handles.edit_min handles.edit_dz handles.edit_max handles.edit_file handles.push_file],'Enable','off')
	set([handles.text_min handles.text_dz handles.text_max handles.text_import],'Enable','off')
	set([handles.edit_ML handles.text_MLcmd],'Enable','on')

% -----------------------------------------------------------------------------------------
function radio_file_Callback(hObject, eventdata, handles)
	if (~get(hObject,'Value')),		set(hObject,'Value',1),		return,		end
	set([handles.radio_MinMax handles.radio_ML],'Val',0)
	set([handles.edit_min handles.edit_dz handles.edit_max handles.edit_ML],'Enable','off')
	set([handles.text_min handles.text_dz handles.text_max handles.text_MLcmd],'Enable','off')
	set([handles.edit_file handles.push_file handles.text_import],'Enable','on')

% -----------------------------------------------------------------------------------------
function edit_min_Callback(hObject, eventdata, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx)),		set(hObject,'String',handles.stairMin),	return,		end
	handles.stairs = escadaria(xx, handles.stairMax, handles.stairDz);
	handles.stairMin = xx;
	guidata(handles.figure1, handles)

% -----------------------------------------------------------------------------------------
function edit_dz_Callback(hObject, eventdata, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx)),		set(hObject,'String',handles.stairDz),	return,		end
	handles.stairs = escadaria(handles.stairMin, handles.stairMax, xx);
	handles.stairDz = xx;
	guidata(handles.figure1, handles)

% -----------------------------------------------------------------------------------------
function edit_max_Callback(hObject, eventdata, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx)),		set(hObject,'String',handles.stairMax),	return,		end
	handles.stairs = escadaria(handles.stairMin, xx, handles.stairDz);
	handles.stairMax = xx;
	guidata(handles.figure1, handles)

% -----------------------------------------------------------------------------------------
function stairs = escadaria(x0, x1, dx)
% Tell contourc to do the job of splitting (x0, x1, dx) in a reasonable way
	n = fix((x1 - x0) / dx);
	Z = zeros(2);
	Z(1,1) = x0;		Z(2,2) = x1;
	Z([2 3]) = (x0 + x1) / 2;
	c = contourc(Z, n);
	stairs = [x0 c(1,1:3:end) x1];

% -----------------------------------------------------------------------------------------
function edit_ML_Callback(hObject, eventdata, handles)

	cmd = get(hObject,'String');
	try
		if (~handles.IamCompiled)
			stairs = eval(cmd);				% This wouldn't work on compiled version
		else
			ind = strfind(str,':');			% Search for a start:inc:end form
			if (numel(ind) ~= 2)
				errordlg('BAD z_min:dz:z_max command','Ignorant'),		return
			end
			start = str2double( str(1:(ind(1)-1)) );
			inc   = str2double( str((ind(1)+1):(ind(2)-1)) );
			fim   = str2double( str((ind(2)+1):end) );
			stairs = start:inc:fim;
		end
	catch
		errordlg(sprintf('Error evaluating command:\n%s',lasterr),'Error'),		return
	end

	stairs = sort(stairs);			% Make sure no inventions
	if (stairs(1) > handles.head(6) || stairs(end) < handles.head(5))
		errordlg('The intervals set by this expression are outside grid limits.','Error')
		set(hObject,'String','')
	else
		handles.stairs = stairs;
		guidata(handles.figure1, handles)
	end

% -----------------------------------------------------------------------------------------
function edit_file_Callback(hObject, eventdata, handles)
	fname = get(hObject,'String');
	push_file_Callback([], [], handles, fname)

% -----------------------------------------------------------------------------------------
function push_file_Callback(hObject, eventdata, handles, opt)
	if (nargin == 3)        % Direct call
		cd(handles.last_dir)
		str1 = {'*.dat;*.DAT;*.txt;*.TXT', 'Data files (*.dat,*.DAT,*.txt,*.TXT)';'*.*', 'All Files (*.*)'};
		[FileName,PathName] = uigetfile(str1,'File with step widths');
		cd(handles.home_dir);
	if isequal(FileName,0),		return,		end
	if (PathName ~= 0),			handles.last_dir = PathName;    end
	else        % File name on input
		[PathName,FNAME,EXT] = fileparts(opt);
		PathName = [PathName filesep];      % To be coherent with the 'if' branch
		FileName = [FNAME EXT];
	end
	fname = [PathName FileName];
    set(handles.edit_file,'String',fname)

	[bin,n_column,multi_seg,n_headers] = guess_file(fname);
	% If error in reading file
	if ( isempty(bin) ),	errordlg(['Error reading file ' fname],'Error'),	return,		end
	if (bin ~= 0),			errordlg('Sorry, reading binary files is not programed','Error'),	return,		end
	if (n_column < 1),		errordlg('This file has nothing usefull','Error'),	return,	end

	[PATH,FNAME,EXT] = fileparts(fname);
	if (strcmpi(EXT,'.cpt'))
		[cmap,z_intervals] = cpt2cmap(['-C' fname]);
		z_intervals = [z_intervals(:,1); z_intervals(end)];			% It was a 2 column array
		if (z_intervals(1) == 0 && z_intervals(end) == 1 )			% A master palette
			stairs = handles.head(5) + z_intervals * (handles.head(6) - handles.head(5));
		end
	else
		if (isempty(n_headers)),    n_headers = NaN;    end
		stairs = text_read(fname,NaN,n_headers);
	end

	if (size(stairs,2) > 1),		stairs = stairs(:,1);	end
	stairs = sort(stairs);			% Make sure no inventions

	if (stairs(1) > handles.head(6) || stairs(end) < handles.head(5))
		errordlg('The intervals in your file are idiot if intended to be applyied to current grid.','Error')
		set(hObject,'String','')
	else
		handles.stairs = stairs;
		guidata(handles.figure1, handles)
	end

% -----------------------------------------------------------------------------------------
function push_ok_Callback(hObject, eventdata, handles)

	handles.figure1 = handles.hMirFig1;		% it won't hurt if we don't save the handles
	[X,Y,Z,head] = load_grd(handles);
	if isempty(Z),  return,		end			% An error message was already issued

	Z_stairs = alloc_mex(size(Z,1), size(Z,2), 'single');
	for (k = 2:numel(handles.stairs))
		ind = (Z >= handles.stairs(k-1) & Z < handles.stairs(k) );		% Indices of a slice 
		Z(ind) = handles.stairs(k-1);
		Z_stairs = (Z_stairs | Z);
	end

	tmp.X = X;		tmp.Y = Y;		tmp.head = head;
	tmp.name = 'Escovinha';
	mirone(Z, tmp)

% -----------------------------------------------------------------------------------------
function new_frame3D(hFig,hText,hFrame)

	% Give a Pro look (3D) to the frame boxes 
	bgcolor = get(0,'DefaultUicontrolBackgroundColor');
	framecolor = max(min(0.65*bgcolor,[1 1 1]),[0 0 0]);
	if (nargin < 3)
		hFrame = findobj(hFig,'Style','Frame');
	end
	for (i = 1:numel(hFrame))
        frame_size = get(hFrame(i),'Position');
        f_bgc = get(hFrame(i),'BackgroundColor');
        usr_d = get(hFrame(i),'UserData');
		frame3D(hFig,frame_size,framecolor,f_bgc,usr_d)
		delete(hFrame(i))
	end

	if (isnan(hText)),		return,		end

	% Recopy the text fields on top of previously created frames (uistack is too slow)
	if (isempty(hText))
		hText = findobj(hFig,'Style','Text');
	end
	for (i = 1:numel(hText))
        usr_d = get(hText(i),'UserData');
        t_size = get(hText(i),'Position');   t_str = get(hText(i),'String');    fw = get(hText(i),'FontWeight');
        fn = get(hText(i),'FontName');	bgc = get (hText(i),'BackgroundColor');   fgc = get (hText(i),'ForegroundColor');
        uicontrol('Parent',hFig, 'Style','text', 'Position',t_size,'String',t_str, ...
            'BackgroundColor',bgc,'ForegroundColor',fgc,'FontWeight',fw, 'FontName',fn, 'UserData',usr_d);
	end
	delete(hText)


% --- Creates and returns a handle to the GUI figure. 
function escadeirar_LayoutFcn(h1,handles);

set(h1,...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'MenuBar','none',...
'Name','Escadeirar',...
'NumberTitle','off',...
'Position',[520 608 321 185],...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1');

uicontrol('Parent',h1,'Position',[10 35 301 51],'Style','frame','Tag','frame3');
uicontrol('Parent',h1,'Position',[10 86 301 46],'Style','frame','Tag','frame2');
uicontrol('Parent',h1,'Position',[10 132 301 48],'Style','frame','Tag','frame1');

uicontrol('Parent',h1,...
'Callback',{@escadeirar_uicallback,h1,'radio_MinMax_Callback'},...
'Position',[20 146 15 15],...
'Style','radiobutton',...
'Value',1,...
'Tag','radio_MinMax');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@escadeirar_uicallback,h1,'edit_min_Callback'},...
'Position',[40 142 81 21],...
'Style','edit',...
'TooltipString','Grid minimum',...
'Tag','edit_min');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@escadeirar_uicallback,h1,'edit_dz_Callback'},...
'Position',[140 142 61 21],...
'Style','edit',...
'TooltipString','"Stair" height',...
'Tag','edit_dz');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@escadeirar_uicallback,h1,'edit_max_Callback'},...
'Position',[220 142 81 21],...
'Style','edit',...
'TooltipString','Grid maximum',...
'Tag','edit_max');

uicontrol('Parent',h1,...
'Callback',{@escadeirar_uicallback,h1,'radio_ML_Callback'},...
'Position',[20 100 15 15],...
'Style','radiobutton',...
'Tag','radio_ML');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@escadeirar_uicallback,h1,'edit_ML_Callback'},...
'Enable','off',...
'Position',[140 97 161 21],...
'Style','edit',...
'TooltipString','An expression of the type [zmin:inc:zmax] or for example [-10 -3 7 10]',...
'Tag','edit_ML');

uicontrol('Parent',h1,...
'Callback',{@escadeirar_uicallback,h1,'radio_file_Callback'},...
'Position',[20 49 15 15],...
'Style','radiobutton',...
'Tag','radio_file');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@escadeirar_uicallback,h1,'edit_file_Callback'},...
'Enable','off',...
'HorizontalAlignment','left',...
'Position',[40 45 241 21],...
'Style','edit',...
'TooltipString','Enter file name here',...
'Tag','edit_file');

uicontrol('Parent',h1,...
'Callback',{@escadeirar_uicallback,h1,'push_file_Callback'},...
'Enable','off',...
'Position',[281 45 21 21],...
'Tag','push_file');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'Position',[54 163 51 15],...
'String','Min',...
'Style','text',...
'Tag','text_min');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'Position',[240 163 51 15],...
'String','Max',...
'Style','text',...
'Tag','text_max');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'Position',[146 163 51 15],...
'String','Dz',...
'Style','text',...
'Tag','text_dz');

uicontrol('Parent',h1,...
'Enable','off',...
'FontName','Helvetica',...
'Position',[41 100 95 16],...
'String','Matlab expression',...
'Style','text',...
'Tag','text_MLcmd');

uicontrol('Parent',h1,...
'Enable','off',...
'FontName','Helvetica',...
'Position',[40 67 120 16],...
'String','Import "stairs" from file',...
'Style','text',...
'Tag','text_import');

uicontrol('Parent',h1,...
'Callback',{@escadeirar_uicallback,h1,'push_ok_Callback'},...
'FontName','Helvetica',...
'FontSize',9,...
'Position',[166 8 66 23],...
'String','OK',...
'Tag','push_ok');

uicontrol('Parent',h1,...
'Callback',{@escadeirar_uicallback,h1,'delete(gcf)'},...
'FontName','Helvetica',...
'FontSize',9,...
'Position',[246 8 66 23],...
'String','Cancel',...
'Tag','push_cancel');

function escadeirar_uicallback(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
feval(callback_name,hObject,[],guidata(h1));

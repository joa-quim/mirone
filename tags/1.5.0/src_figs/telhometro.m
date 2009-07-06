function varargout = telhometro(varargin)
% M-File changed by desGUIDE 

	hObject = figure('Tag','figure1','Visible','off');
	telhometro_LayoutFcn(hObject);
	handles = guihandles(hObject);
	movegui(hObject,'center');

	%#function telha_m choosebox

	if (~isempty(varargin))
		handles.h_calling_fig = varargin{1};
		handles.mironeAxes = get(varargin{1},'CurrentAxes'); 
		if (length(varargin) == 2)          % Called with the line handle in argument
			handles.h_line_orig = varargin{2};
			set(handles.text_activeLine,'String','GOT A LINE TO WORK WITH')
		end
	else
		errordlg('TELHOMETRO: wrong number of arguments.','Error')
		delete(hObject)
		return
	end

	handMir = guidata(handles.h_calling_fig);
	handles.path_data = handMir.path_data;
	handles.path_tmp = handMir.path_tmp;
	handles.path_continent = [handMir.home_dir filesep 'continents' filesep];
	handles.h_line_orig = [];
	handles.p_lon = [];
	handles.p_lat = [];
	handles.p_omega = [];

	guidata(hObject, handles);

	set(hObject,'Visible','on');
	if (nargout),   varargout{1} = hObject;     end

% -------------------------------------------------------------------------------------
function edit_polesFile_Callback(hObject, eventdata, handles)
	fname = get(hObject,'String');
	if isempty(fname)    return;    end
	% Let the pushbutton_readPolesFile_Callback do all the work
	pushbutton_readPolesFile_Callback(hObject,[],guidata(gcbo),fname)

% -------------------------------------------------------------------------------------
function pushbutton_readPolesFile_Callback(hObject, eventdata, handles, opt)
% Get poles file name
	if (nargin == 4)    fname = opt;
	else                opt = [];
	end

	if (isempty(opt))           % Otherwise we already know fname from the 4th input argument
		handMir = guidata(handles.h_calling_fig);
		str1 = {'*.stg;*.dat;*.DAT', 'Data files (*.stg,*.dat,*.DAT)';'*.*', 'All Files (*.*)'};
		[FileName,PathName] = put_or_get_file(handMir,str1,'Select poles file','get');
		if isequal(FileName,0),		return,		end
        fname = [PathName FileName];
	end
	set(handles.edit_polesFile,'String',fname)

% -----------------------------------------------------------------------------------
function edit_timeStart_Callback(hObject, eventdata, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx) || xx < 0)     set(hObject, 'String','0');     end

% -----------------------------------------------------------------------------------
function edit_timeEnd_Callback(hObject, eventdata, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx) || xx < 0)     set(hObject, 'String','');      end

% -------------------------------------------------------------------------------------
function pushbutton_compute_Callback(hObject, eventdata, handles)

	if (isempty(handles.h_line_orig))
        errordlg('Will you be so kind to let me know what line should I rotate?','Unknown target')
        return
	end
	if (get(handles.check_additionalRot,'Val') && (isempty(handles.p_lon) || isempty(handles.p_lat) || isempty(handles.p_omega)) )
        warndlg('Additional rotation is not completely deffined. Don''t you see?','WARNING')
		set(handles.check_additionalRot,'Val',0)
		set([handles.edit_poleLon handles.edit_poleLat handles.edit_poleAng],'Enable','off')
	end

	handMir = guidata(handles.h_calling_fig);
	lt = handMir.DefLineThick;
	lc = handMir.DefLineColor;

	poles_name = get(handles.edit_polesFile,'String');
	if (isempty(poles_name))
        errordlg('No stage poles provided','Error');    return
	end

	axes(handles.mironeAxes)       % Make the Mirone axes the CurrentAxes
	x = get(handles.h_line_orig,'XData');       y = get(handles.h_line_orig,'YData');
	linha = [x(:) y(:)];
	opt_E = ['-E' poles_name];
	if (get(handles.checkbox_revertRot,'Value'))
        opt_I = '-I';
	else
        opt_I = ' ';
	end

	t0 = str2double(get(handles.edit_timeStart,'String'));
	if (t0 > 0)     opt_T = ['-T' num2str(t0)];
	else            opt_T = ' ';
	end

	t1 = str2double(get(handles.edit_timeEnd,'String'));
	if (t1 > 0)     opt_N = ['-N' num2str(t1)];
	else            opt_N = ' ';
	end

	if (t1 < t0)
        errordlg('You are a bit confused with ages. Time End is < Time Start ','Error')
        return
	end

	[out_x,out_y,first_mag,n_flow] = telha_m(linha, opt_E, opt_I, '-B', opt_T, opt_N);

	% if (get(handles.checkbox_ridgeStartTime,'Value'))
	% 	dx = out_x(1,1) - linha(1,1);
	% 	dy = out_y(1,1) - linha(1,2);
	% 	out_x = out_x - dx;
	% 	out_y = out_y - dy;
	% end

	if (get(handles.checkbox_ridgeStartTime,'Value') && t0 > 0)
		% We use a trick here. The -T option of telha only makes it return telhas after age t0, but it still
		% assumes that the line it received was a chunk of Ridge. So the trick is to compute a pole that brings
		% back the first returned telha to the seed line provided in input. 
		% This way, by rotating the whole carpet by that angle, we make like the carpet was beeing generated
		% from that isochron back to old eternity.
        x1 = [linha(1,1) linha(2,1)];       y1 = [linha(1,2) linha(2,2)];		% Lon, Lat vertex of the seed line
        x2 = [out_x(1,1) out_x(2,1)];       y2 = [out_y(1,1) out_y(2,1)];		% Lon, Lat vertex of telha whose age is t0
        [p_lon,p_lat,omega] = calc_bonin_euler_pole (x1,y1,x2,y2);
        [rlon,rlat] = rot_euler(out_x(:),out_y(:),p_lon,p_lat,-omega);
        % Now we have to rebuild the "telhas" matrix format
        n_col = size(out_x,2); 
        out_x = reshape(rlon,4,n_col);
        out_y = reshape(rlat,4,n_col);
        clear rlon rlat;
	end

	% Save this for evental use in write_script
	to_save.line = linha;
	to_save.opt_E = opt_E;      to_save.opt_I = opt_I;
	to_save.opt_T = opt_T;      to_save.opt_N = opt_N;

	% Remove last column, it has only zeros
	out_x(:,end) = [];              out_y(:,end) = [];
	% Split between direct and inverse telhas
	out_x1 = out_x(:,1:2:end);      out_x2 = out_x(:,2:2:end);
	out_y1 = out_y(:,1:2:end);      out_y2 = out_y(:,2:2:end);
	if (first_mag < 0)      % Negative magnetizations correspond to direct polarities in telha.h
        cor_d = 'r';    cor_r = 'b';
	else
        cor_d = 'b';    cor_r = 'r';
	end

	if (get(handles.check_additionalRot,'Val'))
		dims1 = size(out_x1);		dims2 = size(out_x2);
        [out_x1,out_y1] = rot_euler(out_x1,out_y1,handles.p_lon,handles.p_lat,handles.p_omega);
        [out_x2,out_y2] = rot_euler(out_x2,out_y2,handles.p_lon,handles.p_lat,handles.p_omega);
		% Need to reshape because they are all now column vectors
		out_x1=reshape(out_x1,dims1);		out_y1=reshape(out_y1,dims1);
		out_x2=reshape(out_x2,dims2);		out_y2=reshape(out_y2,dims2);
	end

	hd = patch(out_x1,out_y1,cor_d,'EdgeColor','none','Tag','tapete');
	hr = patch(out_x2,out_y2,cor_r,'EdgeColor','none','Tag','tapete_R');
	set(hd,'UserData',to_save)

	%set_telhas_uis(hd);     set_telhas_uis(hr)
	draw_funs(hd,'telhas_patch')
	draw_funs(hr,'telhas_patch')

% --------------------------------------------------------------------
function set_telhas_uis(h)
	cmenuHand = uicontextmenu;
	set(h, 'UIContextMenu', cmenuHand);
	uimenu(cmenuHand, 'Label', 'Delete', 'Callback', 'delete(gco)');

% --------------------------------------------------------------------
function pushbutton_callMagBarCode_Callback(hObject, eventdata, handles)
	MagBarCode([handles.path_data 'Cande_Kent_95.dat'])

% --------------------------------------------------------------------
function pushbutton_polesList_Callback(hObject, eventdata, handles)
	fid = fopen([handles.path_continent 'lista_polos.dat'],'rt');
	c = fread(fid,'*char').';
	fclose(fid);
	s = strread(c,'%s','delimiter','\n');

	[s,v] = choosebox('Name','One Euler list',...
                        'PromptString','List of poles:',...
                        'SelectString','Selected poles:',...
                        'ListSize',[380 300],...
                        'ListString',s);

	if (v == 1)         % Finite pole
		if (~get(handles.check_additionalRot,'Val'))
			warndlg('In order to give use to the single pole selection, the "Apply additional rotation" option must be checked. Ignoring selected pole.','Warning')
			return
		end
		handles.p_lon = s(1);
		handles.p_lat = s(2);
		handles.p_omega = s(3);
		set(handles.edit_poleLon, 'String', s(1))
		set(handles.edit_poleLat, 'String', s(2))
		set(handles.edit_poleAng, 'String', s(3))
		guidata(hObject,handles)
	elseif (v == 2)     % Stage poles
		set(handles.edit_polesFile,'String',s)
	end

% -----------------------------------------------------------------------------------
function pushbutton_pickLine_Callback(hObject, eventdata, handles)
	if (get(hObject,'Value'))
        % Test if we have potential target lines and their type
        h_mir_lines = findobj(handles.h_calling_fig,'Type','line');     % Fish all objects of type line in Mirone figure
        if (isempty(h_mir_lines))                                       % We don't have any lines
            str = ['If you hited this button on purpose, than you deserve the following insult.',...
                    'You #!|"*!%!?~^)--$&.',... 
                    'THERE ARE NO LINES IN THAT FIGURE.'];
            errordlg(str,'Chico Clever');   set(hObject,'Value',0);     return;
        end
        % The above test is not enough. For exemple, coastlines are not eligible neither,
        % but is very cumbersome to test all the possibilities of pure non-eligible lines.
        set(handles.h_calling_fig,'pointer','crosshair')
        h_line = get_polygon(handles.h_calling_fig);          % Get the line handle
        handles.h_line_orig = h_line;
        if (~isempty(h_line))
            % Create a empty line handle that will hold the rotated line
            handles.h_line = line('parent',get(handles.h_calling_fig,'CurrentAxes'),'XData',[],'YData',[], ...
                'LineStyle','-.','LineWidth',2);
            set(handles.text_activeLine,'String','GOT A LINE TO WORK WITH','ForegroundColor',[0 0.6 0])
        else
            handles.h_line_orig = [];
            set(handles.text_activeLine,'String','NO ACTIVE LINE','ForegroundColor',[1 0 0])
            set(hObject,'Value',0)
        end
        set(handles.h_calling_fig,'pointer','arrow')
        set(hObject,'Value',0)
        figure(handles.figure1)                 % Bring this figure to front again
	else        % What should I do?
	end
	guidata(hObject, handles);

% -----------------------------------------------------------------------------------
function edit_poleLon_Callback(hObject, eventdata, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx))
		set(hObject,'String','') 
		handles.p_lon = [];
		return
	end
	handles.p_lon = xx;
	guidata(hObject,handles)

% -----------------------------------------------------------------------------------
function edit_poleLat_Callback(hObject, eventdata, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx))
		set(hObject,'String','')
		handles.p_lat = [];
		return
	end
	handles.p_lat = xx;
	guidata(hObject,handles)

% -----------------------------------------------------------------------------------
function edit_poleAng_Callback(hObject, eventdata, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx))
		set(hObject,'String','')
		handles.p_omega = [];
		return
	end
	handles.p_omega = xx;
	guidata(hObject,handles)

% -----------------------------------------------------------------------------------
function check_additionalRot_Callback(hObject, eventdata, handles)
	if (get(hObject,'Value'))
		set([handles.edit_poleLon handles.edit_poleLat handles.edit_poleAng],'Enable','on')
	else
		set([handles.edit_poleLon handles.edit_poleLat handles.edit_poleAng],'Enable','off')
	end

% -----------------------------------------------------------------------------------
function figure1_KeyPressFcn(hObject, eventdata)
	handles = guidata(hObject);
	if isequal(get(hObject,'CurrentKey'),'escape')
		delete(handles.figure1);
	end


% --- Creates and returns a handle to the GUI figure. 
function telhometro_LayoutFcn(h1)
set(h1,...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'KeyPressFcn',@figure1_KeyPressFcn,...
'MenuBar','none',...
'Name','Telhometro',...
'NumberTitle','off',...
'Position',[520 570 441 221],...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1');

uicontrol('Parent',h1,...
'Position',[10 10 231 71],...
'Style','frame',...
'Tag','frame1');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@telhometro_uicallback,h1,'edit_polesFile_Callback'},...
'HorizontalAlignment','left',...
'Position',[10 141 251 21],...
'Style','edit',...
'Tag','edit_polesFile');

uicontrol('Parent',h1,...
'Callback',{@telhometro_uicallback,h1,'pushbutton_readPolesFile_Callback'},...
'FontSize',10,...
'FontWeight','bold',...
'Position',[260 140 21 21],...
'String','...',...
'Tag','pushbutton_readPolesFile');

uicontrol('Parent',h1,...
'HorizontalAlignment','left',...
'Position',[10 165 90 15],...
'String','Stage poles file',...
'Style','text');

uicontrol('Parent',h1,...
'Callback',{@telhometro_uicallback,h1,'pushbutton_compute_Callback'},...
'FontSize',9,...
'FontWeight','bold',...
'Position',[351 9 80 21],...
'String','Compute',...
'Tag','pushbutton_compute');

uicontrol('Parent',h1,...
'Position',[290 143 145 15],...
'String','Revert sense of rotation',...
'Style','checkbox',...
'TooltipString','Revert the sense of rotation defined by the stages poles',...
'Tag','checkbox_revertRot');

uicontrol('Parent',h1,...
'Callback',{@telhometro_uicallback,h1,'pushbutton_callMagBarCode_Callback'},...
'Position',[300 59 131 23],...
'String','Magnetic Bar Code',...
'TooltipString','Open the magnetic bar code window',...
'Tag','pushbutton_callMagBarCode');

uicontrol('Parent',h1,...
'Callback',{@telhometro_uicallback,h1,'pushbutton_polesList_Callback'},...
'Position',[300 99 131 23],...
'String','Poles selector',...
'TooltipString','Construct a stage poles file',...
'Tag','pushbutton_polesList');

uicontrol('Parent',h1,...
'Callback',{@telhometro_uicallback,h1,'pushbutton_pickLine_Callback'},...
'Position',[10 189 161 23],...
'String','Pick line from Figure',...
'TooltipString','Allows you to mouse select one line from a Mirone figure',...
'Tag','pushbutton_pickLine');

uicontrol('Parent',h1,...
'FontSize',12,...
'FontWeight','bold',...
'ForegroundColor',[1 0 0],...
'Position',[211 191 220 18],...
'HorizontalAlignment','right',...
'String','NO ACTIVE LINE',...
'Style','text',...
'Tag','text_activeLine');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@telhometro_uicallback,h1,'edit_timeStart_Callback'},...
'Position',[10 92 41 21],...
'String','0',...
'Style','edit',...
'Tag','edit_timeStart');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@telhometro_uicallback,h1,'edit_timeEnd_Callback'},...
'Position',[80 92 41 21],...
'Style','edit',...
'Tag','edit_timeEnd');

uicontrol('Parent',h1,...
'Position',[7 115 51 15],...
'String','Start time',...
'Style','text');

uicontrol('Parent',h1,...
'Position',[76 115 51 15],...
'String','Stop time',...
'Style','text');

uicontrol('Parent',h1,...
'Position',[130 92 111 15],...
'String','Ridge at Start time',...
'Style','checkbox',...
'TooltipString','Use when you want and old (Start time) ridge position',...
'Tag','checkbox_ridgeStartTime');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@telhometro_uicallback,h1,'edit_poleLon_Callback'},...
'Position',[30 21 61 21],...
'Style','edit',...
'TooltipString','Longitude of the first Euler pole',...
'Tag','edit_poleLon');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@telhometro_uicallback,h1,'edit_poleLat_Callback'},...
'Position',[100 21 61 21],...
'Style','edit',...
'TooltipString','Latitude of the first Euler pole',...
'Tag','edit_poleLat');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@telhometro_uicallback,h1,'edit_poleAng_Callback'},...
'Position',[170 21 61 21],...
'Style','edit',...
'TooltipString','Angle of rotation of first pole',...
'Tag','edit_poleAng');

uicontrol('Parent',h1, 'Position',[42 44 41 15],...
'String','Lon', 'Style','text');

uicontrol('Parent',h1, 'Position',[113 44 41 15],...
'String','Lat', 'Style','text');

uicontrol('Parent',h1, 'Position',[184 44 41 15],...
'String','Angle', 'Style','text');

uicontrol('Parent',h1,...
'Callback',{@telhometro_uicallback,h1,'check_additionalRot_Callback'},...
'Position',[17 61 150 15],...
'String','Apply additional rotation',...
'Style','checkbox',...
'TooltipString','Apply an additional rotation to the "Stage Blanket"',...
'Tag','check_additionalRot');


function telhometro_uicallback(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
feval(callback_name,hObject,[],guidata(h1));

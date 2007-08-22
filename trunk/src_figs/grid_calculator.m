function varargout = grid_calculator(varargin)
% M-File changed by desGUIDE 
% varargin   command line arguments to grid_calculator (see VARARGIN)

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
	grid_calculator_LayoutFcn(hObject,handles);
	handles = guihandles(hObject);
	movegui(hObject,'east')

	old_vis = get(0,'ShowHiddenHandles');
	set(0,'ShowHiddenHandles','on')
	h_figs = findobj('Type','figure');
	set(0,'ShowHiddenHandles',old_vis)

	handles.grid_patos = [];    % To hold paths of eventual future loaded grids
	handles.loaded_grid = [];   % To hold names of eventual future loaded grids
	handles.name_str = [];
	handles.BL = [];

	if (numel(varargin) >= 1 && ishandle(varargin{1}) )     % Na verdade tem que testar se FIG
		BL = getappdata(varargin{1},'BandList');
		if ( ~isempty(BL) )
            for i=1:ndims(BL{2})        %
                handles.name_str{i} = ['Band_' num2str(i)];
            end
            h_figs = [];        % Tear off the net for the 'h_figs' fishing
            handles.BL = BL{2};
		end
	end

	% Fish whatever arrays are in memory now (hopefully)
	if (~isempty(h_figs))
        n = 1;
        for (i=1:length(h_figs))
            hand_fig = guidata(h_figs(i));
            % Use a try->catch because ML is too dumb to deal correctly with killed figures
            try,    Z = getappdata(hand_fig.figure1,'dem_z');
            catch   Z = [];
            end
            if (isempty(Z))     continue;   end
            name = get(h_figs(i),'Name');
            ind = strfind(name,' @ ');
            if (~isempty(ind))
                name = name(1:ind-1);
            end
            [pathstr,name,ext] = fileparts(name);
            handles.name_str{n} = [name ext];
            handles.h_figs(n) = h_figs(i);      % Save the figure handles
            n = n + 1;
        end
	end

	% Fill the listbox with the names of the in-memory arrays
	if (~isempty(handles.name_str))
        set(handles.listbox_inArrays,'String',handles.name_str)
	end

	% Choose default command line output for grid_calculator_export
	handles.output = hObject;
	guidata(hObject, handles);

	set(hObject,'Visible','on');
	if (nargout),   varargout{1} = hObject;     end

% ------------------------------------------------------------------------
function edit_command_Callback(hObject, eventdata, handles)
    % Nothing to do here

% ------------------------------------------------------------------------
function listbox_inArrays_Callback(hObject, eventdata, handles)
	% if this is a doubleclick, copy the array's name to the edit box
	if strcmp(get(gcbf,'SelectionType'),'open')
        str = get(hObject,'String');
        sel = get(hObject,'Value');
        com = get(handles.edit_command,'String');
        set(handles.edit_command,'String', [com ' &' str{sel}])
	end

% ------------------------------------------------------------------------
function pushbutton_1_Callback(hObject, eventdata, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') '1'])

% ------------------------------------------------------------------------
function pushbutton_2_Callback(hObject, eventdata, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') '2'])

% ------------------------------------------------------------------------
function pushbutton_3_Callback(hObject, eventdata, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') '3'])

% ------------------------------------------------------------------------
function pushbutton_4_Callback(hObject, eventdata, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') '4'])

% ------------------------------------------------------------------------
function pushbutton_5_Callback(hObject, eventdata, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') '5'])

% ------------------------------------------------------------------------
function pushbutton_6_Callback(hObject, eventdata, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') '6'])

% ------------------------------------------------------------------------
function pushbutton_7_Callback(hObject, eventdata, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') '7'])

% ------------------------------------------------------------------------
function pushbutton_8_Callback(hObject, eventdata, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') '8'])

% ------------------------------------------------------------------------
function pushbutton_9_Callback(hObject, eventdata, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') '9'])

% ------------------------------------------------------------------------
function pushbutton_0_Callback(hObject, eventdata, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') '0'])

% ------------------------------------------------------------------------
function pushbutton_dot_Callback(hObject, eventdata, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') '.'])

% ------------------------------------------------------------------------
function pushbutton_equal_Callback(hObject, eventdata, handles)
	%set(handles.edit_command,'String', [get(handles.edit_command,'String') ' = '])
	pushbutton_compute_Callback(hObject, eventdata, handles)

% ------------------------------------------------------------------------
function pushbutton_devide_Callback(hObject, eventdata, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') ' / '])

% ------------------------------------------------------------------------
function pushbutton_mull_Callback(hObject, eventdata, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') ' * '])

% ------------------------------------------------------------------------
function pushbutton_minus_Callback(hObject, eventdata, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') ' - '])

% ------------------------------------------------------------------------
function pushbutton_plus_Callback(hObject, eventdata, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') ' + '])

% ------------------------------------------------------------------------
function pushbutton_loadGrid_Callback(hObject, eventdata, handles)
	% This function doesn't realy loads the grid. It only stores the grid name.
	% True loading is donne in "Compute"
	str1 = {'*.grd;*.GRD;*.nc;*.NC', 'Grid files (*.grd,*.GRD,*.nc,*.NC)';'*.*', 'All Files (*.*)'};
	[FileName,PathName] = uigetfile(str1,'Select GMT grid');
	if isequal(FileName,0);     return;     end
	str = get(handles.listbox_inArrays,'String');
	str{end+1} = FileName;
	set(handles.listbox_inArrays,'String',str);
	handles.grid_patos{end+1} = PathName;
	handles.loaded_grid{end+1} = FileName;
	guidata(hObject,handles)

% ------------------------------------------------------------------------
function pushbutton_sin_Callback(hObject, eventdata, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') ' sin( '])

% ------------------------------------------------------------------------
function pushbutton_cos_Callback(hObject, eventdata, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') ' cos( '])

% ------------------------------------------------------------------------
function pushbutton_tan_Callback(hObject, eventdata, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') ' tan( '])

% ------------------------------------------------------------------------
function pushbutton_log10_Callback(hObject, eventdata, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') ' log10( '])

% ------------------------------------------------------------------------
function pushbutton_log_e_Callback(hObject, eventdata, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') ' log( '])

% ------------------------------------------------------------------------
function pushbutton_exp_Callback(hObject, eventdata, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') ' exp( '])

% ------------------------------------------------------------------------
function pushbutton_sqrt_Callback(hObject, eventdata, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') ' sqrt( '])

% ------------------------------------------------------------------------
function pushbutton_abs_Callback(hObject, eventdata, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') ' abs( '])

% ------------------------------------------------------------------------
function pushbutton_leftPar_Callback(hObject, eventdata, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') '( '])

% ------------------------------------------------------------------------
function pushbutton_rightPar_Callback(hObject, eventdata, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') ' )'])

% ------------------------------------------------------------------------
function pushbutton_compute_Callback(hObject, eventdata, handles)
	com = get(handles.edit_command,'String');
	if (isempty(com))
        oname = get(handles.figure1,'Name');
        set(handles.figure1,'Name','Don''t be IDIOT');    pause(1)
        set(handles.figure1,'Name',oname)
        return;
	end

	if (~isempty(handles.BL))       % We are dealing with Bands arithmetics
        bandArithm(handles, com);   return
	end

	% Those start at 2 because they are meant to be used only when grid names apear repeatedly
	in_g_count = 2;     out_g_count = 2;

	com = move_operator(com);       % Make sure operators are not "glued" to operands
	k = 0;      k = strfind(com,'&');

	try                             % Wrap it here to report any possible error
        if (~isempty(k))            % We have grids
            for (i=1:length(k))     % Loop over grids
                tok = strtok(com(k(i)+1:end));
                if (isempty(handles.name_str))
                    n = [];                 % Here we know that we don't have any pre-loaded grid
                else
                    n = strmatch(tok,handles.name_str);     % n ~= [] when grid is already in memory
                end
            
                if (isempty(n))             % Grid (must be a GMT grid) needs to be loaded
                    load_it = 1;            % Flag to signal that grid must be loaded
                    n_load = strmatch(tok,handles.loaded_grid);
                    if (length(n_load) > 1)         % Grid name comes out more than once
                        n_load = n_load(out_g_count);
                        out_g_count = out_g_count + 1;
                    end
                elseif (length(n) == 1)     % Grid name is not repeated
                    load_it = 0;
                else                        % Grid name comes out more than once
                    n = n(in_g_count);
                    in_g_count = in_g_count + 1;
                    load_it = 0;
                end
            
                if (~load_it)   hand_fig = guidata(handles.h_figs(n));  end
            
                if (i == 1)
                    if (load_it)
                        [tmp.X,tmp.Y,grid_t,tmp.head] = grdread_m([handles.grid_patos{n_load} tok]);
                        grid_t = double(grid_t);      % grdread_m allways outpus singles
                    else
                        grid_t = double(getappdata(hand_fig.figure1,'dem_z'));
                        tmp.X = getappdata(hand_fig.figure1,'dem_x');
                        tmp.Y = getappdata(hand_fig.figure1,'dem_y');
                        tmp.head = getappdata(hand_fig.figure1,'GMThead');
                    end
                    grid.(char(i+96)) = grid_t;
                else
                    if (load_it)
                        [X,Y,grid_t,head] = grdread_m([handles.grid_patos{n_load} tok],'single');
                        grid_t = double(grid_t);      % grdread_m allways outpus singles
                    else
                        grid_t = double(getappdata(hand_fig.figure1,'dem_z'));
                    end
                    grid.(char(i+96)) = grid_t;
                end
            
            end         % Loop over grids
        
            for (i = 1:length(k))
                k = strfind(com,'&');               % We need to recompute '&' positions because they change bellow
                tok = strtok(com(k(i)+1:end));      % Here we get the grid's name
                kf = k(i) + length(tok);            % Find the position of last char of grid's name
                com = [com(1:k(i)) 'grid.' char(i+96) ' ' com(kf+1:end)];
            end
        end
	
        com = strrep(com,'&','');                   % Remove the '&' characters
    
        try             % Try first assuming that we are working from within matlab
            try,        resp = eval(com);                       % See if the Matlab mode worked
            catch,      errordlg(lasterr,'Error');  return      % Shit, it didn't
            end
        catch           % No, we are in standalone mode -- DOESN'T WORK EITHER
            %resp = mexeval(com,['errordlg(lasterr,''Error'')']);
        end
    
        if (length(k) > 0)                  % We had grids in input
            tmp.name = 'Computed grid';
            new_window = mirone(single(resp),tmp);
        elseif (numel(resp) > 1)   % 'resp' is a array. Construct a fake grid
            [m,n] = size(resp);
            resp = single(resp);
            tmp.X = 1:n;        tmp.Y = 1:m;
            [zzz] = grdutils(resp,'-L');  z_min = zzz(1);     z_max = zzz(2);
            tmp.head = [1 n 1 m z_min z_max 0 1 1];
            tmp.name = 'Computed array';
            new_window = mirone(resp,tmp);
        else                % Computations that do not involve grids
            txt = sprintf('%.10f',resp);
            while (txt(end) == '0')     % Remove trailing zeros
                txt(end) = [];
            end
            set(handles.edit_command,'String',txt)
        end
	catch
        errordlg(lasterr,'Error')
	end
    
% ------------------------------------------------------------------------
function bandArithm(handles, com)

	com = move_operator(com);       % Make sure operators are not "glued" to operands
	n = 0;      k = 0;      k = strfind(com,'&');

	try                             % Wrap it here to report any possible error
        if (~isempty(k))            % We have grids
            for (i=1:length(k))     % Loop over bands
                tok = strtok(com(k(i)+1:end));
                n = strmatch(tok,handles.name_str);     % n ~= [] when it is in memory                
                grid.(char(n+96)) = double(handles.BL(:,:,n));
                N(i) = n;
            end

            for (i = 1:length(k))
                tok = strtok(com(k(i)+1:end));      % Here we get the band's name
                kf = k(i) + length(tok);            % Find the position of last char of band's name
                com = [com(1:k(i)) 'grid.' char(N(i)+96) ' ' com(kf+1:end)];
            end
        end

        com = strrep(com,'&','');                   % Remove the '&' characters

        try,        resp = eval(com);                       % See if the Matlab mode worked
        catch,      errordlg(lasterr,'Error');  return      % Shit, it didn't
        end

        if (numel(resp) > 1)   % 'resp' is a array. Construct a fake grid
            if (max(resp(:)) <= 1),      resp = single(resp * 255);      end
            [m,n] = size(resp);
            tmp.X = 1:n;        tmp.Y = 1:m;
            [zzz] = grdutils(resp,'-L');  z_min = zzz(1);     z_max = zzz(2);
            tmp.head = [1 n 1 m z_min z_max 0 1 1];
            tmp.name = 'Computed Band';
            resp = uint8(resp);
            new_window = mirone(resp,tmp);
        else                % Computations that do not involve grids
            txt = sprintf('%.10f',resp);
            while (txt(end) == '0')     % Remove trailing zeros
                txt(end) = [];
            end
            set(handles.edit_command,'String',txt)
        end
	catch
        errordlg(lasterr,'Error')
	end

% ------------------------------------------------------------------------
function str = move_operator(str)
	% Make sure that operatores are not "glued" to operands.
	k = strfind(str,')');
	if (k)  str = strrep(str,')',' ) ');    end
	k = strfind(str,'(');
	if (k)  str = strrep(str,'(',' ( ');    end
	k = strfind(str,'+');
	if (k)  str = strrep(str,'+',' + ');    end
	k = strfind(str,'-');
	if (k)  str = strrep(str,'-',' - ');    end
	k = strfind(str,'*');
	if (k)  str = strrep(str,'*',' * ');    end
	k = strfind(str,'/');
	if (k)  str = strrep(str,'/',' / ');    end
	k = strfind(str,'\');
	if (k)  str = strrep(str,'\',' \ ');    end
	k = strfind(str,'. ');
	if (k)                  % Here we want to have things like '.*' and not '. *'
        while (length(k))
            str = [str(1:k(1)) str(min(k(1)+2,length(str)):end)];
            k = strfind(str,'. ');
        end
	end
	k = strfind(str,' '''); % Hard case of the transpose operator. It may fail often
	if (k)                  % Here we don't want to have things like ")'" turned into ") '"
        while (length(k))
            str = [str(1:k(1)-1) str(min(k(1)+1,length(str)):end)]
            k = strfind(str,' ''');
        end    
	end

% ------------------------------------------------------------------------
function pushbutton_cancel_Callback(hObject, eventdata, handles)
	delete(handles.figure1)

% ------------------------------------------------------------------------
function pushbutton_help_Callback(hObject, eventdata, handles)
	str = sprintf(['This is mainly a grid calculator tool that operates ONLY on GMT grids but\n'...
        'it can also be used in simple scalar calculations.\n\n'...
        'Grids are made available in two ways: (1) all grids loaded in different\n'...
        'Mirone figures show up in the list uppon start of this tool. (2) other\n'...
        'grids may be loaded using the "Load Grid" button.\n\n'...
        'A double click on a grid name moves it to the command zone, where its is\n'...
        'used as any other operand. Complicated expressions can be built, icluding\n'...
        'the use of Matlab commands. Note, however, that there is no error checking\n'...
        'so any comited error will show up as a Matlab error message, which may be\n'...
        'somewhat criptic.\n\n'...
        'You can also generate a grid, for example by entering the rand(200) command.'...
        ' In those cases the generated grid will have coordinates 1:M by 1:N and a grid step of 1']);
	helpdlg(str,'Help')

% --- Creates and returns a handle to the GUI figure. 
function grid_calculator_LayoutFcn(h1,handles);

set(h1,'PaperUnits',get(0,'defaultfigurePaperUnits'),...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'MenuBar','none',...
'Name','Grid calculator',...
'NumberTitle','off',...
'Position',[520 602 669 198],...
'RendererMode','manual',...
'Resize','off',...
'Tag','figure1',...
'UserData',[]);

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@grid_calculator_uicallback,h1,'edit_command_Callback'},...
'HorizontalAlignment','left',...
'Max',3,...
'Position',[10 147 651 41],...
'Style','edit',...
'TooltipString','Enter here any Matlab valid command',...
'Tag','edit_command');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@grid_calculator_uicallback,h1,'listbox_inArrays_Callback'},...
'Position',[10 37 311 91],...
'Style','listbox',...
'TooltipString','Names of currently available arrays',...
'Value',1,'Tag','listbox_inArrays');

uicontrol('Parent',h1,...
'Callback',{@grid_calculator_uicallback,h1,'pushbutton_1_Callback'},...
'FontSize',10,...
'Position',[341 105 23 23],...
'String','1','Tag','pushbutton_1');

uicontrol('Parent',h1,...
'Callback',{@grid_calculator_uicallback,h1,'pushbutton_2_Callback'},...
'FontSize',10,...
'Position',[374 105 23 23],...
'String','2','Tag','pushbutton_2');

uicontrol('Parent',h1,...
'Callback',{@grid_calculator_uicallback,h1,'pushbutton_3_Callback'},...
'FontSize',10,...
'Position',[407 105 23 23],...
'String','3','Tag','pushbutton_3');

uicontrol('Parent',h1,...
'Callback',{@grid_calculator_uicallback,h1,'pushbutton_4_Callback'},...
'FontSize',10,...
'Position',[341 72 23 23],...
'String','4','Tag','pushbutton_4');

uicontrol('Parent',h1,...
'Callback',{@grid_calculator_uicallback,h1,'pushbutton_5_Callback'},...
'FontSize',10,...
'Position',[374 72 23 23],...
'String','5','Tag','pushbutton_5');

uicontrol('Parent',h1,...
'Callback',{@grid_calculator_uicallback,h1,'pushbutton_6_Callback'},...
'FontSize',10,...
'Position',[407 72 23 23],...
'String','6','Tag','pushbutton_6');

uicontrol('Parent',h1,...
'Callback',{@grid_calculator_uicallback,h1,'pushbutton_7_Callback'},...
'FontSize',10,...
'Position',[341 39 23 23],...
'String','7','Tag','pushbutton_7');

uicontrol('Parent',h1,...
'Callback',{@grid_calculator_uicallback,h1,'pushbutton_8_Callback'},...
'FontSize',10,...
'Position',[374 39 23 23],...
'String','8','Tag','pushbutton_8');

uicontrol('Parent',h1,...
'Callback',{@grid_calculator_uicallback,h1,'pushbutton_9_Callback'},...
'FontSize',10,...
'Position',[407 39 23 23],...
'String','9','Tag','pushbutton_9');

uicontrol('Parent',h1,...
'Callback',{@grid_calculator_uicallback,h1,'pushbutton_0_Callback'},...
'FontSize',10,...
'Position',[341 6 23 23],...
'String','0','Tag','pushbutton_0');

uicontrol('Parent',h1,...
'Callback',{@grid_calculator_uicallback,h1,'pushbutton_dot_Callback'},...
'FontSize',10,...
'Position',[374 6 23 23],...
'String','.','Tag','pushbutton_dot');

uicontrol('Parent',h1,...
'Callback',{@grid_calculator_uicallback,h1,'pushbutton_equal_Callback'},...
'FontSize',10,...
'Position',[407 6 23 23],...
'String','=','Tag','pushbutton_equal');

uicontrol('Parent',h1,...
'Callback',{@grid_calculator_uicallback,h1,'pushbutton_devide_Callback'},...
'FontSize',10,...
'Position',[440 105 23 23],...
'String','/','Tag','pushbutton_devide');

uicontrol('Parent',h1,...
'Callback',{@grid_calculator_uicallback,h1,'pushbutton_mull_Callback'},...
'FontSize',10,...
'Position',[440 72 23 23],...
'String','*','Tag','pushbutton_mull');

uicontrol('Parent',h1,...
'Callback',{@grid_calculator_uicallback,h1,'pushbutton_minus_Callback'},...
'FontSize',10,...
'Position',[440 39 23 23],...
'String','-','Tag','pushbutton_minus');

uicontrol('Parent',h1,...
'Callback',{@grid_calculator_uicallback,h1,'pushbutton_plus_Callback'},...
'FontSize',10,...
'Position',[440 6 23 23],...
'String','+','Tag','pushbutton_plus');

uicontrol('Parent',h1,...
'Callback',{@grid_calculator_uicallback,h1,'pushbutton_loadGrid_Callback'},...
'FontSize',10,...
'Position',[11 8 111 24],...
'String','Load Grid',...
'TooltipString','Load ONLY  GMT or Surfer grids',...
'Tag','pushbutton_loadGrid');

uicontrol('Parent',h1,...
'Callback',{@grid_calculator_uicallback,h1,'pushbutton_sin_Callback'},...
'FontSize',10,...
'Position',[491 105 50 23],...
'String','sin','Tag','pushbutton_sin');

uicontrol('Parent',h1,...
'Callback',{@grid_calculator_uicallback,h1,'pushbutton_cos_Callback'},...
'FontSize',10,...
'Position',[551 105 50 23],...
'String','cos','Tag','pushbutton_cos');

uicontrol('Parent',h1,...
'Callback',{@grid_calculator_uicallback,h1,'pushbutton_tan_Callback'},...
'FontSize',10,...
'Position',[611 105 50 23],...
'String','tan','Tag','pushbutton_tan');

uicontrol('Parent',h1,...
'Callback',{@grid_calculator_uicallback,h1,'pushbutton_log10_Callback'},...
'FontSize',10,...
'Position',[491 72 50 23],...
'String','log10','Tag','pushbutton_log10');

uicontrol('Parent',h1,...
'Callback',{@grid_calculator_uicallback,h1,'pushbutton_log_e_Callback'},...
'FontSize',10,...
'Position',[551 72 50 23],...
'String','log e','Tag','pushbutton_log_e');

uicontrol('Parent',h1,...
'Callback',{@grid_calculator_uicallback,h1,'pushbutton_exp_Callback'},...
'FontSize',10,...
'Position',[611 72 50 23],...
'String','exp','Tag','pushbutton_exp');

uicontrol('Parent',h1,...
'Callback',{@grid_calculator_uicallback,h1,'pushbutton_sqrt_Callback'},...
'FontSize',10,...
'Position',[491 39 50 23],...
'String','sqrt','Tag','pushbutton_sqrt');

uicontrol('Parent',h1,...
'Callback',{@grid_calculator_uicallback,h1,'pushbutton_abs_Callback'},...
'FontSize',10,...
'Position',[550 39 50 23],...
'String','abs','Tag','pushbutton_abs');

uicontrol('Parent',h1,...
'Callback',{@grid_calculator_uicallback,h1,'pushbutton_leftPar_Callback'},...
'FontSize',10,...
'Position',[610 39 23 23],...
'String','(','Tag','pushbutton_leftPar');

uicontrol('Parent',h1,...
'Callback',{@grid_calculator_uicallback,h1,'pushbutton_rightPar_Callback'},...
'FontSize',10,...
'Position',[637 39 23 23],...
'String',')','Tag','pushbutton_rightPar');

uicontrol('Parent',h1,...
'Callback',{@grid_calculator_uicallback,h1,'pushbutton_compute_Callback'},...
'FontSize',10,...
'Position',[490 6 71 24],...
'String','Compute','Tag','pushbutton_compute');

uicontrol('Parent',h1,...
'Callback',{@grid_calculator_uicallback,h1,'pushbutton_cancel_Callback'},...
'FontSize',10,...
'Position',[589 6 71 24],...
'String','Cancel','Tag','pushbutton_cancel');

uicontrol('Parent',h1,...
'Callback',{@grid_calculator_uicallback,h1,'pushbutton_help_Callback'},...
'FontSize',10,...
'FontWeight','bold',...
'ForegroundColor',[0 0 1],...
'Position',[221 7 23 23],...
'String','?','Tag','pushbutton_help');

function grid_calculator_uicallback(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
feval(callback_name,hObject,[],guidata(h1));

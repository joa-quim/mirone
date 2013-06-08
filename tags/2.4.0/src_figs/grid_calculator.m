function varargout = grid_calculator(varargin)
% An array calculater 

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

	hObject = figure('Vis','off');
	grid_calculator_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hObject,'right')

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
	handles.home_dir = cd;		% To be able to call put_or_get_file()
	handles.last_dir = cd;		handles.work_dir = cd;
	if (~isempty(h_figs))
        n = 1;
        for (i=1:length(h_figs))
            hand_fig = guidata(h_figs(i));
            % Use a try->catch because ML is too dumb to deal correctly with killed figures            
			try
				Z = getappdata(hand_fig.figure1,'dem_z');
				handles.home_dir = hand_fig.home_dir;		% Since I don't know which is the last good one
				handles.last_dir = hand_fig.last_dir;
				handles.work_dir = hand_fig.work_dir;
			catch
				Z = [];        
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
function edit_command_CB(hObject, handles)
    % Nothing to do here

% ------------------------------------------------------------------------
function listbox_inArrays_CB(hObject, handles)
	% if this is a doubleclick, copy the array's name to the edit box
	if strcmp(get(gcbf,'SelectionType'),'open')
        str = get(hObject,'String');
        sel = get(hObject,'Value');
        com = get(handles.edit_command,'String');
        set(handles.edit_command,'String', [com ' &' str{sel}])
	end

% ------------------------------------------------------------------------
function push_1_CB(hObject, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') '1'])

% ------------------------------------------------------------------------
function push_2_CB(hObject, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') '2'])

% ------------------------------------------------------------------------
function push_3_CB(hObject, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') '3'])

% ------------------------------------------------------------------------
function push_4_CB(hObject, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') '4'])

% ------------------------------------------------------------------------
function push_5_CB(hObject, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') '5'])

% ------------------------------------------------------------------------
function push_6_CB(hObject, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') '6'])

% ------------------------------------------------------------------------
function push_7_CB(hObject, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') '7'])

% ------------------------------------------------------------------------
function push_8_CB(hObject, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') '8'])

% ------------------------------------------------------------------------
function push_9_CB(hObject, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') '9'])

% ------------------------------------------------------------------------
function push_0_CB(hObject, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') '0'])

% ------------------------------------------------------------------------
function push_dot_CB(hObject, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') '.'])

% ------------------------------------------------------------------------
function push_equal_CB(hObject, handles)
	%set(handles.edit_command,'String', [get(handles.edit_command,'String') ' = '])
	push_compute_CB(hObject, handles)

% ------------------------------------------------------------------------
function push_devide_CB(hObject, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') ' / '])

% ------------------------------------------------------------------------
function push_mull_CB(hObject, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') ' * '])

% ------------------------------------------------------------------------
function push_minus_CB(hObject, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') ' - '])

% ------------------------------------------------------------------------
function push_plus_CB(hObject, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') ' + '])

% ------------------------------------------------------------------------
function push_loadGrid_CB(hObject, handles)
	% This function doesn't realy loads the grid. It only stores the grid name.
	% True loading is donne in "Compute"
	str1 = {'*.grd;*.GRD;*.nc;*.NC', 'Grid files (*.grd,*.GRD,*.nc,*.NC)';'*.*', 'All Files (*.*)'};
	[FileName,PathName] = put_or_get_file(handles, str1,'Select GMT grid','get');
	if isequal(FileName,0),		return,		end
	str = get(handles.listbox_inArrays,'String');
	str{end+1} = FileName;
	set(handles.listbox_inArrays,'String',str);
	handles.grid_patos{end+1} = PathName;
	handles.loaded_grid{end+1} = FileName;
	guidata(hObject,handles)

% ------------------------------------------------------------------------
function push_sin_CB(hObject, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') ' sin( '])

% ------------------------------------------------------------------------
function push_cos_CB(hObject, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') ' cos( '])

% ------------------------------------------------------------------------
function push_tan_CB(hObject, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') ' tan( '])

% ------------------------------------------------------------------------
function push_log10_CB(hObject, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') ' log10( '])

% ------------------------------------------------------------------------
function push_log_e_CB(hObject, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') ' log( '])

% ------------------------------------------------------------------------
function push_exp_CB(hObject, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') ' exp( '])

% ------------------------------------------------------------------------
function push_sqrt_CB(hObject, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') ' sqrt( '])

% ------------------------------------------------------------------------
function push_abs_CB(hObject, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') ' abs( '])

% ------------------------------------------------------------------------
function push_leftPar_CB(hObject, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') '( '])

% ------------------------------------------------------------------------
function push_rightPar_CB(hObject, handles)
	set(handles.edit_command,'String', [get(handles.edit_command,'String') ' )'])

% ------------------------------------------------------------------------
function push_compute_CB(hObject, handles)
	com = get(handles.edit_command,'String');
	if (isempty(com))
		oname = get(handles.figure1,'Name');
		set(handles.figure1,'Name','Don''t be FOOL'),		pause(1)
		set(handles.figure1,'Name',oname)
		return
	end

	if (~isempty(handles.BL))       % We are dealing with Bands arithmetics
		bandArithm(handles, com);   return
	end

	% Those start at 2 because they are meant to be used only when grid names apear repeatedly
	in_g_count = 2;     out_g_count = 2;

	%com = move_operator(com);       % Make sure operators are not "glued" to operands (it currently screws names with '-' characters)
	k = strfind(com,'&');

	try                             % Wrap it here to report any possible error
		if (~isempty(k))            % We have grids
			for (i=1:length(k))     % Loop over grids
				tok = strtok(com(k(i)+1:end));
				if (isempty(handles.name_str))
					n = [];					% Here we know that we don't have any pre-loaded grid
				else
					n = strmatch(tok,handles.name_str);     % n ~= [] when grid is already in memory
				end

				if (isempty(n))				% Grid (must be a GMT grid) needs to be loaded
					load_it = 1;			% Flag to signal that grid must be loaded
					n_load = strmatch(tok,handles.loaded_grid);
					if (length(n_load) > 1)         % Grid name comes out more than once
						n_load = n_load(out_g_count);
						out_g_count = out_g_count + 1;
					end
				elseif (numel(n) == 1)		% Grid name is not repeated
					load_it = 0;
				else						% Grid name comes out more than once
					n = n(1);
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
						tmp.head = hand_fig.head;
					end
					grid.(char(i+96)) = grid_t;
				else
					if (load_it)
						[X,Y,grid_t] = grdread_m([handles.grid_patos{n_load} tok],'single');
						grid_t = double(grid_t);      % grdread_m allways outpus singles
					else
						grid_t = double(getappdata(hand_fig.figure1,'dem_z'));
					end
					grid.(char(i+96)) = grid_t;
				end

			end         % Loop over grids

			for (i = 1:length(k))
				k = strfind(com,'&');				% We need to recompute '&' positions because they change bellow
				[tok, r] = strtok(com(k(i)+1:end));	% Here we get the grid's name
				kf = k(i) + numel(tok);				% Find the position of last char of grid's name
				if (strncmp(r,' grid ', 6))			% Patch against those names like "Cropped grid"
					com(kf+1:kf+6) = [];
				end
				com = [com(1:k(i)) 'grid.' char(i+96) ' ' com(kf+1:end)];
			end
		end

		com = strrep(com,'&','');                   % Remove the '&' characters
		com = strrep(com,'e - ','e-');				% base 10 numbers cannot have those spaces
		com = strrep(com,'e + ','e+');

		try             % Try first assuming that we are working from within matlab
			try			resp = eval(com);						% See if the Matlab mode worked
			catch,		errordlg(lasterr,'Error'),	return		% Shit, it didn't
			end
		catch           % No, we are in standalone mode -- DOESN'T WORK EITHER
			%resp = mexeval(com,['errordlg(lasterr,''Error'')']);
		end

		if (~isempty(k))					% We had grids in input
			tmp.head(5:6) = [min(resp(:)) max(resp(:))];
			tmp.name = 'Computed grid';
			mirone(single(resp),tmp);
		elseif (numel(resp) > 1)			% 'resp' is a array. Construct a fake grid
			[m,n] = size(resp);
			resp = single(resp);
			tmp.X = 1:n;        tmp.Y = 1:m;
			[zzz] = grdutils(resp,'-L');  z_min = zzz(1);     z_max = zzz(2);
			tmp.head = [1 n 1 m z_min z_max 0 1 1];
			tmp.name = 'Computed array';
			mirone(resp,tmp);
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
	k = strfind(com,'&');

	try                             % Wrap it here to report any possible error
		if (~isempty(k))            % We have grids
			N = zeros(1, numel(k));
			for (i = 1:numel(k))	% Loop over bands
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

		try			resp = eval(com);                       % See if the Matlab mode worked
		catch,		errordlg(lasterr,'Error');  return      % Shit, it didn't
		end

		if (numel(resp) > 1)   % 'resp' is a array. Construct a fake grid
			if (max(resp(:)) <= 1),      resp = single(resp * 255);      end
			[m,n] = size(resp);
			tmp.X = 1:n;        tmp.Y = 1:m;
			[zzz] = grdutils(resp,'-L');  z_min = zzz(1);     z_max = zzz(2);
			tmp.head = [1 n 1 m z_min z_max 0 1 1];
			tmp.name = 'Computed Band';
			resp = uint8(resp);
			mirone(resp,tmp);
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
	if (k),		str = strrep(str,')',' ) ');    end
	k = strfind(str,'(');
	if (k),		str = strrep(str,'(',' ( ');    end
	k = strfind(str,'+');
	if (k),		str = strrep(str,'+',' + ');    end
	k = strfind(str,'-');
	if (k),		str = strrep(str,'-',' - ');    end
	k = strfind(str,'*');
	if (k),		str = strrep(str,'*',' * ');    end
	k = strfind(str,'/');
	if (k),		str = strrep(str,'/',' / ');    end
	k = strfind(str,'\');
	if (k),		str = strrep(str,'\',' \ ');    end
	k = strfind(str,'. ');
	if (k)					% Here we want to have things like '.*' and not '. *'
        while (~isempty(k))
            str = [str(1:k(1)) str(min(k(1)+2,length(str)):end)];
            k = strfind(str,'. ');
        end
	end
	k = strfind(str,' '''); % Hard case of the transpose operator. It may fail often
	if (k)                  % Here we don't want to have things like ")'" turned into ") '"
        while (~isempty(k))
            str = [str(1:k(1)-1) str(min(k(1)+1,length(str)):end)];
            k = strfind(str,' ''');
        end    
	end

% ------------------------------------------------------------------------
function push_help_CB(hObject, handles)
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
function grid_calculator_LayoutFcn(h1)

set(h1,'PaperUnits',get(0,'defaultfigurePaperUnits'),...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'MenuBar','none',...
'Name','Grid calculator',...
'NumberTitle','off',...
'Position',[520 602 669 198],...
'Resize','off',...
'Tag','figure1');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@grid_calculator_uiCB,...
'HorizontalAlignment','left',...
'Max',3,...
'Position',[10 147 651 41],...
'Style','edit',...
'Tooltip','Enter here any Matlab valid command',...
'Tag','edit_command');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@grid_calculator_uiCB,...
'Position',[10 37 311 91],...
'Style','listbox',...
'Tooltip','Names of currently available arrays',...
'Value',1,'Tag','listbox_inArrays');

uicontrol('Parent',h1,...
'Call',@grid_calculator_uiCB,...
'FontSize',10,...
'Position',[341 105 23 23],...
'String','1','Tag','push_1');

uicontrol('Parent',h1,...
'Call',@grid_calculator_uiCB,...
'FontSize',10,...
'Position',[374 105 23 23],...
'String','2','Tag','push_2');

uicontrol('Parent',h1,...
'Call',@grid_calculator_uiCB,...
'FontSize',10,...
'Position',[407 105 23 23],...
'String','3','Tag','push_3');

uicontrol('Parent',h1,...
'Call',@grid_calculator_uiCB,...
'FontSize',10,...
'Position',[341 72 23 23],...
'String','4','Tag','push_4');

uicontrol('Parent',h1,...
'Call',@grid_calculator_uiCB,...
'FontSize',10,...
'Position',[374 72 23 23],...
'String','5','Tag','push_5');

uicontrol('Parent',h1,...
'Call',@grid_calculator_uiCB,...
'FontSize',10,...
'Position',[407 72 23 23],...
'String','6','Tag','push_6');

uicontrol('Parent',h1,...
'Call',@grid_calculator_uiCB,...
'FontSize',10,...
'Position',[341 39 23 23],...
'String','7','Tag','push_7');

uicontrol('Parent',h1,...
'Call',@grid_calculator_uiCB,...
'FontSize',10,...
'Position',[374 39 23 23],...
'String','8','Tag','push_8');

uicontrol('Parent',h1,...
'Call',@grid_calculator_uiCB,...
'FontSize',10,...
'Position',[407 39 23 23],...
'String','9','Tag','push_9');

uicontrol('Parent',h1,...
'Call',@grid_calculator_uiCB,...
'FontSize',10,...
'Position',[341 6 23 23],...
'String','0','Tag','push_0');

uicontrol('Parent',h1,...
'Call',@grid_calculator_uiCB,...
'FontSize',10,...
'Position',[374 6 23 23],...
'String','.','Tag','push_dot');

uicontrol('Parent',h1,...
'Call',@grid_calculator_uiCB,...
'FontSize',10,...
'Position',[407 6 23 23],...
'String','=','Tag','push_equal');

uicontrol('Parent',h1,...
'Call',@grid_calculator_uiCB,...
'FontSize',10,...
'Position',[440 105 23 23],...
'String','/','Tag','push_devide');

uicontrol('Parent',h1,...
'Call',@grid_calculator_uiCB,...
'FontSize',10,...
'Position',[440 72 23 23],...
'String','*','Tag','push_mull');

uicontrol('Parent',h1,...
'Call',@grid_calculator_uiCB,...
'FontSize',10,...
'Position',[440 39 23 23],...
'String','-','Tag','push_minus');

uicontrol('Parent',h1,...
'Call',@grid_calculator_uiCB,...
'FontSize',10,...
'Position',[440 6 23 23],...
'String','+','Tag','push_plus');

uicontrol('Parent',h1,...
'Call',@grid_calculator_uiCB,...
'FontSize',10,...
'Position',[11 8 111 21],...
'String','Load Grid',...
'TooltipString','Load ONLY  GMT or Surfer grids',...
'Tag','push_loadGrid');

uicontrol('Parent',h1,...
'Call',@grid_calculator_uiCB,...
'FontSize',10,...
'Position',[491 105 50 21],...
'String','sin','Tag','push_sin');

uicontrol('Parent',h1,...
'Call',@grid_calculator_uiCB,...
'FontSize',10,...
'Position',[551 105 50 21],...
'String','cos','Tag','push_cos');

uicontrol('Parent',h1,...
'Call',@grid_calculator_uiCB,...
'FontSize',10,...
'Position',[611 105 50 21],...
'String','tan','Tag','push_tan');

uicontrol('Parent',h1,...
'Call',@grid_calculator_uiCB,...
'FontSize',10,...
'Position',[491 72 50 21],...
'String','log10','Tag','push_log10');

uicontrol('Parent',h1,...
'Call',@grid_calculator_uiCB,...
'FontSize',10,...
'Position',[551 72 50 21],...
'String','log e','Tag','push_log_e');

uicontrol('Parent',h1,...
'Call',@grid_calculator_uiCB,...
'FontSize',10,...
'Position',[611 72 50 21],...
'String','exp','Tag','push_exp');

uicontrol('Parent',h1,...
'Call',@grid_calculator_uiCB,...
'FontSize',10,...
'Position',[491 39 50 21],...
'String','sqrt','Tag','push_sqrt');

uicontrol('Parent',h1,...
'Call',@grid_calculator_uiCB,...
'FontSize',10,...
'Position',[550 39 50 21],...
'String','abs','Tag','push_abs');

uicontrol('Parent',h1,...
'Call',@grid_calculator_uiCB,...
'FontSize',10,...
'Position',[610 39 23 23],...
'String','(','Tag','push_leftPar');

uicontrol('Parent',h1,...
'Call',@grid_calculator_uiCB,...
'FontSize',10,...
'Position',[637 39 23 23],...
'String',')','Tag','push_rightPar');

uicontrol('Parent',h1,...
'Call',@grid_calculator_uiCB,...
'FontSize',10,...
'Position',[589 6 71 21],...
'String','Compute','Tag','push_compute');

uicontrol('Parent',h1,...
'Call',@grid_calculator_uiCB,...
'FontSize',10,...
'FontWeight','bold',...
'ForegroundColor',[0 0 1],...
'Position',[221 7 23 23],...
'String','?','Tag','push_help');

function grid_calculator_uiCB(hObject, eventdata)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));

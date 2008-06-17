function varargout = focal_meca(varargin)
% M-File changed by desGUIDE 
% varargin   command line arguments to focal_meca (see VARARGIN)

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

	if (isempty(varargin))
		errordlg('FOCAL MECA: wrong number of arguments.','Error')
		return
    end

	hObject = figure('Tag','figure1','Visible','off');
	focal_meca_LayoutFcn(hObject);
	handles = guihandles(hObject);
	movegui(hObject,'center');
   
    handles.mirone_fig = varargin{1};
    handMir = guidata(handles.mirone_fig);
	handles.no_file = handMir.no_file;
        
    if (~handMir.no_file && ~handMir.is_projected && ~handMir.geog)
        errordlg('This operation is only possible for geographic data OR when the Map Projection is known','ERROR')
        delete(hObject);    return
    end
    
    handles.is_projected = handMir.is_projected;
    handles.defCoordsIn = handMir.defCoordsIn;
    handles.mironeAxes = handMir.axes1;
    zz = get(handles.mironeAxes,'XLim');
    handles.x_min = zz(1);    handles.x_max = zz(2);
    zz = get(handles.mironeAxes,'YLim');
    handles.y_min = zz(1);    handles.y_max = zz(2);

	handles.date = [];
	
	% Import icons
	f_data = [handMir.home_dir filesep 'data' filesep];
	load([f_data 'mirone_icons.mat'],'Mfopen_ico');
	set(handles.pushbutton_readFile,'CData',Mfopen_ico)
	clear Mfopen_ico;
	
	% Fill the listbox fields with the currently available reading filters
	%str = {'lon,lat,dep,strike,dip,rake,mag,[lon0,lat0,title]'; 'ISF formated catalog (ascii)';};
	str = {'ISF formated catalog (ascii)';
			'Aki & Richard''s convention file';
			'Harvards''s CMT convention file';
			'Harvards''s CMT .ndk file'};
	set(handles.listbox_readFilter,'String',str);
	set(handles.checkbox_plotDate,'Enable','off')

    handles_fake.figure1 = handles.mirone_fig;              % Create a fake handles only for
    handles_fake.axes1 = handles.mironeAxes;                % geog2projected_pts() satisfaction
    handles_fake.geog = handMir.geog;
    handles.handles_fake = handles_fake;
    
    % See what about projection 
    if (handles.is_projected && handles.defCoordsIn > 0)        % We need a proj job here
        tmp = [handles.x_min handles.y_min; handles.x_max handles.y_max];
        lims = [handles.x_min handles.x_max handles.y_min handles.y_max 0];
        [tmp, msg] = geog2projected_pts(handles.handles_fake,tmp, lims);
        x_min = tmp(1,1);           x_max = tmp(2,1);
        y_min = tmp(1,2);           y_max = tmp(2,2);
        handles.lims_geogs = [x_min x_max y_min y_max];     % We'll need this if reading an external file
    end

	% ------------- Give a Pro look (3D) to the frame boxes --------------------
	bgcolor = get(0,'DefaultUicontrolBackgroundColor');
	framecolor = max(min(0.65*bgcolor,[1 1 1]),[0 0 0]);
	h_f = [handles.frame1 handles.frame2];
	for i=1:length(h_f)
		frame_size = get(h_f(i),'Position');
		f_bgc = get(h_f(i),'BackgroundColor');
		usr_d = get(h_f(i),'UserData');
		if abs(f_bgc(1)-bgcolor(1)) > 0.01           % When the frame's background color is not the default's
			frame3D(hObject,frame_size,framecolor,f_bgc,usr_d)
		else
			frame3D(hObject,frame_size,framecolor,'',usr_d)
			delete(h_f(i))
		end
	end
	
	% Add this figure handle to the carra?as list
	plugedWin = getappdata(handMir.figure1,'dependentFigs');
	plugedWin = [plugedWin hObject];
	setappdata(handMir.figure1,'dependentFigs',plugedWin);

	% Choose default command line output for focal_meca
	set(hObject,'Visible','on');
	if (nargout),	varagout{1} = hObject;		end
	guidata(hObject, handles);

% -------------------------------------------------------------------------------------
function listbox_readFilter_Callback(hObject, eventdata, handles)
	switch get(hObject,'Value')
        case 1
			str = sprintf(['Read an ISF formated catalog file (like the ones\n'...
                'you can get from www.isc.ac.uk) and extract\n'...
                'the included (if any) focal mechanisms.']);
        case 2		% Aki & Richards
			str = sprintf(['ASCII file with lon,lat,depth,strike,dip,rake,mag.\n'...
                '8th and 9th columns are optional. If present, they\n'...
                'will determine where the beach ball will be ploted.']);
        case 3		% Simple Harvard's CMT
			str = sprintf(['ASCII file with lon,lat,depth,strike1,dip1,rake1,\n'...
                'strike2,dip2,rake2,mantissa and exponent of moment in N-m.\n'...
                '12th and 13th columns are optional. If present, they\n'...
                'will determine where the beach ball will be ploted.']);
        case 4		% Harvards's CMT .ndk format
			str = sprintf(['Read an CMT .ndk formated catalog file (like the ones\n'...
                'you can get from http://www.globalcmt.org/CMTfiles.html) and extract\n'...
                'the included (if any) focal mechanisms.']);
	end
	set(hObject,'TooltipString',str)

% -------------------------------------------------------------------------------------
function pushbutton_readFile_Callback(hObject, eventdata, handles)
	% OK. Now read the earthquakes_export file and retain only the requested interval
	item = get(handles.listbox_readFilter,'Value');     % Get the reading filter number
	switch item
        case 1      % Read a formated ISF catalog
            str1 = {'*.isf;*.ISF', 'Data files (*.isf,*.ISF)';'*.*', 'All Files (*.*)'};
			filtro = 'isf';
        case 2		% Aki & R
            str1 = {'*.dat;*.DAT', 'Data files (*.dat,*.DAT)';'*.*', 'All Files (*.*)'};
			filtro = 'aki';
        case 3		% simple CMT
            str1 = {'*.dat;*.DAT', 'Data files (*.dat,*.DAT)';'*.*', 'All Files (*.*)'};
			filtro = 'cmt';
        case 4      % Read a file formated with the CMT convention
            str1 = {'*.ndk;*.NDK', 'Data files (*.ndk,*.NDK)';'*.*', 'All Files (*.*)'};
			filtro = 'ndk';
	end

	% Get file name
	[FileName,PathName] = uigetfile(str1,'Select focal file');	pause(0.05)
	if isequal(FileName,0),		return,		end
	fname = [PathName,FileName];
	handles.date = [];			% Allways reset

	try
        set(handles.figure1,'Pointer','watch')
		if (strcmp(filtro,'aki') || strcmp(filtro,'cmt'))      % Aki & Richard or CMT file
            [numeric_data,n_column,error] = read_file(fname);
            if (error),		return,		end
			
			if (~handles.no_file)					% If we know where we are
                if (handles.is_projected && handles.defCoordsIn > 0)        % Image is projected, we need to use this
					x_min = handles.lims_geogs(1);      x_max = handles.lims_geogs(2);
					y_min = handles.lims_geogs(3);      y_max = handles.lims_geogs(4);
                else
					x_min = handles.x_min;      x_max = handles.x_max;
					y_min = handles.y_min;      y_max = handles.y_max;
                end
                
                % Get rid of events that are outside the map limits
				ind = (numeric_data(:,1) < x_min | numeric_data(:,1) > x_max);
                numeric_data(ind,:) = [];
				ind = (numeric_data(:,2) < y_min | numeric_data(:,2) > y_max);
                numeric_data(ind,:) = [];
                if (all(isempty(numeric_data))),	return,		end     % Nothing inside region
				
			else				% If we have a nothing window
				region = [min(numeric_data(:,1)) max(numeric_data(:,1)) min(numeric_data(:,2)) max(numeric_data(:,2)) 1];
				handMir = guidata(handles.mirone_fig);
				mirone('FileNewBgFrame_CB', handMir, region + [-1 1 -1 1 0]*.1);		% Create a background
			end

			if (strcmp(filtro,'aki'))			% Aki & Richard
                if (~(n_column == 7 || n_column == 9 || n_column == 10))
					errordlg('Wrong number of columns for an A&R file','Error');    return
                end
                % [lon lat depth str1 dip1 rake1 mag]
                handles.data = numeric_data(:,1:7);
                mag = numeric_data(:,7);
                switch n_column
                    case 7,			handles.plot_pos = numeric_data(:,1:2);
                    case 9,			handles.plot_pos = numeric_data(:,8:9);
                    case 10,		handles.plot_pos = numeric_data(:,8:9);
				end
			else                % CMT convention
                if (~(n_column == 11 || n_column == 13 || n_column == 14))
                    errordlg('Wrong number of columns for an CMT file','Error');    return
                end
                handles.mantiss_exp = numeric_data(:,10:11);
                mag = (log10(numeric_data(:,10)) + numeric_data(:,11) - 9.1) * 2 / 3;    % In fact Mw
                % [lon lat depth str1 dip1 rake1 str2 dip2 rake2 mag]
                handles.data = [numeric_data(:,1:9) mag];
                switch n_column
					case 11,		handles.plot_pos = numeric_data(:,1:2);
					case 13,		handles.plot_pos = numeric_data(:,12:13);
					case 14,		handles.plot_pos = numeric_data(:,12:13);
				end
			end
	
        elseif (strcmp(filtro,'isf'))				% Read a ISF formated catalog
			opt_R = '-';							% When no image at all
			if (~handles.no_file)					% If we know where we are
                if (handles.is_projected)			% Image is projected, we need this
                    opt_R = sprintf('-R%f/%f/%f/%f',handles.lims_geogs(1:4));
                else
                    opt_R = sprintf('-R%f/%f/%f/%f', handles.x_min, handles.x_max, handles.y_min, handles.y_max);
                end
			end
            [out_d,out_i] = read_isf(fname,opt_R,'-M');
            if (isempty(out_d))		% Nothing inside region
				warndlg('Nope. No mechanisms in this file/region','Warning')
				set(handles.figure1,'Pointer','arrow')
				return
			end

			if (handles.no_file)			% If we have a nothing window
				region = [min(out_d(1,:)) max(out_d(1,:)) min(out_d(2,:)) max(out_d(2,:)) 1];
				handMir = guidata(handles.mirone_fig);
				mirone('FileNewBgFrame_CB', handMir, region + [-1 1 -1 1 0]*.1);		% Create a background
			end

			handles.mantiss_exp = [out_d(10,:)' out_d(11,:)'];
            mag = (log10(out_d(10,:)) + out_d(11,:) - 9.1) * 2 / 3;    % In fact Mw
            handles.data = [out_d(1,:)' out_d(2,:)' out_d(3,:)' out_d(4,:)' out_d(5,:)' out_d(6,:)' out_d(7,:)' ...
							out_d(8,:)' out_d(9,:)' mag'];
            handles.plot_pos = [out_d(1,:)' out_d(2,:)'];		clear out_d;
            n = size(out_i,2);
            handles.date = cell(n,1);
            for (k=1:n)
				handles.date{k} = sprintf('%d/%d/%d',double(out_i(3,k)), double(out_i(2,k)), double(out_i(1,k)));
            end
            clear out_i;
            set(handles.checkbox_plotDate,'Enable','on')
			
        elseif (strcmp(filtro,'ndk'))				% CMT .ndk formated catalog
			[handles.data, handles.mantiss_exp, handles.date, error] = readHarvardCMT(fname);
			if (error),		return,		end
			if (handles.no_file)			% If we have a nothing window
				region = [min(handles.data(:,1)) max(handles.data(:,1)) min(handles.data(:,2)) max(handles.data(:,2)) 1];
				handMir = guidata(handles.mirone_fig);
				mirone('FileNewBgFrame_CB', handMir, region + [-1 1 -1 1 0]*.1);		% Create a background
			end
			handles.plot_pos = handles.data(:,1:2);
			set(handles.checkbox_plotDate,'Enable','on')

		end

		handles.got_userFile = 1;
		handles.usr_DepthMin   = min(handles.data(:,3));    handles.usr_DepthMax = max(handles.data(:,3));
		if (strcmp(filtro,'aki'))		% This is stupid
			handles.usr_MagMin = min(mag);	handles.usr_MagMax = max(mag);
		else
			handles.usr_MagMin = min(handles.data(:,10));	handles.usr_MagMax = max(handles.data(:,10));
		end
		set(handles.edit_MagMin,'String',num2str(floor(handles.usr_MagMin)))
		set(handles.edit_MagMax,'String',num2str(ceil(handles.usr_MagMax)))
		set(handles.edit_DepthMin,'String',num2str(floor(handles.usr_DepthMin)))
		set(handles.edit_DepthMax,'String',num2str(ceil(handles.usr_DepthMax)))
		guidata(hObject,handles)
		set(handles.figure1,'Pointer','arrow')
	catch   % In case of error, set the pointer back to "normal" 
		set(handles.figure1,'Pointer','arrow')
		w{1} = 'An error occured while reading file. Check that it has the apropriate format.';
		w{2} = '';
		w{3} = ['The error message was: ' lasterr];
		w{4} = '';
		w{5} = ['Alternatively, if you are sure that the format is correct check that there ' ...
				'are no empty spaces at the end of your data lines. This may cause an error in decoding the ascii file.'];
		warndlg(w,'Warning')
	end

% -------------------------------------------------------------------------------------
function edit_MagMin_Callback(hObject, eventdata, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx) || xx < 1 || xx > 10)    set(hObject,'String','1');     end

% -------------------------------------------------------------------------------------
function edit_MagMax_Callback(hObject, eventdata, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx) || xx < 1 || xx > 10)    set(hObject,'String','10');    end

% -------------------------------------------------------------------------------------
function edit_Mag5_Callback(hObject, eventdata, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx) || xx < 0)    set(hObject,'String','1');   end

% -------------------------------------------------------------------------------------
function edit_DepthMin_Callback(hObject, eventdata, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx) || xx < 0)    set(hObject,'String','0');   end

% -------------------------------------------------------------------------------------
function edit_DepthMax_Callback(hObject, eventdata, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx) || xx > 900)    set(hObject,'String','900');end

% -------------------------------------------------------------------------------------
function checkbox_depSlices_Callback(hObject, eventdata, handles)
	if (get(hObject,'Value'))
		set(handles.popup_dep0_33,'Enable','on');       set(handles.popup_dep33_70,'Enable','on')
		set(handles.popup_dep70_150,'Enable','on');     set(handles.popup_dep150_300,'Enable','on')
		set(handles.popup_dep300,'Enable','on')
	else
		set(handles.popup_dep0_33,'Enable','off');      set(handles.popup_dep33_70,'Enable','off')
		set(handles.popup_dep70_150,'Enable','off');    set(handles.popup_dep150_300,'Enable','off')
		set(handles.popup_dep300,'Enable','off')
	end

% -------------------------------------------------------------------------------------
function pushbutton_OK_Callback(hObject, eventdata, handles)
	MagMin = str2double(get(handles.edit_MagMin,'String'));
	MagMax = str2double(get(handles.edit_MagMax,'String'));
	DepthMin = str2double(get(handles.edit_DepthMin,'String'));
	DepthMax = str2double(get(handles.edit_DepthMax,'String'));
	item = get(handles.listbox_readFilter,'Value');             % Get the reading filter number
	switch item
        case 1,		filtro = 'isf';
        case 2,		filtro = 'aki';
        case 3,		filtro = 'cmt';
        case 4,		filtro = 'ndk';
	end
	
	if (isnan(MagMin)),			MagMin = 1;         end
	if (isnan(MagMax)),			MagMax = 10;        end
	if (isnan(DepthMin)),		DepthMin = 0;       end
	if (isnan(DepthMax)),		DepthMax = 900;     end

	if (~handles.got_userFile)
        errordlg('Plot What? Your christmas ballons?','Chico Clever');  return;
	end

	% Retain only the requested interval
	ind1 = find(handles.data(:,3) < DepthMin | handles.data(:,3) > DepthMax);
    handles.data(ind1,:) = [];      handles.plot_pos(ind1,:) = [];
    if (~isempty(handles.date))     handles.date(ind1,:) = [];  end
	if (strcmp(filtro,'aki'))
		ind2 = (handles.data(:,7) < MagMin | handles.data(:,7) > MagMax);
	else							% ISF catalog, CMT, CMT .ndk 
        handles.mantiss_exp(ind1,:) = [];   % This risked to heve been left behind
		ind2 = (handles.data(:,10) < MagMin | handles.data(:,10) > MagMax);
        handles.mantiss_exp(ind2,:) = [];
    end
    handles.data(ind2,:) = [];      handles.plot_pos(ind2,:) = [];
    if (~isempty(handles.date))     handles.date(ind2,:) = [];  end

	if (all(isempty(handles.data)))
        warndlg('There were no events left.','Warning');  return;
	end

	if (get(handles.checkbox_depSlices,'Value'))    % We have a depth slice request
		do_depSlices = 1;
		contents = get(handles.popup_dep0_33,'String');     cor_str{1} = contents{get(handles.popup_dep0_33,'Value')};
		contents = get(handles.popup_dep33_70,'String');    cor_str{2} = contents{get(handles.popup_dep33_70,'Value')};
		contents = get(handles.popup_dep70_150,'String');   cor_str{3} = contents{get(handles.popup_dep70_150,'Value')};
		contents = get(handles.popup_dep150_300,'String');  cor_str{4} = contents{get(handles.popup_dep150_300,'Value')};
		contents = get(handles.popup_dep300,'String');      cor_str{5} = contents{get(handles.popup_dep300,'Value')};
	else
		do_depSlices = 0;
	end
    
    % See if we need to project
    if (handles.is_projected && handles.defCoordsIn > 0)        % We need a proj job here
		lims = [handles.x_min handles.x_max handles.y_min handles.y_max];
		[tmp, msg] = geog2projected_pts(handles.handles_fake,handles.data(:,1:2), lims);
		handles.data(:,1:2) = tmp;
		[handles.plot_pos, msg] = geog2projected_pts(handles.handles_fake,handles.plot_pos, lims);
    end

% ------------ OK, now we are ready to plot the mechanisms
oldunit = get(handles.mironeAxes,'Units');		set(handles.mironeAxes,'Units','centimeters')
pos = get(handles.mironeAxes,'Position');		set(handles.mironeAxes,'Units',oldunit)
y_lim = get(handles.mironeAxes,'YLim');
handles.size_fac = (y_lim(2) - y_lim(1)) / (pos(4) - pos(2)) * 0.4;  % Scale facor
Mag5 = get(handles.edit_Mag5,'String');			% Size (cm) of a mag 5 event
handles.Mag5 = str2double(Mag5);
setappdata(handles.mirone_fig,'MecaMag5',Mag5)	% For eventual use in 'write_script'
n_meca = size(handles.data(:,1),1);
axes(handles.mironeAxes)
h_pat = zeros(n_meca,3);
plot_text = get(handles.checkbox_plotDate,'Value');
for (k=1:n_meca)
	if (strcmp(filtro,'aki'))
        [c,d] = patch_meca(handles.data(k,4), handles.data(k,5), handles.data(k,6));
        mag = handles.data(k,7);
	else							% ISF catalog, CMT, CMT .ndk 
        [c,d] = patch_meca(handles.data(k,4), handles.data(k,5), handles.data(k,6), ...
            handles.data(k,7), handles.data(k,8), handles.data(k,9));
        mag = handles.data(k,10);
	end
    dim = handles.size_fac * mag / 5 * handles.Mag5;    % Scale the balls against the selected Mag 5 size
    c = c * dim;    d = d * dim;
    cx = c(:,1) + handles.plot_pos(k,1);
    cy = c(:,2) + handles.plot_pos(k,2);
    dx = d(:,1) + handles.plot_pos(k,1);
    dy = d(:,2) + handles.plot_pos(k,2);
    h_pat(k,3) = line('Parent',handles.mironeAxes,'XData',[handles.data(k,1) handles.plot_pos(k,1)], ...
        'YData',[handles.data(k,2) handles.plot_pos(k,2)], 'Linestyle','-', 'Marker','o', ...
        'MarkerSize',6, 'MarkerFaceColor','k', 'Tag','FocalMecaAnchor');
    if (~do_depSlices)      % Paint all compressive quadrants with black
        h_pat(k,1) = patch(cx,cy, [0 0 0],'Tag','FocalMeca');
    else
        cor = find_color(handles.data(k,3), cor_str);
        h_pat(k,1) = patch(cx,cy, cor,'Tag','FocalMeca');
    end
    h_pat(k,2) = patch(dx,dy, [1 1 1],'Tag','FocalMeca');
	ht = [];
    if (plot_text)          % Plot event text identifier (normaly its date)
        offset = handles.size_fac * mag / 5 * (handles.Mag5 + 0.2);  % text offset regarding the beach ball (2 mm) 
        ht = text(handles.plot_pos(k,1),handles.plot_pos(k,2)+offset,handles.date{k},'HorizontalAlignment', ...
				'Center','VerticalAlignment','Bottom','FontSize',8,'Tag','TextMeca');
        draw_funs(ht,'DrawText');
    end
	
	if (strcmp(filtro,'aki'))
		setappdata(h_pat(k,1),'psmeca_com',[handles.data(k,1:7) handles.plot_pos(k,1:2) ht]);
	else							% ISF catalog, CMT, CMT .ndk 
		setappdata(h_pat(k,1),'psmeca_com',[handles.data(k,1:9) handles.mantiss_exp(k,:) handles.plot_pos(k,1:2) ht]);
    end
    setappdata(h_pat(k,1),'other_hand',[h_pat(k,2) h_pat(k,3) ht]);		% For using in the uiedit
    setappdata(h_pat(k,2),'other_hand',[h_pat(k,1) h_pat(k,3) ht]);		% For using in the uiedit
	
	lim_x = [handles.plot_pos(k,1) handles.plot_pos(k,1) handles.plot_pos(k,1) handles.plot_pos(k,1)] + [-1 -1 1 1]*dim;
	lim_y = [handles.plot_pos(k,2) handles.plot_pos(k,2) handles.plot_pos(k,2) handles.plot_pos(k,2)] + [-1 1 1 -1]*dim;
	setappdata(h_pat(k,1),'Limits',[lim_x(:) lim_y(:)]);            % For using in the uiedit
	setappdata(h_pat(k,2),'Limits',[lim_x(:) lim_y(:)]);            % For using in the uiedit
	set_uicontext(h_pat(k,1));    set_uicontext(h_pat(k,2));
end
hand = guidata(handles.mirone_fig);     % Get the Mirone's handles structure
hand.have_focal = handles.Mag5;         % Signal that we have focal mechanisms and store the Mag5 size symbol
guidata(handles.mirone_fig,hand)        % Save the updated Mirone handles

% -------------------------------------------------------------------------------------
function cor = find_color(z, id)
	if (z < 33),					cor = id{1};
	elseif (z >= 33 && z < 70)		cor = id{2};
	elseif (z >= 70 && z < 150)		cor = id{3};
	elseif (z >= 150 && z < 300)	cor = id{4};
	else							cor = id{5};
	end

% -------------------------------------------------------------------------------------
function [numeric_data,n_column,error] = read_file(fname)
	error = 0;
	hFig = gcf;
	[bin,n_column,multi_seg,n_headers] = guess_file(fname);
	% If msgbox exist we have to move it from behind the main window. So get it's handle
	hMsgFig = gcf;
	if (hFig ~= hMsgFig)        figure(hMsgFig);   end   % If msgbox exists, bring it forward
	% If error in reading file
	if isempty(bin) && isempty(n_column) && isempty(multi_seg) && isempty(n_headers)
		errordlg(['Error reading file ' fname],'Error');
		error = 1;  return
	elseif (bin)			% BINARY
		errordlg('Sorry, reading binary files is not programed','Error');
		error = 1;  return
	end
	
	if (isempty(n_headers))     n_headers = NaN;    end
	if (multi_seg)
		[numeric_data,multi_segs_str,headerlines] = text_read(fname,NaN,n_headers,'>');
	else
		[numeric_data,multi_segs_str,headerlines] = text_read(fname,NaN,n_headers);
	end
	if (hFig ~= hMsgFig);       figure(hFig);   end     % gain access to the drawing figure
	if (iscell(numeric_data))
		n_segments = length(numeric_data);
	else
		n_segments = 1;
	end

% -------------------------------------------------------------------------------------
function set_uicontext(h)
	% Set uicontexts to the Meca patches
	cmenuHand = uicontextmenu;
	set(h, 'UIContextMenu', cmenuHand);
	uimenu(cmenuHand, 'Label', 'Delete this', 'Callback', {@del_Meca,h,'this'});
	uimenu(cmenuHand, 'Label', 'Delete all', 'Callback', {@del_Meca,h,'all'});
	%uimenu(cmenuHand, 'Label', 'Resize', 'Callback', {@resize_Meca,h});
	ui_edit_patch_special(h)

% -------------------------------------------------------------------------------------
function del_Meca(obj,eventdata,h,opt)
	% Delete one or all focal mechanisms
	if (strcmp(opt,'this'))
		delete(getappdata(h,'other_hand'))
		delete(h)
	else
		delete(findobj('Type','patch','Tag','FocalMeca'));
		delete(findobj('Type','line','Tag','FocalMecaAnchor'));
		delete(findobj('Type','text','Tag','TextMeca'));
	end

% -------------------------------------------------------------------------------------
function resize_Meca(obj,eventdata,h)
% 	% Resize the focal mechanisms
% 	handles = guidata(gcf);
% 	h_all = findobj('Type','patch','Tag','FocalMeca');
% 	n_meca = length(h_all) / 2;     % Each ball has two patches
% 	mag = (log10(numeric_data(:,10)) + numeric_data(:,11) - 9.1) * 2 / 3;    % In fact Mw
% 	for (k=1:n_meca)
% 		meca_com = getappdata(h_all(k),'psmeca_com');
% 		dim = handles.size_fac * mag / 5 * handles.Mag5;    % Scale the balls against the selected Mag 5 size
% 		c = c * dim;    d = d * dim;
% 		cx = c(:,1) + handles.plot_pos(k,1);
% 		cy = c(:,2) + handles.plot_pos(k,2);
% 		dx = d(:,1) + handles.plot_pos(k,1);
% 		dy = d(:,2) + handles.plot_pos(k,2);
% 	end

% ----------------------------------------------------------------------------------
function [data, mantiss_exp, eventDate, error] = readHarvardCMT(fname)

	data = [];		mantiss_exp = [];		eventDate = [];		error = 0;
	fid = fopen(fname, 'r');
	if (fid < 0)
		error = 1;		error(['Error opening file ' fname],'Error');
		return
	end

	c = fread(fid,'*char');     fclose(fid);
    todos = strread(c,'%s','delimiter','\n');   clear c fid;
	nEvents = numel(todos)/5;
	
	% Since we are not using info from line 2 & 3 o each event (whicha has a 5 lines descriptor)
	% the best is simple get rid of those unused lines. That even simplifies the parsing
	ind = [(-3 + cumsum(repmat(5,1,nEvents)))' (-2 + cumsum(repmat(5,1,nEvents)))']';
	ind = ind(:);			% [2     3     7     8    12    13    17    18 ...]
	todos(ind) = [];
	
	data = zeros(nEvents, 10);
	mantiss_exp = zeros(nEvents, 2);
	eventDate = cell(nEvents,1);
	
	for (k = 1:nEvents)
		n = (k - 1) * 3 + 1;		% To read the data line numbers correctly
		[data(k,2) data(k,1) data(k,3)] = strread(todos{n}(28:47),'%f %f %f');		% lat lon dep
		eventDate{k} = [todos{n}(15:15) todos{n}(10:13) todos{n}(6:9)];		% day/month/year is the civilized way of displaying dates
		mantiss_exp(k,2) = str2double(todos{n+1}(1:2)) - 7;					% -7 because I want SI units
		[data(k,4) data(k,5) data(k,6) data(k,7) data(k,8) data(k,9)] = strread(todos{n+2}(58:80),'%f %f %f %f %f %f');		% str1 dip1 rake1 str2 dip2 rake2
		mantiss_exp(k,1) = str2double(todos{n+2}(52:56));					% M0 mantissa
		M0  = mantiss_exp(k,1) * 10.0^mantiss_exp(k,2);
		data(k,10) = 2/3 * (log10(M0) - 9.1);			% Mw
	end

% ----------------------------------------------------------------------------------
function figure1_KeyPressFcn(hObject, eventdata)
	if isequal(get(hObject,'CurrentKey'),'escape')
		delete(hObject);
	end

% --- Creates and returns a handle to the GUI figure. 
function focal_meca_LayoutFcn(h1)

set(h1,...
'PaperUnits',get(0,'defaultfigurePaperUnits'),...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'KeyPressFcn',@figure1_KeyPressFcn,...
'MenuBar','none',...
'Name','Focal mechanisms',...
'NumberTitle','off',...
'Position',[520 445 390 355],...
'Resize','off',...
'Tag','figure1');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@focal_meca_uicallback,h1,'listbox_readFilter_Callback'},...
'Position',[50 290 251 61],...
'String',{  'Listbox' },...
'Style','listbox',...
'Value',1,...
'Tag','listbox_readFilter');

uicontrol('Parent',h1,...
'Callback',{@focal_meca_uicallback,h1,'pushbutton_readFile_Callback'},...
'FontWeight','bold',...
'Position',[300 310 23 23],...
'TooltipString','Browse for wanted file',...
'Tag','pushbutton_readFile');

uicontrol('Parent',h1, 'Position',[10 196 371 80], 'Style','frame', 'Tag','frame1');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@focal_meca_uicallback,h1,'edit_MagMin_Callback'},...
'Position',[71 243 47 21],...
'Style','edit',...
'TooltipString','Do not plot events weeker than this',...
'Tag','edit_MagMin');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@focal_meca_uicallback,h1,'edit_MagMax_Callback'},...
'Position',[195 243 47 21],...
'Style','edit',...
'TooltipString','Do not plot events stronger than this',...
'Tag','edit_MagMax');

uicontrol('Parent',h1, 'Position',[18 236 51 30],...
'String',{'Minimum'; 'magnitude' }, 'Style','text');

uicontrol('Parent',h1, 'Position',[143 237 51 30],...
'String',{  'Maximum'; 'magnitude' }, 'Style','text');

uicontrol('Parent',h1,'Position',[10 44 371 141],'Style','frame','Tag','frame2');

uicontrol('Parent',h1, 'Position',[15 140 51 30],...
'String',{'Minimum'; 'depth'}, 'Style','text');

uicontrol('Parent',h1, 'Position',[135 139 51 30],...
'String',{'Maximum'; 'depth' }, 'Style','text');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@focal_meca_uicallback,h1,'edit_Mag5_Callback'},...
'Position',[323 244 47 21],...
'String','0.8',...
'Style','edit',...
'TooltipString','The beach balls will be scaled to this value',...
'Tag','edit_Mag5');

uicontrol('Parent',h1,...
'Position',[69 208 101 15],...
'String','Plot event date',...
'Style','checkbox',...
'TooltipString','Plot time information',...
'Tag','checkbox_plotDate');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@focal_meca_uicallback,h1,'edit_DepthMin_Callback'},...
'Position',[70 145 47 21],...
'String','0',...
'Style','edit',...
'TooltipString','Do not plot events shalower than this',...
'Tag','edit_DepthMin');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@focal_meca_uicallback,h1,'edit_DepthMax_Callback'},...
'Position',[187 145 47 21],...
'Style','edit',...
'TooltipString','Do not plot events deeper than this',...
'Tag','edit_DepthMax');

uicontrol('Parent',h1,...
'Callback',{@focal_meca_uicallback,h1,'checkbox_depSlices_Callback'},...
'Position',[19 110 211 15],...
'String','Use different colors for depth intervals',...
'Style','checkbox',...
'TooltipString','Destinguish the epicenter depths by color',...
'Tag','checkbox_depSlices');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Enable','off',...
'Position',[18 53 62 22],...
'String',{'red'; 'green'; 'blue'; 'cyan'; 'yellow'; 'magenta'; 'kblak' },...
'Style','popupmenu',...
'TooltipString','Symbol color for this depth interval',...
'Value',1,...
'Tag','popup_dep0_33');

uicontrol('Parent',h1,...
'FontSize',10,...
'Position',[21 76 47 16],...
'String','0-33 km',...
'Style','text');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Enable','off',...
'Position',[91 53 62 22],...
'String',{  'red'; 'green'; 'blue'; 'cyan'; 'yellow'; 'magenta'; 'kblak' },...
'Style','popupmenu',...
'TooltipString','Symbol color for this depth interval',...
'Value',2,...
'Tag','popup_dep33_70');

uicontrol('Parent',h1, 'Position',[94 76 54 16],...
'FontSize',10,...
'String','33-70 km',...
'Style','text');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Enable','off',...
'Position',[165 53 62 22],...
'String',{  'red'; 'green'; 'blue'; 'cyan'; 'yellow'; 'magenta'; 'kblak' },...
'Style','popupmenu',...
'TooltipString','Symbol color for this depth interval',...
'Value',3,...
'Tag','popup_dep70_150');

uicontrol('Parent',h1,...
'FontSize',10,...
'Position',[166 76 61 16],...
'String','70-150 km',...
'Style','text');

uicontrol('Parent',h1,...
'FontSize',10,...
'Position',[312 76 55 16],...
'String','> 300 km',...
'Style','text');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Enable','off',...
'Position',[238 53 62 22],...
'String',{  'red'; 'green'; 'blue'; 'cyan'; 'yellow'; 'magenta'; 'kblak' },...
'Style','popupmenu',...
'TooltipString','Symbol color for this depth interval',...
'Value',4,...
'Tag','popup_dep150_300');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Enable','off',...
'Position',[310 53 62 22],...
'String',{'red'; 'green'; 'blue'; 'cyan'; 'yellow'; 'magenta'; 'kblak' },...
'Style','popupmenu',...
'TooltipString','Symbol color for this depth interval',...
'Value',5,...
'Tag','popup_dep300');

uicontrol('Parent',h1,...
'FontSize',10,...
'Position',[235 76 68 16],...
'String','150-300 km',...
'Style','text');

uicontrol('Parent',h1,...
'Callback',{@focal_meca_uicallback,h1,'pushbutton_OK_Callback'},...
'Position',[315 10 66 21],...
'String','OK',...
'Tag','pushbutton_OK');

uicontrol('Parent',h1,...
'Position',[255 238 68 30],...
'String',{  'Magnitude 5'; 'size (cm)' },...
'Style','text');

function focal_meca_uicallback(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
feval(callback_name,hObject,[],guidata(h1));

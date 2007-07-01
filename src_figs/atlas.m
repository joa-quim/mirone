function varargout = atlas(varargin)
% M-File changed by desGUIDE 
% varargin   command line arguments to atlas (see VARARGIN)

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
atlas_LayoutFcn(hObject,handles);
handles = guihandles(hObject);
movegui(hObject,'center');
 
	if (isempty(varargin))
        delete(hObject);        return
	end

    handMir = varargin{1};
    projGMT = getappdata(handMir.figure1,'ProjGMT');
    projWKT = getappdata(handMir.axes1,'ProjWKT');
    if (isempty(projGMT) && isempty(projWKT) && ~handMir.geog)
        errordlg('This operation is only possible for geographic data OR when the Map Projection is known','ERROR')
        delete(hObject);    return
    end

    handles.mirone_fig = handMir.figure1;
    handles.mironeAxes = handMir.axes1;
    handles.is_projected = handMir.is_projected;
    handles.d_path = handMir.path_data;
    handles_fake.figure1 = handles.mirone_fig;              % Create a fake handles only for
    handles_fake.axes1 = handles.mironeAxes;                % geog2projected_pts() satisfaction
    handles.handles_fake = handles_fake;
    handles.h_calling_lims = getappdata(handles.mironeAxes,'ThisImageLims');
    if (isempty(handles.h_calling_lims))
        handles.h_calling_lims = [get(handles.mironeAxes,'Xlim') get(handles.mironeAxes,'Ylim')];
        handles.CeateBG = 1;
    end

	if (isequal(handles.h_calling_lims,[0 1 0 1]))
        handles.CeateBG = 1;
	else
        handles.CeateBG = 0;
	end

	handles.minArea = 0;
	handles.fontSize = 10;
	handles.colors = 1;
	handles.transparency = 0;
	handles.atlas_file = [handles.d_path 'countries_dp5.bin'];
	handles.continents = {'All'};
	handles.atlas = {'All'};
	handles.got_uisetfont = 0;
	
	tmp = dir([handles.d_path filesep 'countries*.bin']);
	if (length(tmp) == 2)
        set(handles.popup_resolution,'String',{'lower' 'higher'})
	else
        set(handles.popup_resolution,'String','lower')
	end

	% It should be a better solution, but for know I have to read from atlas_export.h
	%file = textread([pwd filesep 'mex' filesep 'countries.h'],'%s','delimiter','\n','whitespace','');
	file = dataread('file',[pwd filesep 'mex' filesep 'countries.h'],'%s','delimiter','\n','whitespace','');
	res=findcell('*continent_list',file);
	tmp = cell(9,1);
	tmp {1} = '';
	tmp {2} = 'All';
	m = 2;
	for (k=res.cn+1:res.cn+7)
        tmp{m+1} = file{k}(3:end-2);
        m = m + 1;
	end
	set(handles.listbox_continents,'String',tmp,'Value',2)
	
	% Now fill the atlas_export list box
	res=findcell('*country_list',file);
	n_to_read = length(file) - res.cn;
	tmp = cell(n_to_read,1);
	tmp {1} = '';
	tmp {2} = 'All';
	m = 2;
	for (k=res.cn+1:res.cn+n_to_read-1)
        tmp{m+1} = file{k}(3:end-2);
        m = m + 1;
	end
	set(handles.listbox_allCountries,'String',tmp,'Value',2)
    
	set(hObject,'Visible','on');
	
	% Choose default command line output for atlas_export
	if (nargout),       varargout{1} = hObject;    end
	guidata(hObject, handles);

% ------------------------------------------------------------------------------
function listbox_allCountries_Callback(hObject, eventdata, handles)
	% Hints: contents = get(hObject,'String') returns listbox_allCountries contents as cell array
	%        contents{get(hObject,'Value')} returns selected item from listbox_allCountries
	contents = get(hObject,'String');
	country = contents{get(hObject,'Value')};
	if (~strcmp(country,'All'))     % Set the continents listbox into a "NULL" selection
        set(handles.listbox_continents,'Value',1);
	else                            % Set the continents listbox to "All" as well
        set(handles.listbox_continents,'Value',2);
	end

% ------------------------------------------------------------------------------
function listbox_continents_Callback(hObject, eventdata, handles)
	contents = get(hObject,'String');
	continent = contents{get(hObject,'Value')};
	if (~strcmp(continent,'All'))   % Set the atlas listbox into a "NULL" selection
        set(handles.listbox_allCountries,'Value',1);
	else                            % Set the atlas listbox to "All" as well
        set(handles.listbox_allCountries,'Value',2);
	end

% ------------------------------------------------------------------------------
function popup_resolution_Callback(hObject, eventdata, handles)
	contents = get(hObject,'String');
	if (strcmp(contents{get(hObject,'Value')},'lower'))
        handles.atlas_file = [handles.d_path 'countries_dp5.bin'];
	elseif (strcmp(contents{get(hObject,'Value')},'higher'))
        handles.atlas_file = [handles.d_path 'countries.bin'];
	end
	%handles.atlas_export_file = [handles.d_path contents{get(hObject,'Value')}];
	guidata(hObject,handles)

% ------------------------------------------------------------------------------
function edit_minArea_Callback(hObject, eventdata, handles)
	handles.minArea = str2double(get(hObject,'String'));
	if (isnan(handles.minArea) | handles.minArea < 0)
        set(hObject, 'String','0')
        handles.minArea = 0;
	end
	guidata(hObject,handles)

% ------------------------------------------------------------------------------
function edit_fontSize_Callback(hObject, eventdata, handles)
	handles.fontSize = str2double(get(hObject,'String'));
	if (isnan(handles.fontSize) | handles.fontSize <= 2)
        set(hObject, 'String','10')
        handles.fontSize = 10;
	end
	guidata(hObject,handles)

% ------------------------------------------------------------------------------
function pushbutton_selectFont_Callback(hObject, eventdata, handles)
	handles.Font = uisetfont('Select Font');
	if (~isstruct(handles.Font) & handles.Font == 0)    return;   end
	set(handles.edit_fontSize,'String',num2str(handles.Font.FontSize))
	handles.fontSize = handles.Font.FontSize;
	handles.got_uisetfont = 1;
	guidata(hObject,handles)

% ------------------------------------------------------------------------------
function slider_transparency_Callback(hObject, eventdata, handles)
	handles.transparency = get(hObject,'Value');
	set(handles.text_Transparency,'String',['Transparency = ' num2str(handles.transparency) ' %'])
	guidata(hObject,handles)

% ------------------------------------------------------------------------------
function checkbox_plotNames_Callback(hObject, eventdata, handles)
	if (get(hObject,'Value'))
        set(handles.edit_fontSize,'Enable','on')
        set(handles.pushbutton_selectFont,'Enable','on')
	else
        set(handles.edit_fontSize,'Enable','off')
        set(handles.pushbutton_selectFont,'Enable','off')
	end

% ------------------------------------------------------------------------------
function pushbutton_OK_Callback(hObject, eventdata, handles)
	list_v = get(handles.listbox_continents,'Value');
	list_s = get(handles.listbox_continents,'String');
	continents = list_s(list_v);
	
	list_v = get(handles.listbox_allCountries,'Value');
	list_s = get(handles.listbox_allCountries,'String');
	atlas = list_s(list_v);
	% [h_p, h_t] = my_worldmap([lon(1) lon(2)],[lat(1) lat(2)],'patch');
	
	if (handles.CeateBG)   % We DO NOT have a background to plot. It must be created later
        opt_R = ' ';
	else                         % We have a background map where to plot
        opt_R = ['-R' sprintf('%f/%f/%f/%f',handles.h_calling_lims(1:4))];
	end
	
	if (strcmp(atlas{1},'All') | strcmp(continents{1},'All'))     % Required to plot all atlas (well, the ones who fit in)
            opt_P = ' ';
            opt_T = ' ';
	else
        if (~strcmp(continents{1},'All') & ~strcmp(continents{1},'')) % If a continent was selected
            opt_T = ['-T' continents{1}];
            opt_P = ' ';
        else                                % A country was selected
            opt_P = ['-P' atlas{1}];
            opt_T = ' ';
        end
	end
	paises.ct = country_select(handles.atlas_file,opt_R,opt_P,opt_T,['-A' num2str(handles.minArea)]);

	% Clean up the empty fields in the ct struct (given I could not do it at mex level)
	id = false(length(paises.ct),1);
	for (k = 1:length(paises.ct))
        if (isempty(paises.ct(k).Country)),     id(k) = true;    end
	end
	paises.ct(id) = [];
	
	if (isempty(paises.ct))
        warndlg('There is nothing to plot inside this region','Warning')
        return
	end
	
	handles.transparency = handles.transparency / 100;
	
	% See if user wants country names
	if (get(handles.checkbox_plotNames,'Value'))
        handles.plot_fontSize = handles.fontSize;
	else
        handles.plot_fontSize = [];
	end
	
	% See if user wants uicontexts
	if (get(handles.checkbox_setUicontrols,'Value'))
        handles.uicontrols = 1;
	else
        handles.uicontrols = 0;
	end

    handles.projection = 0;

if (handles.CeateBG)    % Find out the limits off all polygons
    min_x = 1e20;    max_x = -1e20;    min_y = min_x;   max_y = max_x;
    for (k = 1:length(paises.ct))
        min_x0 = min(paises.ct(k).Country(1,:));   
        max_x0 = max(paises.ct(k).Country(1,:));
        min_y0 = min(paises.ct(k).Country(2,:));   
        max_y0 = max(paises.ct(k).Country(2,:));
        min_x = min(min_x0,min_x);          max_x = max(max_x0,max_x);
        min_y = min(min_y0,min_y);          max_y = max(max_y0,max_y);
    end
    region = [min_x max_x min_y max_y];
    lon = [min_x max_x];
    lat = [min_y max_y];
%     % Snap to 0.5 degree increments, with a 0.25 buffer
%     inc = 0.5;
%     buf = 0.25;
%     lat = [inc*floor(min(lat-buf)/inc) inc*ceil(max(lat+buf)/inc)];
%     lon = [inc*floor(min(lon-buf)/inc) inc*ceil(max(lon+buf)/inc)];
    if (handles.projection)
        proj_str = choose_projection(lon,lat);
        min_x = 1e20;    max_x = -1e20;    min_y = min_x;   max_y = max_x;
        for (k = 1:length(paises.ct))
            %tmp = [handles.output.ct(k).Country(1,1:end-1)' handles.output.ct(k).Country(2,1:end-1)'], ...
            %tmp = tmp(1:end-1,:);
            %tmp = mapproject_m(tmp, ...
            tmp = mapproject_m([paises.ct(k).Country(1,1:end-1)' paises.ct(k).Country(2,1:end-1)'], ...
                proj_str{1}, proj_str{2}, '-F', '-C');
            tmp(end+1,1:2) = NaN;       % Reset the last line to NaNs
            % Now we have to find the new limits. Se here we go again
            min_x0 = min(tmp(:,1));         max_x0 = max(tmp(:,1));
            min_y0 = min(tmp(:,2));         max_y0 = max(tmp(:,2));
            min_x = min(min_x0,min_x);      max_x = max(max_x0,max_x);
            min_y = min(min_y0,min_y);      max_y = max(max_y0,max_y);
            paises.ct(k).Country = [tmp(:,1)'; tmp(:,2)'];
        end
        handles.region = [min_x max_x min_y max_y];
    else
        handles.region = [lon lat];
    end
else
    handles.region = [];
end

guidata(hObject, handles);
do_ploting(handles,paises)

% --------------------------------------------------------------------
function do_ploting(handles,paises)
% Plot atlas_export as patches with pre-set colors

huelims = [0 1];    satlims = [.25 .5];     vallims = [1 1];
% rand('state',2);    randomvalues = rand(length(ct),1);   randomhues = huelims(1) + randomvalues*diff(huelims);
% rand('state',2);    randomvalues = rand(length(ct),1);   randomsats = satlims(1) + randomvalues*diff(satlims);
% rand('state',2);    randomvalues = rand(length(ct),1);   randomvals = vallims(1) + randomvalues*diff(vallims);
% rand('state',2);    randomvalues = rand(length(out.ct),1);
% randomhues = huelims(1) + randomvalues*diff(huelims);
% randomsats = satlims(1) + randomvalues*diff(satlims);
% randomvals = vallims(1) + randomvalues*diff(vallims);
% hsv = [randomhues randomsats randomvals];
% pcm = hsv2rgb(hsv);
if (handles.colors == 1)
    pcm = rand(length(paises.ct),3);
elseif (handles.colors == 0)
    rand('state',2);    pcm = rand(length(paises.ct),3);
else    % treta de imitacao de nao cor (MUDAR)
    pcm = repmat([1 1 1], length(paises.ct), 1);
end

if (handles.transparency < 0.01)    no_alfa = 1;
else                                no_alfa = 0;    alfa = handles.transparency; end

% See if we have to create a BG map
if (handles.CeateBG)      % Yes
    region = [handles.region 1];
    if (abs(region(2) - region(1)) > 360 || abs(region(4) - region(3)) > 180),   region(5) = 0;   end
    mirone('FileNewBgFrame_CB',handles.mirone_fig,[],guidata(handles.mirone_fig), region)
end

setappdata(handles.mirone_fig,'AtlasResolution',handles.atlas_file);    % Save this for use in write_gmt_script
    
axes(handles.mironeAxes)       % Make Mirone axes active here

for (k = 1:length(paises.ct))
    id = find(isnan(paises.ct(k).Country(1,:)));
    h = zeros(numel(id),1);
    for (m=1:numel(id))
        if (m == 1)     ini = 1;
        else            ini = id(m-1)+1;    end
        fim = id(m)-1;
        xx = paises.ct(k).Country(1,ini:fim);
        yy = paises.ct(k).Country(2,ini:fim);
        
        if (handles.is_projected)        % We need a proj job here
            [tmp, msg] = geog2projected_pts(handles.handles_fake,[xx; yy]', handles.h_calling_lims);
            xx = tmp(:,1);           yy = tmp(:,2);
        end
        
        if (no_alfa)
            h(m) = patch('Parent',handles.mironeAxes,'XData',xx,'YData', yy,'FaceColor',pcm(k,:), ...
                'Tag','Atlas','UserData',paises.ct(k).Tag);
        else
            h(m) = patch('Parent',handles.mironeAxes,'XData',xx,'YData', yy,'FaceColor',pcm(k,:), ...
                'FaceAlpha',alfa,'Tag','Atlas','UserData',paises.ct(k).Tag);
        end
    end
    if (handles.uicontrols)           % Set patch's uicontextmenu
        draw_funs(h,'country_patch')
    end
end

if (~isempty(handles.plot_fontSize))
	for (k = 1:length(paises.ct))
        str = strrep(paises.ct(k).Tag,'_',' ');
        str(1) = upper(str(1));
        xx = paises.ct(k).Centroide(1);
        yy = paises.ct(k).Centroide(2);
        if (handles.is_projected)        % We need a proj job here
            [tmp, msg] = geog2projected_pts(handles.handles_fake,[xx yy], handles.h_calling_lims);
            xx = tmp(:,1);           yy = tmp(:,2);
        end
        if (~handles.got_uisetfont)
            h = text('Parent',handles.mironeAxes, 'Position',[xx, yy], 'String',str, ...
                'HorizontalAlignment','center', 'FontSize',handles.plot_fontSize,'FontWeight','bold');
        else
            h = text('Parent',handles.mironeAxes, 'Position',[xx, yy], 'String',str, 'HorizontalAlignment','center', ...
                'FontSize',handles.Font.FontSize,'FontWeight',handles.Font.FontWeight,'FontAngle',handles.Font.FontAngle,...
                'FontName',handles.Font.FontName,'FontUnits',handles.Font.FontUnits);
        end
        draw_funs(h,'DrawText')
	end
end

delete(handles.figure1)

% for (k = 1:length(h_p))
%     if (no_alfa)        h = patch('XData',h_p{k,1},'YData', h_p{k,2},'FaceColor',pcm(k,:));
%     else                h = patch('XData',h_p{k,1},'YData', h_p{k,2},'FaceColor',pcm(k,:),'FaceAlpha',alfa);
%     end
%     draw_funs(h,'line_uicontext')      % Set patch's uicontextmenu
% end
% for (k = 1:size(h_t,1))
%     h = text(h_t{k,1},h_t{k,2},h_t{k,3});
%     draw_funs(h,'DrawText')
% end

% [p,name] = fileparts(get(handles.mirone_fig,'Name'));
% k = strfind(name,'_');
% seg = str2num(name(1:k-1));
% out = degree2dms(seg/3600,'DDMM',2,'str');
% str = ['Time  ' out.dd ':' out.mm];
% text(88,-11.5,str,'FontSize',18,'FontWeight','Bold');
% h = text(70,-12.5,'J.Luis','FontSize',12);
% draw_funs(h,'DrawText')

% ------------------------------------------------------------------------------
function proj_str = choose_projection(lon,lat)
opt_R = ['-R' num2str(lon(1)) '/' num2str(lon(2)) '/' num2str(lat(1)) '/' num2str(lat(2))];
if isequal(lat,[-90 90]) & diff(lon) >= 360 % entire globe
    proj_str{1} = '-Rd';     proj_str{2} = '-Jn0/1';     % Robinson
elseif max(abs(lat)) < 20   % straddles equator, but doesn't extend into extreme latitudes
    proj_str{1} = opt_R;     proj_str{1} = '-Jm/1';     % Mercator;
elseif abs(diff(lat)) <= 90 & abs(sum(lat)) > 20  & max(abs(lat)) < 90 % doesn't extend to the pole, not stradling equator
    par1 = lat(2) - mean(lat)/2;
    par2 = lat(1) + mean(lat)/2;
    opt = [num2str(mean(lon)) '/' num2str(mean(lat)) '/' num2str(par2) '/' num2str(par1)];
    proj_str{1} = opt_R;     proj_str{2} = ['-Jd' opt '/1'];     % Equal distace conic
elseif abs(diff(lat)) < 85  & max(abs(lat)) < 90 % doesn't extend to the pole, not stradling equator
    proj_str{1} = opt_R;     proj_str{2} = ['-Ji' num2str(mean(lon)) '/1'];     % Sinusoidal
elseif max(lat) == 90 & min(lat) >= 84
    opt = [num2str(mean(lon)) '/90/1/' num2str(mean(lat))];
    proj_str{1} = opt_R;     proj_str{2} = ['-Js' opt];     % Polar Strographic (N);
elseif min(lat) == -90 & max(lat) <= -80
    opt = [num2str(mean(lon)) '/-90/1/' num2str(mean(lat))];
    proj_str{1} = opt_R;     proj_str{1} = ['-Js' opt];     % Polar Strographic (S);
% elseif max(abs(lat)) == 90 & abs(diff(lon)) < 180
%    proj_str = 'polycon';
elseif max(abs(lat)) == 90 
    opt = [num2str(mean(lon)) '/' num2str(mean(lat)) '/1'];
    proj_str{1} = opt_R;     proj_str{2} = ['-Je' opt];     % Azimuthal Equidistant
else
    proj_str{1} = opt_R;     proj_str{2} = ['-Jj' num2str(mean(lon)) '/1'];     % Miller Cylindrical
end

% ------------------------------------------------------------------------------
function pushbutton_cancel_Callback(hObject, eventdata, handles)
delete(handles.figure1);

% -------------------------------------------------------------------------------
function radiobutton_randColors_Callback(hObject, eventdata, handles)
if (get(hObject,'Value'))
    set(handles.radiobutton_fixedColors,'Value',0)
    set(handles.radiobutton_noColors,'Value',0)
    handles.colors = 1;
elseif (~get(handles.radiobutton_fixedColors) & ~get(handles.radiobutton_noColors))
        set(hObject,'Value',1)
end
guidata(hObject, handles);

% -------------------------------------------------------------------------------
function radiobutton_fixedColors_Callback(hObject, eventdata, handles)
if (get(hObject,'Value'))
    set(handles.radiobutton_randColors,'Value',0)
    set(handles.radiobutton_noColors,'Value',0)
    handles.colors = 0;
elseif (~get(handles.radiobutton_randColors) & ~get(handles.radiobutton_noColors))
        set(hObject,'Value',1)
end
guidata(hObject, handles);

% -------------------------------------------------------------------------------
function radiobutton_noColors_Callback(hObject, eventdata, handles)
if (get(hObject,'Value'))
    set(handles.radiobutton_randColors,'Value',0)
    set(handles.radiobutton_fixedColors,'Value',0)
    handles.colors = -1;
elseif (~get(handles.radiobutton_randColors) & ~get(handles.radiobutton_fixedColors))
        set(hObject,'Value',1)
end
guidata(hObject, handles);

% -------------------------------------------------------------------------------
function checkbox_setUicontrols_Callback(hObject, eventdata, handles)
str = sprintf(['Give the possibility of change colors,\n'...
    'transparency and other attributes. Be awere,\n'...
    'however, that is highly memory consuming.\n'...
    'Particularly with the high definition file.']);
set(hObject,'TooltipString',str)


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
delete(hObject);


% --- Executes on key press over figure1 with no controls selected.%
function figure1_KeyPressFcn(hObject, eventdata, handles)
if isequal(get(hObject,'CurrentKey'),'escape')
    delete(handles.figure1);
end

% --- Creates and returns a handle to the GUI figure. 
function atlas_LayoutFcn(h1,handles);

set(h1,...
'CloseRequestFcn',{@figure1_CloseRequestFcn,handles},...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'KeyPressFcn',{@figure1_KeyPressFcn,handles},...
'MenuBar','none',...
'Name','Atlas',...
'NumberTitle','off',...
'Position',[520 446 411 354],...
'Renderer',get(0,'defaultfigureRenderer'),...
'RendererMode','manual',...
'Resize','off',...
'Tag','figure1',...
'UserData',[]);

h2 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@atlas_uicallback,h1,'listbox_allCountries_Callback'},...
'Max',2,...
'Position',[250 130 151 91],...
'String',{  'Listbox' },...
'Style','listbox',...
'Value',1,...
'Tag','listbox_allCountries');

h3 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@atlas_uicallback,h1,'listbox_continents_Callback'},...
'Max',2,...
'Position',[250 258 151 71],...
'String',{  'Listbox' },...
'Style','listbox',...
'Value',1,...
'Tag','listbox_continents');

h4 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@atlas_uicallback,h1,'popup_resolution_Callback'},...
'Position',[10 308 121 22],...
'String','lower',...
'Style','popupmenu',...
'Value',1,...
'Tag','popup_resolution');

h5 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@atlas_uicallback,h1,'edit_minArea_Callback'},...
'Position',[10 149 47 21],...
'String','0',...
'Style','edit',...
'TooltipString','Polygons with area inferior to this are not drawn',...
'Tag','edit_minArea');

h6 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@atlas_uicallback,h1,'edit_fontSize_Callback'},...
'Enable','off',...
'Position',[10 89 47 21],...
'String','10',...
'Style','edit',...
'TooltipString','Font size in points',...
'Tag','edit_fontSize');

h7 = uicontrol('Parent',h1,...
'Callback',{@atlas_uicallback,h1,'pushbutton_selectFont_Callback'},...
'Enable','off',...
'FontSize',10,...
'FontWeight','bold',...
'Position',[57 87 24 24],...
'TooltipString','Select font annotation of country names',...
'Tag','pushbutton_selectFont');

h8 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@atlas_uicallback,h1,'slider_transparency_Callback'},...
'Max',100,...
'Position',[170 59 231 15],...
'String',{  '' },...
'Style','slider',...
'TooltipString','Use color transparency',...
'Tag','slider_transparency');

h9 = uicontrol('Parent',h1,...
'HorizontalAlignment','left',...
'Position',[170 85 120 15],...
'String','Transparency = 0 %',...
'Style','text',...
'Tag','text_Transparency');

h10 = uicontrol('Parent',h1,...
'Callback',{@atlas_uicallback,h1,'checkbox_plotNames_Callback'},...
'Position',[10 115 110 15],...
'String','Plot country names',...
'Style','checkbox',...
'Tag','checkbox_plotNames');

h11 = uicontrol('Parent',h1,...
'Callback',{@atlas_uicallback,h1,'pushbutton_OK_Callback'},...
'FontSize',10,...
'Position',[220 16 71 24],...
'String','OK',...
'Tag','pushbutton_OK');

h12 = uicontrol('Parent',h1,...
'Callback',{@atlas_uicallback,h1,'pushbutton_cancel_Callback'},...
'FontSize',10,...
'Position',[330 16 71 24],...
'String','Cancel',...
'Tag','pushbutton_cancel');

h13 = uicontrol('Parent',h1,...
'Callback',{@atlas_uicallback,h1,'radiobutton_randColors_Callback'},...
'Position',[10 274 90 15],...
'String','Random colors',...
'Style','radiobutton',...
'Value',1,...
'Tag','radiobutton_randColors');

h14 = uicontrol('Parent',h1,...
'Callback',{@atlas_uicallback,h1,'radiobutton_fixedColors_Callback'},...
'Position',[10 252 84 15],...
'String','Fixed colors',...
'Style','radiobutton',...
'Tag','radiobutton_fixedColors');

h15 = uicontrol('Parent',h1,...
'Callback',{@atlas_uicallback,h1,'radiobutton_noColors_Callback'},...
'Position',[10 228 84 15],...
'String','No color',...
'Style','radiobutton',...
'Tag','radiobutton_noColors');

h16 = uicontrol('Parent',h1,...
'Callback',{@atlas_uicallback,h1,'checkbox_setUicontrols_Callback'},...
'Position',[9 185 155 15],...
'String','Provide controls to polygons',...
'Style','checkbox',...
'Tag','checkbox_setUicontrols');

h17 = uicontrol('Parent',h1,...
'HorizontalAlignment','left',...
'Position',[290 334 71 15],...
'String','By Continents',...
'Style','text',...
'Tag','text2');

h18 = uicontrol('Parent',h1,...
'HorizontalAlignment','left',...
'Position',[290 224 71 15],...
'String','By Countries',...
'Style','text',...
'Tag','text3');

h19 = uicontrol('Parent',h1,...
'HorizontalAlignment','left',...
'Position',[11 334 61 15],...
'String','Resolution',...
'Style','text',...
'Tag','text4');

function atlas_uicallback(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
feval(callback_name,hObject,[],guidata(h1));

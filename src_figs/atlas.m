function varargout = atlas(varargin)
% Helper window to choose countries or continents and plot them

%	Copyright (c) 2004-2015 by J. Luis
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

% $Id$

	if (isempty(varargin))		return,		end

	hObject = figure('Vis','off');
	atlas_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hObject,'center');

    handMir = varargin{1};
    if ( ~handMir.no_file && (~handMir.is_projected && ~handMir.geog) )
		errordlg('This operation is only possible for geographic data OR when the Map Projection is known','ERROR')
		delete(hObject);    return
    end

	handles.mirone_fig = handMir.figure1;
	handles.mironeAxes = handMir.axes1;
	handles.is_projected = handMir.is_projected;
	handles.d_path = handMir.path_data;
	handles_fake.figure1 = handles.mirone_fig;              % Create a fake handles only for
	handles_fake.axes1 = handles.mironeAxes;                % proj2proj_pts() satisfaction
	handles_fake.geog = handMir.geog;
	handles.handles_fake = handles_fake;
	handles.path_tmp = handMir.path_tmp;

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

	% It should find a better solution, but for know I have to read from countries.h
	file = dataread('file',[cd filesep 'mex' filesep 'countries.h'],'%s','delimiter','\n','whitespace','');
	res = findcell('*continent_list',file);
	tmp = cell(10,1);
	tmp {1} = '';
	tmp {2} = 'All';
	m = 2;
	for (k=res.cn+1:res.cn+8)
		tmp{m+1} = file{k}(3:end-2);
		m = m + 1;
	end
	set(handles.listbox_continents,'String',tmp,'Value',2)
	
	% Now fill the atlas_export list box
	res = findcell('*country_list',file);
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
function listbox_allCountries_CB(hObject, handles)
	contents = get(hObject,'String');
	val = get(hObject,'Value');
	country = contents{val(1)};		% This permits multiple selections
	if (~strcmp(country,'All'))		% Set the continents listbox into a "NULL" selection
		set(handles.listbox_continents,'Value',1);
	else							% Set the continents listbox to "All" as well
		set(handles.listbox_continents,'Value',2);
	end

% ------------------------------------------------------------------------------
function listbox_continents_CB(hObject, handles)
	contents = get(hObject,'String');
	val = get(hObject,'Value');
	continent = contents{val(1)};		% This permits multiple selections
	if (~strcmp(continent,'All'))		% Set the atlas listbox into a "NULL" selection
		set(handles.listbox_allCountries,'Value',1);
	else								% Set the atlas listbox to "All" as well
		set(handles.listbox_allCountries,'Value',2);
	end

% ------------------------------------------------------------------------------
function popup_resolution_CB(hObject, handles)
	contents = get(hObject,'String');
	if (strcmp(contents{get(hObject,'Value')},'lower'))
        handles.atlas_file = [handles.d_path 'countries_dp5.bin'];
	elseif (strcmp(contents{get(hObject,'Value')},'higher'))
        handles.atlas_file = [handles.d_path 'countries.bin'];
	end
	%handles.atlas_export_file = [handles.d_path contents{get(hObject,'Value')}];
	guidata(hObject,handles)

% ------------------------------------------------------------------------------
function edit_minArea_CB(hObject, handles)
	handles.minArea = str2double(get(hObject,'String'));
	if (isnan(handles.minArea) || handles.minArea < 0)
        set(hObject, 'String','0')
        handles.minArea = 0;
	end
	guidata(hObject,handles)

% ------------------------------------------------------------------------------
function edit_fontSize_CB(hObject, handles)
	handles.fontSize = str2double(get(hObject,'String'));
	if (isnan(handles.fontSize) || handles.fontSize <= 2)
        set(hObject, 'String','10')
        handles.fontSize = 10;
	end
	guidata(hObject,handles)

% ------------------------------------------------------------------------------
function push_selectFont_CB(hObject, handles)
	handles.Font = uisetfont('Select Font');
	if (~isstruct(handles.Font) && handles.Font == 0),		return,		end
	set(handles.edit_fontSize,'String',num2str(handles.Font.FontSize))
	handles.fontSize = handles.Font.FontSize;
	handles.got_uisetfont = 1;
	guidata(hObject,handles)

% ------------------------------------------------------------------------------
function slider_transparency_CB(hObject, handles)
	handles.transparency = get(hObject,'Value');
	set(handles.text_Transparency,'String',['Transparency = ' num2str(handles.transparency) ' %'])
	guidata(hObject,handles)

% ------------------------------------------------------------------------------
function check_plotNames_CB(hObject, handles)
	if (get(hObject,'Value'))
        set(handles.edit_fontSize,'Enable','on')
        set(handles.push_selectFont,'Enable','on')
	else
        set(handles.edit_fontSize,'Enable','off')
        set(handles.push_selectFont,'Enable','off')
	end

% ------------------------------------------------------------------------------
function push_OK_CB(hObject, handles)
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
        opt_R = sprintf('-R%f/%f/%f/%f',handles.h_calling_lims(1:4));
	end

	if (strcmp(atlas{1},'All') || strcmp(continents{1},'All'))     % Required to plot all atlas (well, the ones who fit in)
            opt_P = ' ';
            opt_T = ' ';
	else
        if (~strcmp(continents{1},'All') && ~strcmp(continents{1},'')) % If a continent was selected
            opt_T = ['-T' continents{1}];
            opt_P = ' ';
        else                                % A country was selected
            opt_P = ['-P' atlas{1}];
            opt_T = ' ';
        end
	end
	
	% See if user wants uicontexts. Test here so that it can be overriden below
	if (get(handles.check_setUicontrols,'Value')),		handles.uicontrols = 1;
	else												handles.uicontrols = 0;
	end

	if (numel(atlas) > 1)
		fname = [handles.path_tmp 'paises.txt'];
		fid = fopen(fname,'w');
		for (k = 1:numel(atlas))
			fprintf(fid,'%s\n', atlas{k});
		end
		fclose(fid);
		opt_P = ['-P' fname];
		paises.ct = country_select(handles.atlas_file,opt_R,opt_P,opt_T,['-A' num2str(handles.minArea)]);
		builtin('delete',fname);
		if ( numel(atlas) < 50 ),	handles.uicontrols = 1;		end			% Override uicontexts choice 
	else
		paises.ct = country_select(handles.atlas_file,opt_R,opt_P,opt_T,['-A' num2str(handles.minArea)]);
		if (opt_T(1) == ' '),	handles.uicontrols = 1;		end				% Override uicontexts choice
	end

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
	if (get(handles.check_plotNames,'Value'))
        handles.plot_fontSize = handles.fontSize;
	else
        handles.plot_fontSize = [];
	end

    handles.projection = 0;

if (handles.CeateBG)    % Find out the limits off all polygons
	min_x = 1e20;    max_x = -1e20;    min_y = min_x;   max_y = max_x;
	for (k = 1:length(paises.ct))
		min_x0 = min(paises.ct(k).Country(1,:));   
		max_x0 = max(paises.ct(k).Country(1,:));
		min_y0 = min(paises.ct(k).Country(2,:));   
		max_y0 = max(paises.ct(k).Country(2,:));
		min_x = min(min_x0,min_x);			max_x = max(max_x0,max_x);
		min_y = min(min_y0,min_y);			max_y = max(max_y0,max_y);
	end
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
            %tmp = c_mapproject(tmp, ...
            tmp = c_mapproject([paises.ct(k).Country(1,1:end-1)' paises.ct(k).Country(2,1:end-1)'], ...
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

% huelims = [0 1];    satlims = [.25 .5];     vallims = [1 1];
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
		rand('state',2);		pcm = rand(length(paises.ct),3);
	else		% treta de imitacao de nao cor (MUDAR)
		pcm = repmat([1 1 1], length(paises.ct), 1);
	end

	if (handles.transparency < 0.01)	no_alfa = 1;
	else								no_alfa = 0;    alfa = handles.transparency;
	end

	% See if we have to create a BG map
	if (handles.CeateBG)      % Yes
		region = [handles.region 1];
		if (abs(region(2) - region(1)) > 360 || abs(region(4) - region(3)) > 180),   region(5) = 0;   end
		if (numel(paises.ct) == 1),		figTitle = paises.ct.Tag;
		else							figTitle = 'Atlas';
		end
		mirone('FileNewBgFrame_CB',guidata(handles.mirone_fig), region, figTitle)
	end

	setappdata(handles.mirone_fig,'AtlasResolution',handles.atlas_file);    % Save this for use in write_gmt_script

	axes(handles.mironeAxes)       % Make Mirone axes active here

for (k = 1:length(paises.ct))
    id = find(isnan(paises.ct(k).Country(1,:)));
    h = zeros(numel(id),1);
    for (m = 1:numel(id))
		if (m == 1),	ini = 1;
		else			ini = id(m-1)+1;
		end
        fim = id(m)-1;
        xx = paises.ct(k).Country(1,ini:fim);
        yy = paises.ct(k).Country(2,ini:fim);
        
        if (handles.is_projected)        % We need a proj job here
			tmp = proj2proj_pts(handles.handles_fake,[xx; yy]', 'srcProj4','+proj=longlat','lim',handles.h_calling_lims);
			xx = tmp(:,1);			yy = tmp(:,2);
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
			tmp = proj2proj_pts(handles.handles_fake,[xx yy], 'srcProj4','+proj=longlat','lim',handles.h_calling_lims);
            xx = tmp(:,1);           yy = tmp(:,2);
        end
        if (~handles.got_uisetfont)
            h = text('Parent',handles.mironeAxes, 'Pos',[xx, yy], 'String',str, ...
                'HorizontalAlignment','center', 'FontSize',handles.plot_fontSize,'FontWeight','bold');
        else
            h = text('Parent',handles.mironeAxes, 'Pos',[xx, yy], 'String',str, 'HorizontalAlignment','center', ...
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
if isequal(lat,[-90 90]) && diff(lon) >= 360 % entire globe
    proj_str{1} = '-Rd';     proj_str{2} = '-Jn0/1';     % Robinson
elseif max(abs(lat)) < 20   % straddles equator, but doesn't extend into extreme latitudes
    proj_str{1} = opt_R;     proj_str{1} = '-Jm/1';     % Mercator;
elseif abs(diff(lat)) <= 90 && abs(sum(lat)) > 20  && max(abs(lat)) < 90 % doesn't extend to the pole, not stradling equator
    par1 = lat(2) - mean(lat)/2;
    par2 = lat(1) + mean(lat)/2;
    opt = [num2str(mean(lon)) '/' num2str(mean(lat)) '/' num2str(par2) '/' num2str(par1)];
    proj_str{1} = opt_R;     proj_str{2} = ['-Jd' opt '/1'];     % Equal distace conic
elseif abs(diff(lat)) < 85 && max(abs(lat)) < 90 % doesn't extend to the pole, not stradling equator
    proj_str{1} = opt_R;     proj_str{2} = ['-Ji' num2str(mean(lon)) '/1'];     % Sinusoidal
elseif max(lat) == 90 && min(lat) >= 84
    opt = [num2str(mean(lon)) '/90/1/' num2str(mean(lat))];
    proj_str{1} = opt_R;     proj_str{2} = ['-Js' opt];     % Polar Strographic (N);
elseif min(lat) == -90 && max(lat) <= -80
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

% -------------------------------------------------------------------------------
function radio_randColors_CB(hObject, handles)
if (get(hObject,'Value'))
    set(handles.radio_fixedColors,'Value',0)
    set(handles.radio_noColors,'Value',0)
    handles.colors = 1;
elseif (~get(handles.radio_fixedColors) && ~get(handles.radio_noColors))
        set(hObject,'Value',1)
end
guidata(hObject, handles);

% -------------------------------------------------------------------------------
function radio_fixedColors_CB(hObject, handles)
if (get(hObject,'Value'))
    set(handles.radio_randColors,'Value',0)
    set(handles.radio_noColors,'Value',0)
    handles.colors = 0;
elseif (~get(handles.radio_randColors) && ~get(handles.radio_noColors))
        set(hObject,'Value',1)
end
guidata(hObject, handles);

% -------------------------------------------------------------------------------
function radio_noColors_CB(hObject, handles)
if (get(hObject,'Value'))
    set(handles.radio_randColors,'Value',0)
    set(handles.radio_fixedColors,'Value',0)
    handles.colors = -1;
elseif (~get(handles.radio_randColors) && ~get(handles.radio_fixedColors))
        set(hObject,'Value',1)
end
guidata(hObject, handles);

% -------------------------------------------------------------------------------
function check_setUicontrols_CB(hObject, handles)
str = sprintf(['Give the possibility of change colors,\n'...
    'transparency and other attributes. Be awere,\n'...
    'however, that is highly memory consuming.\n'...
    'Particularly with the high definition file.']);
set(hObject,'Tooltip',str)


% --- Executes on key press over figure1 with no controls selected.%
function figure1_KeyPressFcn(hObject, eventdata)
	if isequal(get(hObject,'CurrentKey'),'escape')
        delete(hObject);
	end

% -------------------------------------------------------------------------------
function varargout = findcell(k,c,varargin)
% [cixn,cixo] = findcell(key,ca,'-opt')
% cix         = findcell(key,ca,'-opt')
%
%	returns cell number:offset
%	of any occurrence of <key> in the cell array <ca>
%
% key	search pattern
%	format
%		string
%		numeric
%
% ca	a row or col array of cell(s)
%	format
%		string(s)
%		numeric
%	classes must not be mixed
%
% cix	result vector(s)/struct
%	format
%		struct					[def]
%		[cell_nr offset]			[opt:     -p]
%		[cell_nr] [offset]			[nargout:  2]
%
% opt	parameter	action
%  -a :	-		search ACROSS cells
%  -p :	-		output <cix> as vector(s)	[def: struct]
%
% examples
%	ac={'xa';'xxa';'foo';'xaxa';'goo'};
%	findcell('a',ac)
%		% cell nr:offset   1   2
%		% cell nr:offset   2   3
%		% cell nr:offset   4   2
%		% cell nr:offset   4   4
%	nc={[1:10];pi;pi*[1:10]};
%	findcell(pi,nc)
%		% cell nr:offset   2   1
%		% cell nr:offset   3   1
%	findcell([pi pi],nc);
%		% findcell> key not found
%	findcell([pi pi],nc,'-a');
%		% cell nr:offset   2   1
%	zc={['this' 0 'is0'];[0 0 0 'this' 0 'is0']};
%	findcell(['is' 0],zc);
%		% cell nr:offset   1   3
%		% cell nr:offset   2   6
%
% remark
%	program based on a problem by
%	CSSM member <michael robbins>
%	<Vectorizing FINDSTR> (3/11/03)

% search engine
%	since <findstr> does not find <nan>s, we use <nan>s as
%	stop markers after each cell to prevent <across> results,
%	which is the default behavior
%	cellarray
%		cell1 {nan}
%		cell2 {nan}
%		...
%		cellN {nan}
%
%	since <nan>s cannot be converted to <char>s, we must
%	convert cells containing string arrays to double

% created:
%	us	12-Mar-2003
% modified:
%	us	14-Mar-2003 20:40:37	/ TMW

if	nargout
	varargout = cell(nargout,1);
end
if	nargin < 2 || isempty(k) || isempty(c)
    error('Wrong number of input args')
end
[p,k,c] = chk_input(k,c);
if (p.res), return;     end

% get/set options
mat=[];
p = get_opt(p,varargin{:});

% precondition key/cells
c=c(:);
cs=length(c);
cl=length([c{:}]);

if (p.cflg)			% ... do NOT search across <cell>s
	if (ischar(k) && all(p.isc))
		k=double(k);
		inn=cellfun('length',c)';
		c=[c num2cell(0*ones(cs,1))];
		c=reshape(c.',1,2*cs);
		in=cellfun('length',c)';
		c=double([c{:}]);
		c(cumsum(inn+1))=nan;
	elseif	isnumeric(k) && all(p.isn)
		c=[c num2cell(nan*ones(cs,1))];
		c=reshape(c.',1,2*cs);
		in=cellfun('length',c)';
		c=[c{:}];
	else
		fprintf('findcell> unexpected error');
		return
	end
else				% ...	do search across <cell>s
	in=cellfun('length',c)';
	c=[c{:}];
end

% find indices
ix=findstr(k,c);
if	~isempty(ix)
	[mx,mn]=meshgrid(ix,cumsum(in)-in);
	crow=sum(mx>mn,1)';
	if	p.cflg
		crow=.5*(crow-1)+1;
	end
	cix=find(crow==fix(crow));
	ccol=mx-mn;
	ccol(ccol<=0)=nan;
	ccol=min(ccol,[],1)';
	mat=[crow(cix) ccol(cix)];
end

% prepare output
if (isempty(mat))
	%disp('findcell> key not found');
	return
end
if (nargout == 1)
	if	p.pflg
		varargout{1}.par=p;
		varargout{1}.csiz=cl;
		varargout{1}.cn=mat(:,1);
		varargout{1}.co=mat(:,2);
		varargout{1}.cno=mat;
	else
		varargout{1}=mat;
	end
	elseif	nargout == 2
		varargout{1}=mat(:,1);
		varargout{2}=mat(:,2);
	else	% no output requested
		fprintf('cell nr:offset %8d %8d\n',mat.');
end

%--------------------------------------------------------------------------------
function [p,k,c] = chk_input(k,c)
% must do extensive checking ...
p.res=0;
p.key=k;
p.cs=numel(c);
p.cl=length(c);
if	p.cs ~= p.cl
	txt=sprintf('%d/',size(c));
	fprintf('findcell> input must be a ROW or COL array of cells [%s]',txt(1:end-1))
	p.res=1;
	return
end
if	~iscell(c)
	fprintf('findcell> input must be a CELL array [%s]',class(c));
	p.res=2;
	return
end
p.isc=cellfun('isclass',c,'char');
p.isn=cellfun('isclass',c,'double');
if	sum(p.isc) ~= p.cs && sum(p.isn) ~= p.cs
	fprintf('findcell> input contains invalid or mixed classes');
	p.res=3;
	return
end
if	any(isnan(k))
	fprintf('findcell> numeric key must NOT include <NaN>s');
	p.res=4;
	return
end
if	ischar(k) && all(p.isn)
	p.res=5;
	fprintf('findcell> input mismatch: key(char) / cells(double)');
	return
end
if	isnumeric(k) && all(p.isc)
	k=char(k);
end
%--------------------------------------------------------------------------------
function p = get_opt(p,varargin)
	p.cflg=1;	% do NOT search across <cell>s
	p.pflg=1;	% output is a <struct>

	if	nargin
		if	any(strcmp('-a',varargin))
			p.cflg=0;
		end
		if	any(strcmp('-p',varargin))
			p.pflg=0;
		end
	end

% --- Creates and returns a handle to the GUI figure. 
function atlas_LayoutFcn(h1)

set(h1, 'Pos',[520 446 411 354],...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'KeyPressFcn',@figure1_KeyPressFcn,...
'MenuBar','none',...
'Name','Atlas',...
'NumberTitle','off',...
'Resize','off',...
'Tag','figure1');

uicontrol('Parent',h1, 'Pos',[250 130 151 91],...
'BackgroundColor',[1 1 1],...
'Call',@atlas_uiCB,...
'Max',2,...
'Style','listbox',...
'Value',1,...
'Tag','listbox_allCountries');

uicontrol('Parent',h1, 'Pos',[250 258 151 71],...
'BackgroundColor',[1 1 1],...
'Call',@atlas_uiCB,...
'Max',2,...
'Style','listbox',...
'Value',1,...
'Tag','listbox_continents');

uicontrol('Parent',h1, 'Pos',[10 308 121 22],...
'BackgroundColor',[1 1 1],...
'Call',@atlas_uiCB,...
'String','lower',...
'Style','popupmenu',...
'Value',1,...
'Tag','popup_resolution');

uicontrol('Parent',h1, 'Pos',[10 149 47 21],...
'BackgroundColor',[1 1 1],...
'Call',@atlas_uiCB,...
'String','0',...
'Style','edit',...
'Tooltip','Polygons with area inferior to this are not drawn',...
'Tag','edit_minArea');

uicontrol('Parent',h1, 'Pos',[10 89 47 21],...
'BackgroundColor',[1 1 1],...
'Call',@atlas_uiCB,...
'Enable','off',...
'String','10',...
'Style','edit',...
'Tooltip','Font size in points',...
'Tag','edit_fontSize');

uicontrol('Parent',h1, 'Pos',[57 87 24 24],...
'Call',@atlas_uiCB,...
'Enable','off',...
'FontSize',10,...
'FontWeight','bold',...
'Tooltip','Select font annotation of country names',...
'Tag','push_selectFont');

uicontrol('Parent',h1, 'Pos',[170 59 231 15],...
'BackgroundColor',[1 1 1],...
'Call',@atlas_uiCB,...
'Max',100,...
'Style','slider',...
'Tooltip','Use color transparency',...
'Tag','slider_transparency');

uicontrol('Parent',h1, 'Pos',[170 85 120 15],...
'HorizontalAlignment','left',...
'String','Transparency = 0 %',...
'Style','text',...
'Tag','text_Transparency');

uicontrol('Parent',h1, 'Pos',[10 115 130 15],...
'Call',@atlas_uiCB,...
'String','Plot country names',...
'Style','checkbox',...
'Tag','check_plotNames');

uicontrol('Parent',h1, 'Pos',[330 16 71 21],...
'Call',@atlas_uiCB,...
'FontSize',10,...
'String','OK',...
'Tag','push_OK');

uicontrol('Parent',h1, 'Pos',[10 274 110 15],...
'Call',@atlas_uiCB,...
'String','Random colors',...
'Style','radiobutton',...
'Value',1,...
'Tag','radio_randColors');

uicontrol('Parent',h1, 'Pos',[10 252 95 15],...
'Call',@atlas_uiCB,...
'String','Fixed colors',...
'Style','radiobutton',...
'Tag','radio_fixedColors');

uicontrol('Parent',h1, 'Pos',[10 228 84 15],...
'Call',@atlas_uiCB,...
'String','No color',...
'Style','radiobutton',...
'Tag','radio_noColors');

uicontrol('Parent',h1, 'Pos',[9 185 180 15],...
'Call',@atlas_uiCB,...
'String','Provide controls to polygons',...
'Style','checkbox',...
'Tag','check_setUicontrols');

uicontrol('Parent',h1, 'Pos',[290 334 71 15],...
'HorizontalAlignment','left',...
'String','By Continents',...
'Style','text');

uicontrol('Parent',h1, 'Pos',[290 224 71 15],...
'HorizontalAlignment','left',...
'String','By Countries',...
'Style','text');

uicontrol('Parent',h1, 'Pos',[11 334 61 15],...
'HorizontalAlignment','left',...
'String','Resolution',...
'Style','text');

function atlas_uiCB(hObject, eventdata)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));

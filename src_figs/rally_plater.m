function varargout = rally_plater(varargin)
% Display plate motions

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
	rally_plater_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hObject,'center')

	% Import icons
	mir_dirs = getappdata(0,'MIRONE_DIRS');
	if (~isempty(mir_dirs))
		d_path = [mir_dirs.home_dir filesep 'data' filesep];
	else
		d_path = [cd filesep 'data' filesep];
	end
	load([d_path 'mirone_icons.mat'],'zoom_ico','Mfnew_ico','refresh_ico','color_ico','help_ico','Mplay_ico');

	h_toolbar = uitoolbar('parent',hObject,'Clipping', 'on', 'BusyAction','queue','HandleVisibility','on',...
		'Interruptible','on','Tag','FigureToolBar','Visible','on');
	uipushtool('parent',h_toolbar,'Click',@LoadFile_clickedcallback,'cdata',Mfnew_ico,'Tooltip','Load file');
	uitoggletool('parent',h_toolbar,'Click',@zoom_clickedcallback,'Tag','zoom','cdata',zoom_ico,'Tooltip','Zoom');
	uipushtool('parent',h_toolbar,'Click',@newColors_clickedcallback,'cdata',color_ico,...
		'Tooltip','Don''t like these colors. Try others');
	uipushtool('parent',h_toolbar,'Click',@animate_clickedcallback,'cdata',Mplay_ico,'Sep','on','Tooltip','Run Animation');
	uipushtool('parent',h_toolbar,'Click',@reset_clickedcallback,'cdata',refresh_ico,'Tooltip','Reset');
	uipushtool('parent',h_toolbar,'Click',@help_clickedcallback,'Tooltip','Help','cdata',help_ico,'Sep','on');

	handles.path_continent = [pwd filesep 'continents' filesep];
	handles.hCallingFig = [];

	n_bodies = 11;					% This is the number of current default "travelling bodies"
	c_map = rand(n_bodies,3);
	handles.plates_color = c_map;

	% Load default plates
	h_africa = read_plate_bodies(handles,'africa.dat',c_map(1,:),'default','single_seg','Africa');
	h_antartida = read_plate_bodies(handles,'antartica.dat',c_map(2,:),'default','single_seg','Antartica');
	h_arabia = read_plate_bodies(handles,'arabia_dp.dat',c_map(3,:),'default','single_seg','Arabia');
	h_australia = read_plate_bodies(handles,'australia.dat',c_map(4,:),'default','single_seg','Australia');
	h_eurasia = read_plate_bodies(handles,'eu_dp.dat',c_map(5,:),'default','whatever','Eurasia');
	h_greenland = read_plate_bodies(handles,'greenland.dat',c_map(6,:),'default','single_seg','Greenland');
	h_iberia = read_plate_bodies(handles,'iberia.dat',c_map(7,:),'default','single_seg','Iberia');
	h_india = read_plate_bodies(handles,'india_dp.dat',c_map(8,:),'default','single_seg','India');
	h_madagascar = read_plate_bodies(handles,'madagascar.dat',c_map(9,:),'default','single_seg','Madagascar');
	h_Namerica = read_plate_bodies(handles,'na_dp.dat',c_map(10,:),'default','whatever','NorthAmerica');
	h_Samerica = read_plate_bodies(handles,'south_america.dat',c_map(11,:),'default','single_seg','SouthAmerica');

	handles.def_tags = {'Iberia' 'Africa' 'Eurasia' 'NorthAmerica' 'SouthAmerica' 'Antartica' 'Arabia' 'Australia' 'Greenland' 'India'};

	%set(handles.listbox_stages,'String',{'africa2nam_double.stg'; 'eurasia2nam_double.stg'; 'iberia2nam_double.stg'});

	set(handles.listbox_stages,'String',{'AUS_torsvik.stg'; 'EANT_torsvik.stg'; 'EUR_torsvik.stg'; 
										'GRE_torsvik.stg'; 'IBE_Saf.stg'; 'IND_torsvik.stg';
										'MAD_torsvik.stg'; 'NAM_torsvik.stg'; 'SAC_torsvik.stg';
										'africa2nam_double.stg'; 'eurasia2nam_double.stg'; 'iberia2nam_double.stg'});

	str_stgs{1} = [handles.path_continent 'EANT_torsvik.stg'];
	str_stgs{2} = [handles.path_continent 'AUS_torsvik.stg'];
	str_stgs{3} = [handles.path_continent 'EUR_torsvik.stg'];
	str_stgs{4} = [handles.path_continent 'GRE_torsvik.stg'];
	str_stgs{5} = [handles.path_continent 'IBE_Saf.stg'];
	str_stgs{6} = [handles.path_continent 'IND_torsvik.stg'];
	str_stgs{7} = [handles.path_continent 'MAD_torsvik.stg'];
	str_stgs{8} = [handles.path_continent 'NAM_torsvik.stg'];
	str_stgs{9} = [handles.path_continent 'SAC_torsvik.stg'];
	str_stgs{10} = [handles.path_continent 'africa2nam_double.stg'];
	str_stgs{11} = [handles.path_continent 'eurasia2nam_double.stg'];
	str_stgs{12} = [handles.path_continent 'iberia2nam_double.stg'];

	nDefPoles = numel(str_stgs);
	handles.stgs_path = cell(1, nDefPoles);
	for (k = 1:nDefPoles)
		handles.stgs_path{k} = handles.path_continent;
	end

	stg2plate(handles,h_africa,str_stgs,0,1)		% Initialize with no associated poles
	stg2plate(handles,h_antartida,str_stgs,1,0)		% 
	stg2plate(handles,h_arabia,str_stgs,0,1)		% Initialize with no associated poles
	stg2plate(handles,h_australia,str_stgs,2,0)		% 
	stg2plate(handles,h_eurasia,str_stgs,3,0)		% Eurasia
	stg2plate(handles,h_greenland,str_stgs,4,0)		% 
	stg2plate(handles,h_iberia,str_stgs,5,0)		% Iberia
	stg2plate(handles,h_india,str_stgs,6,0)			% 
	stg2plate(handles,h_madagascar,str_stgs,7,0)	% 
	stg2plate(handles,h_Namerica,str_stgs,8,0)		% 
	stg2plate(handles,h_Samerica,str_stgs,9,0)		% 

	handles.plates{1} = h_africa;
	handles.plates{2} = h_antartida;
	handles.plates{3} = h_arabia;
	handles.plates{4} = h_australia;
	handles.plates{5} = h_eurasia;
	handles.plates{6} = h_greenland;
	handles.plates{7} = h_iberia;
	handles.plates{8} = h_india;
	handles.plates{9} = h_madagascar;
	handles.plates{10} = h_Namerica;
	handles.plates{11} = h_Samerica;

	% Make a copy of the default plates and set them invisible
	handles.plates_bak{1} = copyobj(h_africa,handles.axes1);	set(handles.plates_bak{1},'Visible','off');
	handles.plates_bak{2} = copyobj(h_antartida,handles.axes1);	set(handles.plates_bak{2},'Visible','off');
	handles.plates_bak{3} = copyobj(h_arabia,handles.axes1);	set(handles.plates_bak{3},'Visible','off');
	handles.plates_bak{4} = copyobj(h_australia,handles.axes1);	set(handles.plates_bak{4},'Visible','off');
	handles.plates_bak{5} = copyobj(h_eurasia,handles.axes1);	set(handles.plates_bak{5},'Visible','off');
	handles.plates_bak{6} = copyobj(h_greenland,handles.axes1);	set(handles.plates_bak{6},'Visible','off');
	handles.plates_bak{7} = copyobj(h_iberia,handles.axes1);	set(handles.plates_bak{7},'Visible','off');
	handles.plates_bak{8} = copyobj(h_india,handles.axes1);		set(handles.plates_bak{8},'Visible','off');
	handles.plates_bak{9} = copyobj(h_madagascar,handles.axes1);set(handles.plates_bak{9},'Visible','off');
	handles.plates_bak{10} = copyobj(h_Namerica,handles.axes1);	set(handles.plates_bak{10},'Visible','off');
	handles.plates_bak{11} = copyobj(h_Samerica,handles.axes1);	set(handles.plates_bak{11},'Visible','off');
	handles.moved_body = zeros(100,1);

	% Create a circle for use in the Orthographic projection
	handles.h_circ = line(cos(linspace(1,2*pi,180)),sin(linspace(1,2*pi,180)), 'Color', [0 0 0], 'Vis', 'off','LineWidth',2);

	set(handles.edit_ageSlider,'String','200')
	set(handles.slider_age,'Max',200)

	if (~isempty(handles.hCallingFig))						% If we know the handle to the calling fig
		cfig_handles = guidata(handles.hCallingFig);		% get handles of the calling fig    
		handles.last_dir = cfig_handles.last_dir;
		handles.home_dir = cfig_handles.home_dir;
		handles.work_dir = cfig_handles.work_dir;
	elseif (~isempty(mir_dirs))
		handles.home_dir = mir_dirs.home_dir;
		handles.work_dir = mir_dirs.work_dir;
		handles.last_dir = mir_dirs.last_dir;
	else
		handles.home_dir = cd;
		handles.last_dir = cd;
		handles.work_dir = cd;
	end

	%------------ Give a Pro look (3D) to the frame boxes  --------
	new_frame3D(hObject, [handles.txt_Pc handles.txt_Prj])
	%------------- END Pro look (3D) ------------------------------

	guidata(hObject, handles);
	set(hObject,'Visible','on');
	if (nargout),	varargout{1} = hObject;		end

% -----------------------------------------------------------------------------------------
function stg2plate(handles,h_plate,str_stgs,n_pole,opt)
% Set uicontrols, and populate them with the list of available poles contained in STR_STGS
% H_PLATE  -> the handles vector of the plate
% STR_STGS -> cell array of strings with the name (and full path) of the poles file
% N_POLE   -> position on the STR_STGS cell of the plate's default stage poles
% OPT      -> If == 1, Just initialize the uimenus, but do not associate a pole to the plate

	if (nargin == 4)    opt = 0;   end

	cmenuHand = uicontextmenu('Parent',handles.figure1);
	set(h_plate, 'UIContextMenu', cmenuHand);

	uimenu(cmenuHand, 'Label', 'Freeze this plate','Tag','FreezePlate','Call',{@freeze_plate,h_plate});

	if (opt)            % Just initialize
		label = 'Poles --> None set';
	else                % Called with a stage pole to set
		[PATH,FNAME,EXT] = fileparts(str_stgs{n_pole});
		label = ['Poles --> ' FNAME EXT];
	end

	uimenu(cmenuHand, 'Label', label,'Tag','ActiveStg','Separator','on');
	other_stg = uimenu(cmenuHand, 'Label', 'Other Poles', 'Tag', 'OtherPoles');
	for (k = 1:numel(str_stgs))
		[PATH,FNAME,EXT] = fileparts(str_stgs{k});
		name_stg = [FNAME EXT];
		uimenu(other_stg,'Label',name_stg,'Call',{@set_stg,h_plate,str_stgs{k},opt})
	end

	% Set the delete uicontext
	tag = get(h_plate,'Tag');
	if (iscell(tag))    tag = tag{1};   end
	uimenu(cmenuHand,'Label',['Delete this ' tag ' element'],'Tag','DeleteSingle','Sep','on',...
		'Call',{@delete_element,h_plate,0});
	uimenu(cmenuHand,'Label',['Delete all ' tag ' family'],'Tag','DeleteAll','Call',{@delete_element,h_plate,1});

	% set the "change color" uicontext
	uimenu(cmenuHand,'Label',['Change ' tag ' color'],'Separator','on','Call',@change_color);

	if (opt)            % Signal that this plate has no associated poles
		set(h_plate,'UserData',0)
	else                % Set the default (to this plate) stage poles
		set(h_plate,'UserData',str_stgs{n_pole})
	end

% -----------------------------------------------------------------------------------------
function set_stg(obj,eventdata,patchHand,poles_file,opt)
% Associate the stage pole POLES_FILE (it includes its full path name) to the plate PATCHHAND
%
% Now the problems arrive when a new element was added lately (relatively to time of the
% 'patchHand' creation). In this case the "old" patchHand vector was not updated. That is
% why we are going to fish all handles that share the same 'Tag' as the one of old patchHand

	handles = guidata(obj);
	this_body_tag = get(patchHand,'Tag');
	if (iscell(this_body_tag))		this_body_tag = this_body_tag{1};	end
	h = findobj(handles.figure1,'Tag',this_body_tag);
	patchHand = h;

	for (k = 1:numel(patchHand))
		set(patchHand(k),'UserData',poles_file)
	end
	[PATH,FNAME,EXT] = fileparts(poles_file);
	% Update the 'Active pole' label
	h = get(patchHand,'UIContextMenu');
	if (iscell(h))
		h = cat(1,h{:});
	end
	hA = findobj(h,'Tag','ActiveStg');
	set(hA,'Label',['Poles --> ' FNAME EXT])

	if (opt)        % Plate was initialized without any associated poles. So do it now
		set(patchHand,'UserData',poles_file)
	end

% -----------------------------------------------------------------------------------------
function freeze_plate(obj,eventdata,patchHand)
% 
pole_file = get(patchHand,'UserData');
if (isequal(pole_file,0))               % This plate has currently no associated poles
    try                                 % See if we have a backup
        pole_file = getappdata(patchHand,'BackupPoles');
    catch                               % No, we do not. So there is nothing to freeze
        return
    end
else
    for (k=1:length(patchHand))
        setappdata(patchHand(k),'BackupPoles',pole_file)
    end
end

if (iscell(pole_file))      pole_file = pole_file{1};   end

if (strcmp(get(obj,'Checked'),'on'))
    set(obj,'Checked','off')
    %set(patchHand,'UserData',pole_file);% Note that, in this case, pole_file was retrieved from appdata
else
    set(obj,'Checked','on')
    pole_file = 0;        % This is what realy freezes the plate
end

for (k=1:length(patchHand))
    set(patchHand(k),'UserData',pole_file)
end

% -----------------------------------------------------------------------------------------
function delete_element(obj,eventdata,h_element,mode)
% MODE = 0 Delete only the selected family element
% MODE = 1 Delete all family elements

handles = guidata(obj);
h_current = gco;

% The killing is easy, but we have also to clean all traces of the victim
tag = get(h_current,'Tag');
if (iscell(tag))    tag = tag{1};   end

% Again as above, problems arrive when a new element was added lately (relatively to time of the
% 'h_element' creation). In this case the "old" h_element vector was not updated. That is
% why we going to fish all handles that share the same 'Tag' as the one of old h_element
h = findobj(handles.figure1,'Tag',tag,'Visible','on');
h_element = sort(h);
h_element = h_element(:)';		% Make sure it is a sorted row vector

n_pos = 0;
for (k = 1:numel(handles.plates))
    if (isequal(h_element,handles.plates{k}))
        n_pos = k;  % OK, we found the position on the variable list
        break
    end
end

if (n_pos == 0)     % Souldn't occur. Semething wrong happened before
    warndlg('Sorry, due to a previous unknown error I cannot delete the element(s)','Warning')
    return
end

id = strmatch(tag,handles.def_tags);

if (mode)           % Delete all family
    handles.plates(n_pos) = [];
    handles.plates_bak(n_pos) = [];
    handles.def_tags(id) = [];
    handles.plates_color(id,:) = [];
    delete(h_element);
else                % Delete the element. But, and if it is the only one in the family?
    for (k=1:length(handles.plates{n_pos}))
        if (isequal(h_current,handles.plates{n_pos}(k)))
            handles.plates{n_pos}(k) = [];
            handles.plates_bak{n_pos}(k) = [];
            delete(h_current);
            break
        end
    end
    % We still have to check if it was the only one in the family.
    % The stupid thing is that I found no way of searching for an empty cell other than this
    for (l=1:numel(handles.plates))
        if (isempty(handles.plates{l}))     % It was
            handles.plates(l) = [];
            handles.plates_bak(l) = [];
            handles.def_tags(id) = [];
            handles.plates_color(id,:) = [];
            break
        end
    end
end
 
guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------------------
function newColors_clickedcallback(obj,eventdata)
% Randomly change colors of all families

handles = guidata(obj);
for (k=1:length(handles.plates))
    tag = get(handles.plates{k},'Tag');
    if (iscell(tag))    tag = tag{1};   end
    h_patch = findobj(handles.figure1,'Tag',tag,'Type','patch','Visible','on');
    h_line = findobj(handles.figure1,'Tag',tag,'Type','line','Visible','on');
    c = rand(1,3);
    if (~isempty(h_patch))
        set(h_patch,'FaceColor',c);
    end
    if (~isempty(h_line))
        set(h_line,'Color',c);
    end
    handles.plates_color(k,1:3) = c;    % Update the color list
end
refresh;

% -----------------------------------------------------------------------------------------
function change_color(obj,eventdata)
% Change the color of all family
handles = guidata(obj);
tag = get(gco,'Tag');
if (iscell(tag))    tag = tag{1};   end
h_patch = findobj(handles.figure1,'Tag',tag,'Type','patch','Visible','on');
h_line = findobj(handles.figure1,'Tag',tag,'Type','line','Visible','on');

c = uisetcolor;
if length(c) > 1            % That is, if a color was selected
    if (~isempty(h_patch))
        set(h_patch,'FaceColor',c,'EdgeColor',[0 0 0]);
    end
    if (~isempty(h_line))
        set(h_line,'Color',c);
    end
    refresh;
end

% -----------------------------------------------------------------------------------------
function [h_plate,out_tag] = read_plate_bodies(handles,plate_body,c_map,opt1,opt2,opt3)
% OPT1 == 'default' or ...
% OPT2 == 'single_seg' or 'whatever'
% OPT3 =  Tag of patch/line

	out_tag = [];				% To be used if an external file is read and has an identification in 1st line
	if (nargin == 1)			% Read an external file
		fname = handles;
		opt1 = 'whatever';
	end

	if (strcmp(opt1,'default'))
		if (strcmp(opt2,'single_seg'))
			plate = text_read([handles.path_continent plate_body]);
			h_plate = patch(plate(:,1),plate(:,2),c_map,'FaceAlpha',0.5,'Tag',opt3);
		else		% Multiseg plate
			plate = text_read([handles.path_continent plate_body],NaN,NaN,'>');
			n_segments = length(plate);
			h_plate = zeros(1, n_segments);
			for (k = 1:n_segments)
				h_plate(k) = patch(plate{k}(:,1),plate{k}(:,2),c_map,'FaceAlpha',0.5,'Tag',opt3);
			end
		end
	else
		[bin,n_column,multi_seg,n_headers] = guess_file(fname);
		if isempty(bin) && isempty(n_column) && isempty(multi_seg) && isempty(n_headers)
			errordlg(['Error reading file ' fname],'Error');    return
		end
		if (n_column < 2)
			errordlg('File error. Your file doesn''t have at least 2 columns','Error'); return
		end
		if (multi_seg)
			[h_plate,multi_segs_str,headerlines,hdr_txt] = text_read(fname,NaN,NaN,'>');
		else
			[h_plate,multi_segs_str,headerlines,hdr_txt] = text_read(fname,NaN,NaN);
		end
		if (~isempty(hdr_txt))
			[t,r]=strtok(hdr_txt{1}(2:end));
			if (isempty(r))         % Only one word. We interpret this as the file tag
				out_tag = t;
			end
		end
	end

% -----------------------------------------------------------------------------------------
function [poles,p_name] = read_stgs(poles_file)
% Read an stage poles file and store it in a cell array
% When the first line in file has the form # Name   (no blanks in 'Name')
% it is assumed that 'Name' is the moving plate's name and the value will
% be stored in P_NAME. Though I do not use it yet, I have some ideas for future uses

fid = fopen(poles_file,'r');
c = fread(fid,'*char').';
fclose(fid);
s=strread(c,'%s','delimiter','\n');
ix = strmatch('#',s);

hdr = s(ix);
n_hdr = length(hdr);
n_stgs = length(s)-n_hdr;
poles = zeros(n_stgs,5);
try
	for i=1:n_stgs
         tmp = sscanf(s{i+n_hdr}','%f',5);
         poles(i,1:5) = tmp';
	end
catch
    errordlg(['The file ' poles_file 'is not a properly formated Stage poles file.'],'Error');
    poles = [];     p_name = '';
    return
end

p_name = '';
if (~isempty(ix))           % There are header lines in the file
    [t,r]=strtok(hdr{1}(2:end));
    if (isempty(r))         % Only one word. We interpret this as the plate's name
        p_name = t;
    end
end

% --------------------------------------------------------------------------------------------------
function zoom_clickedcallback(obj,eventdata)
	if (strcmp(get(obj,'State'),'on'))
		zoom_j('on')
	else
		zoom_j('off')
	end

% --------------------------------------------------------------------
function animate_clickedcallback(obj, eventdata)
% Animate
handles = guidata(obj);
do_proj = get(handles.radio_projOrtho,'Val');		% Are we working with projected coords?
orig = [get(handles.slider_projOrigLat,'Val') get(handles.slider_projOrigLon,'Val') ...
        get(handles.slider_projOrigPitch,'Val')];	% In case we are, we will need this
opt_L = ['-L' get(handles.edit_ageStart,'String') '/' get(handles.edit_ageStop,'String') ...
        '/' get(handles.edit_ageStep,'String')];
n_flow = [];

nPlates = numel(handles.plates);
out = cell(1, nPlates);		n_data = zeros(1, nPlates);		n_seg = zeros(1, nPlates);

for (k = 1:nPlates)					% Loop over number of plates with associated poles
    pole_file = get(handles.plates{k},'UserData');
    if (isequal(pole_file,0))				% This plate has no associated poles
        continue
    end
    x = get(handles.plates_bak{k},'XData');
    y = get(handles.plates_bak{k},'YData');
    if (iscell(pole_file))					% This occurs with plates made of more than one element
        pole_file = pole_file{1};			% In this case, they are all equal
        if (isequal(pole_file,0))	continue,	end
        nseg = length(x);
        xx = [];    yy = [];
        for (l = 1:nseg)					% Loop over number of segments of this active plate
            xx = [xx; NaN; x{l}];
            yy = [yy; NaN; y{l}];
        end
        x = xx;     y = yy;     clear xx yy;
    end
    opt_E = ['-E' pole_file];
    [out{k},n_data(k),n_seg(k),n_flow] = telha_m([x y], opt_E, opt_L);

    if (do_proj)							% Project coords
        [x,y] = orthographic(out{k}(:,1), out{k}(:,2), orig);
        out{k} = [x y];
    else
        [out{k}(:,2),out{k}(:,1)] = trimpatch(out{k}(:,2), [-Inf 89], noJumpLong(out{k}(:,1)), [-180 180]);
    end
    
    handles.moved_body(k) = k;				% Keep track of which body was moved
end
if (isempty(n_flow))	return,		end		% There wasn't any active plate
clear x y;

% Now do the animation

frm_step = str2double(get(handles.edit_frameInterval,'String'));    % Get frame interval
t_step = str2double(get(handles.edit_ageStep,'String'));

% Get the animation direction. Forward (from past to present) or Backward (the oposit)
if (get(handles.radio_animForward,'Value'))
    j_dir = n_flow:-1:1;		% Forward
else
    j_dir = 1:n_flow;			% Backward
end

id = cell(1, nPlates);		ini = cell(1, nPlates);		fim = cell(1, nPlates);

for (j = j_dir)									% For each time increment
	for (m = 1:numel(out))						% Loop over the number of moved plates
        if (isempty(out{m}))	continue,	end % Freezed body
        [x{m},y{m}] = get_time_slice(out{m},n_data(m),n_seg(m),j);
        x{m}(length(x{m})+1) = NaN;     y{m}(length(y{m})+1) = NaN;    % Needed for processing multiple patches.
        id{m} = find(isnan(x{m}));
	
		for (i = 1:length(id{m}) )				% Cycle through and display each element
			if (i == 1)		ini{m}(i) = 1;
			else			ini{m}(i) = id{m}(i-1)+1;
			end
			fim{m}(i) = id{m}(i)-1;
			try									% This is crutial when working on proj coords
				set(handles.plates{m}(i),'XData',x{m}(ini{m}(i):fim{m}(i)),'YData',y{m}(ini{m}(i):fim{m}(i)))
			end
		end
	end
	pause(frm_step);
	set(handles.figure1,'Name',['Rally Plater    ' num2str((j-1)*t_step) '  Ma'])
end
guidata(handles.figure1,handles);

% --------------------------------------------------------------------
function [x,y] = get_time_slice(data,n_data,n_seg,n,first)
	i1 = (n-1)*(n_data + n_seg) + 2;
	i2 = i1 + n_data - 1 + n_seg - 1;
	x = data(i1:i2,1);		y = data(i1:i2,2);

% --------------------------------------------------------------------
function reset_clickedcallback(hObject, eventdata)
handles = guidata(hObject);
set(handles.figure1,'Name','Rally Plater')
set(handles.edit_ageSlider,'String','0')
set(handles.slider_age,'Value',0)

if (get(handles.radio_projOrtho,'Value'))
    swap_proj(handles)                          % Just reset to the initial Ortho projected coords and return
    return
end

n = find(handles.moved_body);					% Get the number of moved bodies (it is a vector)
if (isempty(n))		return,		end				% Nothing moved yet (the guy is plaing with the buttons)
for (k = 1:length(n))							% Loop over moved bodies
    for (m = 1:length(handles.plates{n(k)}))	% Loop over each element of the outer loop moved body
        set(handles.plates{n(k)}(m),'XData',get(handles.plates_bak{n(k)}(m),'XData'), ...
            'YData',get(handles.plates_bak{n(k)}(m),'YData'));
    end
end
handles.moved_body = zeros(100,1);				% Reset to bodies_moved = 0
guidata(handles.figure1,handles);

% --------------------------------------------------------------------------------------------------
function help_clickedcallback(obj,eventdata)
str = sprintf(['This tool is a mix between a Plate reconstruction demo and.\n'...
    'something you can use to do real work.\n\n'...
    'There are several available tools on the web that do plate reconstructions\n'...
    'but this one has a big advantage. It is ready to work when you start it,\n'...
    'and is stupidly easy to use\n\n'...
    'I wont go for the trouble to explain what the buttons do. They self\n'...
    'explain with their tooltips. To access other options right-click over the "plates".\n\n'...
    'Another matter is to use your own data. There are two types of data: the\n'...
    'poles files, which control how the plates move; and the plates files.\n'...
    'The poles used in this program are Stage Poles (not finite rotation\n'...
    'poles). They have a particular format described in the "rotconverter" GMT\n'...
    'program, but you can look at them at the "continents" directory of\n'...
    'Mirone''s installation. To create new stage poles files from finite\n'...
    'rotations poles hit the "Make stage poles" and follow instructions there.\n\n'...
    'The plates files are simple ascii files with two columns: LON & LAT\n'...
    'Now, those files may be multisegments, with character ">" serving as\n'...
    'separator. The individual segments may be closed (equal first and last\n'...
    'point). If they are closed, the contour is filled, otherwise they are\n'...
    'drawn as lines. Use this for drawing your invention of the COB.\n\n'...
    'Very important: The first line in file must be of the form ''# Name''\n'...
    '(without the quotes). This is how I know how to associate a plate and\n'...
    'its companion lines with a stage pole.\n'...
    'You can see the default plates name when you right-click it, so if you\n'...
    'want to add it another element (e.g. the 3000 m bathymetric) you must\n'...
    'give that file the exact same name as the plate you want it to join.\n']);
helpdlg(str,'Help')

% --------------------------------------------------------------------
function listbox_stages_CB(hObject, handles)
% Hints: contents = get(hObject,'String') returns listbox_stages contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_stages

% --------------------------------------------------------------------
function push_unSet_CB(hObject, handles)
% Unlink all plates to their associated poles
	for (k = 1:numel(handles.plates))
		set(handles.plates{k},'UserData',0)
	end

% --------------------------------------------------------------------
function push_setABS_CB(hObject, handles)
% Load the "absolute" poles
	str_stgs{1} = [handles.path_continent 'SAfr_ABS_hybrid.stg'];
	str_stgs{2} = [handles.path_continent 'EANT_ABS_hybrid.stg'];
	str_stgs{3} = 0;			% Arabia. We don't currently have poles for it
	str_stgs{4} = [handles.path_continent 'AUS_ABS_hybrid.stg'];
	str_stgs{5} = [handles.path_continent 'EUR_ABS_hybrid.stg'];
	str_stgs{6} = [handles.path_continent 'GRE_ABS_hybrid.stg'];
	str_stgs{7} = [handles.path_continent 'IBE_ABS_hybrid.stg'];
	str_stgs{8} = [handles.path_continent 'IND_ABS_hybrid.stg'];
	str_stgs{9} = [handles.path_continent 'MAD_ABS_hybrid.stg'];
	str_stgs{10} = [handles.path_continent 'NAM_ABS_hybrid.stg'];
	str_stgs{11} = [handles.path_continent 'SAC_ABS_hybrid.stg'];
	
	for (k = 1:numel(handles.plates))
		set(handles.plates{k},'UserData',str_stgs{k})
	end
	
	str = {'SAfr_ABS_hybrid.stg'; 'EANT_ABS_hybrid.stg'; 'AUS_ABS_hybrid.stg'; 'EUR_ABS_hybrid.stg';
		'GRE_ABS_hybrid.stg'; 'IBE_ABS_hybrid.stg'; 'IND_ABS_hybrid.stg'; 'MAD_ABS_hybrid.stg';
		'NAM_ABS_hybrid.stg'; 'SAC_ABS_hybrid.stg'};
	set(handles.listbox_stages, 'Str', str)

% --------------------------------------------------------------------
function push_loadStages_CB(hObject, handles)
% Get poles file name

    [FileName, PathName, handles] = put_or_get_file(handles,{ ...
			'*.stg;*.dat;*.DAT', 'Stage poles (*.stg,*.dat,*.DAT)';'*.*', 'All Files (*.*)'},'Select poles file','get');
    if isequal(FileName,0),		return,		end
	fname = [PathName FileName];

	% Check that it is a goodly formated stage poles file. We don't want surprises in mexs
	dumb = read_stgs(fname);
	if (isempty(dumb))		return,		end			% Bad stage poles file

	max_age = max(dumb(1,3), get(handles.slider_age,'Max'));
	set(handles.slider_age,'Max',max_age)

	str_stgs = get(handles.listbox_stages,'String');
	if (~iscell(str_stgs))		str_stgs = {str_stgs; FileName};
	else						str_stgs{end+1} = FileName;
	end

	[str_stgs,id] = sort(str_stgs);					% We want them sorted
	handles.stgs_path{end+1} = PathName;
	handles.stgs_path = handles.stgs_path(id);		% Uppdate also the poles path variable

	set(handles.listbox_stages,'String',str_stgs)

	hA = findobj(handles.figure1,'Tag','OtherPoles');
	for (k = 1:numel(hA))					% Loop over uimenus of activated plates
		hC = sort(get(hA(k),'Children'));
		for (m = 1:numel(hC))				% Loop over uimenus offering alternative stage poles
			cb = get(hC(m),'Call');
			cb{3} = [handles.stgs_path{m} str_stgs{m}];
			set(hC(m),'Label',str_stgs{m},'Call',cb)
		end
		% But we still need to add another uimenu to account for the newly imported stage poles
		uimenu(hA(k),'Label',str_stgs{end},'Call',{@set_stg,cb{2},[handles.stgs_path{end} str_stgs{end}],1})
	end

	guidata(handles.figure1, handles)

% --------------------------------------------------------------------
function push_makeStages_CB(hObject, handles)
% 
	fid = fopen([handles.path_continent 'lista_polos.dat'],'rt');
	c = fread(fid,'*char').';
	fclose(fid);
	s = strread(c,'%s','delimiter','\n');

	choosebox('Name','One Euler list',...
				'PromptString','List of poles:',...
				'SelectString','Selected poles:',...
				'ListSize',[380 300],'ListString',s);

% --------------------------------------------------------------------
function LoadFile_clickedcallback(obj, eventdata)
% Get the external file and draw it. All closed polygons are drawn as patches

	handles = guidata(obj);
    [FileName,PathName] = put_or_get_file(handles,{ ...
			'*.dat;*.DAT', 'Data file (*.dat,*.DAT)';'*.*', 'All Files (*.*)'},'Select input xy file name','get');
    if isequal(FileName,0),		return,		end

	[xy,tag] = read_plate_bodies([PathName FileName]);
	if (isempty(tag))
		tag = 'unknown';        % New body/plate
	end
	
	if (iscell(xy))
		for (k=1:length(xy))
			draw_element(handles,xy{k}(:,1),xy{k}(:,2),tag)
			handles = guidata(handles.figure1);     % The handles was updated inside draw_element, but this loop
		end                                         % does not know it. So we have to update also here.
	else
		draw_element(handles,xy(:,1),xy(:,2),tag)
	end

% --------------------------------------------------------------------
function draw_element(handles,x,y,tag)
% This is what I expect to happen. If the object shares an existing tag, it will be added
% to the family of those handles that share the tag.
% If it doesn't, a new family is created. This means that all untaged objects, which were
% given the tag 'unknown', will be added to that family. As a consequence, ...

	dd = [x(1) y(1)] - [x(end) y(end)];
	if (~all(dd))                   % Closed polygon
		h = patch('XData',x,'YData',y,'FaceColor','none','EdgeColor','k','FaceAlpha',0.5,'Tag',tag);
		prop_str = 'EdgeColor';
		is_patch = 1;
	else                            % A polyline (open)
		h = line(x,y,'Tag',tag,'Color',[0 0 0]);
		prop_str = 'Color';
		is_patch = 0;
	end

	if (is_patch)
		uistack(h,'bottom')			% Needed in order to not hide previous elements, but dangerous
	end

	id = strmatch(tag,handles.def_tags);
	if (~isempty(id))				% The tag matches one existent body. Add its handle to that body list
		set(h,prop_str,handles.plates_color(id,:));     % Give it the family color
		stg = get(handles.plates{id}(1),'UserData');
		uictx = get(handles.plates{id}(1),'Uicontext');
		set(h,'UserData',stg,'Uicontext',uictx);
		handles.plates{id}(end+1) = h;			% Store the new handle in its family cell
		handles.plates_bak{id}(end+1) = copyobj(h,handles.axes1);
		set(handles.plates_bak{id}(end),'Visible','off')
	else						% We have a new element with a new Tag - It means, a new plate or a stray element
		c_map = rand(1,3);						% Create a new color for this new family
		set(h,prop_str,c_map)					% Set it to the new color
		str_stgs = get(handles.listbox_stages,'String');
		if (~iscell(str_stgs))      str_stgs = {str_stgs};      end
		stg2plate(handles,h,str_stgs,0,1)		% Initialize with no associated poles
		handles.plates{end+1} = h;				% Increase the plates counter
		handles.plates_bak{end+1} = copyobj(h,handles.axes1);   % Make a copy of the new object
		set(handles.plates_bak{end},'Visible','off')
		handles.def_tags{end+1} = tag;			% Update the tag list
		handles.plates_color(end+1,1:3) = c_map;% Update the color list
	end

	guidata(handles.figure1,handles)

% --------------------------------------------------------------------
function slider_age_CB(hObject, handles)
% Reconstruct to a particular age determined by the slider value

do_proj = get(handles.radio_projOrtho,'Value');
orig = [get(handles.slider_projOrigLat,'Value') get(handles.slider_projOrigLon,'Value') ...
        get(handles.slider_projOrigPitch,'Value')];

n_flow = [];
age = get(hObject,'Value');
set(handles.edit_ageSlider,'String',num2str(age));

nPlates = numel(handles.plates);
out = cell(1, nPlates);		n_data = zeros(1, nPlates);		n_seg = zeros(1, nPlates);

for (k = 1:nPlates)				% Loop over number of plates with associated poles
    pole_file = get(handles.plates{k},'UserData');
    if (isequal(pole_file,0))				% This plate has currently no associated poles
        continue
    end
    x = get(handles.plates_bak{k},'XData');
    y = get(handles.plates_bak{k},'YData');
    if (iscell(pole_file))
        pole_file = pole_file{1};			% In this case, they are all equal
        if (isequal(pole_file,0))		continue,	end
        nseg = length(x);
        xx = [];    yy = [];
        for (l = 1:nseg)					% Loop over number of segments of this active plate
            xx = [xx; NaN; x{l}];
            yy = [yy; NaN; y{l}];
        end
    else
        xx = [NaN; x];
        yy = [NaN; y];
    end
    x = xx;     y = yy;     clear xx yy;
    
    opt_E = ['-E' pole_file];
    [out{k},n_data(k),n_seg(k),n_flow] = telha_m([x y], age, opt_E, '-P');
    if (isempty(out{k}))	continue,	end % This hapens when the one of the stage poles reach its oldest age
    
    if (do_proj)							% Project coords
        [x,y] = orthographic(out{k}(:,1), out{k}(:,2), orig);
        out{k} = [x y];
    else
        [out{k}(:,2),out{k}(:,1)] = trimpatch(out{k}(:,2), [-Inf 89], noJumpLong(out{k}(:,1)), [-180 180]);
    end
    
    if (n_seg(k) == 1)						% This is because telha outputs an extra line with (0,0)
        out{k}(end,:) = [];					% in the case of a single segment with one rotation only
    end
    handles.moved_body(k) = k;				% Keep track of which body was moved
end
if (isempty(n_flow))	return,		end		% There wasn't any active plate
clear x y;

% Plot the plates at their new positions
id = cell(1, nPlates);		ini = cell(1, nPlates);		fim = cell(1, nPlates);

for (m = 1:length(out))						% Loop over the number of moved plates
    if (isempty(out{m}))	continue,	end % Freezed body
    if (n_seg(m) > 1)
        out{m}(1,:) = [];					% Mata a primeira linha de anoes
    end
    x{m} = out{m}(:,1);             y{m} = out{m}(:,2);
    x{m}(length(x{m})+1) = NaN;     y{m}(length(y{m})+1) = NaN;    % Needed for processing multiple patches.
    id{m} = find(isnan(x{m}));

	for (i = 1:length(id{m}) )				% Cycle through and display each element
        if (i == 1)     ini{m}(i) = 1;
        else            ini{m}(i) = id{m}(i-1)+1;    end
        fim{m}(i) = id{m}(i)-1;
        try									% We realy nead this when working on proj coords
            set(handles.plates{m}(i),'XData',x{m}(ini{m}(i):fim{m}(i)),'YData',y{m}(ini{m}(i):fim{m}(i)));
        end
	end
end
guidata(handles.figure1,handles)

% --------------------------------------------------------------------
function lon = noJumpLong(lon)
% Finds jumps in longitude (crossing of date line) and removes them
% By using this function before calling trimpatch, the result is that
% The plates will desapear at < -180 or > 180. It seams resonable.
dif = diff(lon);
id = find(abs(dif) > 181);

while ~isempty(id)
	lon(id+1:end) = lon(id+1:end)-sign(dif(id(1)))*360;
	dif = diff(lon);
	id = find(abs(dif) > 181);
end	

% --------------------------------------------------------------------
function edit_ageSlider_CB(hObject, handles)
val = str2double(get(hObject,'String'));
if (val > get(handles.slider_age,'Max'))
    set(handles.slider_age,'Max',val)
end
set(handles.slider_age,'Value',val)

% Call the slider callback to do the rest of the work
slider_age_CB(handles.slider_age, [], handles)

% --------------------------------------------------------------------
function [x,y] = orthographic(lon, lat, origin)
% Project geographic coords into spherical orthographic projection

D2R = pi / 180;
lat = lat * D2R;   %  Convert to radians
lon = lon * D2R;
origin = origin * D2R;

% Rotate sphere only if it is needed
if (sum(origin))
    [lon,lat] = rotate(lon,lat,origin);
end
% Compute the azims and distances on the sphere
orig0 = zeros(size(lat));
azim  = azimuth_geo(orig0, orig0, lat, lon,'radians');
rng = acos(cos(lat).*cos(lon));         % NOTE: the formula simplifies because we are computing from a [0 0] origin

%  Trim data exceeding the visible part of the sphere
[rng,azim] = trimpatch(rng, [-Inf 89]*D2R, azim, [-inf inf]);

x = sin(rng) .* sin(azim);
y = sin(rng) .* cos(azim);

%--------------------------------------------------------------------------------------------------
function [lon1,lat1] = rotate(lon,lat,orig)
%ROTATE  Rotate data for specified orig and orientation (angles are in radians)

rot1 = [cos(orig(2)) sin(orig(2))  0        % Rotation matrix about x axis
       -sin(orig(2)) cos(orig(2))  0
	    0            0             1];
rot2 = [cos(orig(1)) 0 sin(orig(1))         % Rotation matrix about y axis
        0            1 0
	   -sin(orig(1)) 0 cos(orig(1))];
rot3 = [1  0            0
        0  cos(orig(3)) sin(orig(3))        % Rotation matrix about z axis
        0 -sin(orig(3)) cos(orig(3))];

rot = rot3 * rot2 * rot1;                   % Euler rotation matrix

%  Move pi/2 points epsilon inward to prevent round-off problems with pi/2 points.
epsilon = 1e-6;
indx = find(abs(pi/2 - abs(lat)) <= epsilon);
if ~isempty(indx)
	lat(indx) = (pi/2 - epsilon) * sign(lat(indx));
end

%  Prevent possible confusion with points at +180 or -180 degrees
lon = atan2(sin(lon*(1 - 1e-6)),cos(lon*(1 - 1e-6)));

%  Compute the new x,y,z point in cartesian space
xyz = ( rot * ([cos(lat).*cos(lon) cos(lat).*sin(lon) sin(lat)]') )';% We want column vectors

% epsilon = 1.0e-8;
% indx = find(abs(xyz(:,1)) <= epsilon & abs(xyz(:,2)) <= epsilon);   % Be careful with x & y nearely 0 in atan2
% if ~isempty(indx);   x(indx) = 0;  y(indx) = 0;   end

[lon1, lat1] = cart2sph(xyz(:,1),xyz(:,2),xyz(:,3));  % Transform to spherical coordinates

% --------------------------------------------------------------------
function radio_proj_CB(hObject, handles)
% Set the map projection to either Linear or Orthographic

	do_Ortho = false;
	if ( strcmp(get(hObject,'tag'), 'radio_projOrtho') )	do_Ortho = true;	end

	set_Ortho = false;
	if ( (get(hObject,'Value') && do_Ortho) || (~get(hObject,'Value') && ~do_Ortho) )
		set_Ortho = true;
	end

	if (set_Ortho)
		set(handles.axes1,'xlim',[-1 1], 'ylim',[-1 1],'DataAspectRatio',[1 1 1], 'Vis', 'off')
		set(handles.radio_projLinear,'Val',0),				set(handles.radio_projOrtho,'Val',1)
		set(handles.slider_projOrigLon,'Enable','on')
		set(handles.slider_projOrigLat,'Enable','on')
		set(handles.slider_projOrigPitch,'Enable','on')
		set(handles.h_circ,'Visible','on')
	else
		set(handles.axes1,'xlim',[-180 180], 'ylim',[-90 90],'XtickMode','auto', 'YtickMode','auto', 'Vis', 'on')
		set(handles.radio_projLinear,'Val',1),				set(handles.radio_projOrtho,'Val',0)
		set(handles.slider_projOrigLon,'Enable','off')
		set(handles.slider_projOrigLat,'Enable','off')
		set(handles.slider_projOrigPitch,'Enable','off')
		set(handles.h_circ,'Visible','off')
	end
	swap_proj(handles)

	set(handles.edit_ageSlider,'String','0')
	set(handles.slider_age,'Val',0)

% --------------------------------------------------------------------
function slider_projOrigLon_CB(hObject, handles)
	set(handles.text_projOrigLon,'String',['Lon ' num2str(get(hObject,'Value'))])
	swap_proj(handles)		% In this case there is no projection swapping, only origin updating

% --------------------------------------------------------------------
function slider_projOrigLat_CB(hObject, handles)
	set(handles.text_projOrigLat,'String',['Lat ' num2str(get(hObject,'Value'))])
	swap_proj(handles)		% In this case there is no projection swapping, only origin updating

% --------------------------------------------------------------------
function slider_projOrigPitch_CB(hObject, handles)
	set(handles.text_projOrigPitch,'String',['Pitch ' num2str(get(hObject,'Value'))])
	swap_proj(handles)		% In this case there is no projection swapping, only origin updating

% --------------------------------------------------------------------
function swap_proj(handles)

if (get(handles.radio_projOrtho,'Value'))
    orig = [get(handles.slider_projOrigLat,'Value') get(handles.slider_projOrigLon,'Value') ...
            get(handles.slider_projOrigPitch,'Value')];
    for (k = 1:numel(handles.plates_bak))
        for (l = 1:numel(handles.plates_bak{k}))
            x = get(handles.plates_bak{k}(l),'XData');
            y = get(handles.plates_bak{k}(l),'YData');
            [x,y] = orthographic(x, y, orig);
            set(handles.plates{k}(l),'XData',x,'YData',y)
        end
    end
else
    for (k = 1:numel(handles.plates_bak))
        for (l = 1:numel(handles.plates_bak{k}))
            x = get(handles.plates_bak{k}(l),'XData');
            y = get(handles.plates_bak{k}(l),'YData');
            set(handles.plates{k}(l),'XData',x,'YData',y)
        end
    end    
end

% --------------------------------------------------------------------
function radio_animForward_CB(hObject, handles)
% Just make sure that only one of this radiobuttons pair is on. The animation
% callback will check the status of it and decide on the reconstruction direction.
	if (get(hObject,'Val'))
		set(handles.radio_animBackward,'Val',0)
	else
		set(handles.radio_animBackward,'Val',1)
	end

% --------------------------------------------------------------------
function radio_animBackward_CB(hObject, handles)
	if (get(hObject,'Val'))
		set(handles.radio_animForward,'Val',0)
	else
		set(handles.radio_animForward,'Val',1)
	end

% --- Creates and returns a handle to the GUI figure. 
function rally_plater_LayoutFcn(h1)

set(h1,...
'PaperUnits',get(0,'defaultfigurePaperUnits'),...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'DoubleBuffer','on',...
'MenuBar','none',...
'Name','Rally Plater',...
'NumberTitle','off',...
'Position',[520 171 969 629],...
'Renderer',get(0,'defaultfigureRenderer'),...
'RendererMode','manual',...
'Tag','figure1');

axes('Parent',h1,...
'Units','pixels',...
'Color',get(0,'defaultaxesColor'),...
'ColorOrder',get(0,'defaultaxesColorOrder'),...
'Position',[40 138 911 483],...
'XColor',get(0,'defaultaxesXColor'),...
'XLim',[-180 180],...
'XLimMode','manual',...
'YColor',get(0,'defaultaxesYColor'),...
'YLim',[-90 90],...
'YLimMode','manual',...
'Tag','axes1');

uicontrol('Parent',h1, 'Position',[390 8 271 101], 'Style','frame');

uicontrol('Parent',h1, 'Position',[10 8 221 100],...
'BackgroundColor',[1 1 1],...
'Call',{@main_uiCB,h1,'listbox_stages_CB'},...
'Style','listbox',...
'Tooltip','List of currently available satge poles',...
'Value',1,...
'Tag','listbox_stages');

uicontrol('Parent',h1, 'Position',[10 108 80 21],...
'Call',{@main_uiCB,h1,'push_unSet_CB'},...
'String','Unset All',...
'Tooltip','Unlink all plates of their currently set poles',...
'Tag','push_unSet');

uicontrol('Parent',h1, 'Position',[140 108 80 21],...
'Call',{@main_uiCB,h1,'push_setABS_CB'},...
'String','Set ABS',...
'Tooltip','Load the Absolute plate models',...
'Tag','push_setABS');

uicontrol('Parent',h1, 'Position',[240 53 111 21],...
'Call',{@main_uiCB,h1,'push_loadStages_CB'},...
'String','Load stage poles',...
'Tooltip','Load a file with your own stage poles',...
'Tag','push_loadStages');

uicontrol('Parent',h1, 'Position',[700 9 47 21],...
'BackgroundColor',[1 1 1],...
'String','0',...
'Style','edit',...
'Tooltip','Time of animation start (Ma)',...
'Tag','edit_ageStart');

uicontrol('Parent',h1, 'Position',[824 8 47 21],...
'BackgroundColor',[1 1 1],...
'String','200',...
'Style','edit',...
'Tooltip','Time of animation stop (Ma)',...
'Tag','edit_ageStop');

uicontrol('Parent',h1, 'Position',[764 8 47 21],...
'BackgroundColor',[1 1 1],...
'String','5',...
'Style','edit',...
'Tooltip','Animation step (Ma)',...
'Tag','edit_ageStep');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Position',[900 9 47 21],...
'String','0.5',...
'Style','edit',...
'Tooltip','Frame interval (seconds)',...
'Tag','edit_frameInterval');

uicontrol('Parent',h1,...
'BackgroundColor',[0.899999976158142 0.899999976158142 0.899999976158142],...
'Call',{@main_uiCB,h1,'slider_age_CB'},...
'Position',[700 85 201 16],...
'Style','slider',...
'Tooltip','Slide to select a certain age of reconstruction',...
'Tag','slider_age');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',{@main_uiCB,h1,'edit_ageSlider_CB'},...
'Position',[900 83 47 22],...
'Style','edit',...
'Tooltip','Age of reconstruction (Ma)',...
'Tag','edit_ageSlider');

uicontrol('Parent',h1, 'Position',[703 32 41 15],...
'String','T start','Style','text');

uicontrol('Parent',h1, 'Position',[767 31 41 15],...
'String','T step',...
'Style','text');

uicontrol('Parent',h1,...
'Position',[826 32 41 15],...
'String','T end',...
'Style','text');

uicontrol('Parent',h1,...
'Call',{@main_uiCB,h1,'push_makeStages_CB'},...
'Position',[240 13 111 21],...
'String','Make stage poles',...
'Tooltip','Create stage poles from a finite rotation poles list',...
'Tag','push_makeStages');

uicontrol('Parent',h1, 'Position',[902 32 41 15],...
'String','Delay',...
'Style','text');

uicontrol('Parent',h1, 'Position',[902 105 41 15],...
'String','Time',...
'Style','text');

uicontrol('Parent',h1,...
'BackgroundColor',[0.899999976158142 0.899999976158142 0.899999976158142],...
'Call',{@main_uiCB,h1,'slider_projOrigLon_CB'},...
'Max',180,...
'Min',-180,...
'Position',[500 76 101 14],...
'Style','slider',...
'SliderStep',[0.00277777777777778 0.0138888888888889],...
'Tag','slider_projOrigLon');

uicontrol('Parent',h1,...
'BackgroundColor',[0.899999976158142 0.899999976158142 0.899999976158142],...
'Call',{@main_uiCB,h1,'slider_projOrigLat_CB'},...
'Max',90,...
'Min',-90,...
'Position',[500 48 101 14],...
'Style','slider',...
'SliderStep',[0.00555555555555556 0.0277777777777778],...
'Tag','slider_projOrigLat');

uicontrol('Parent',h1,...
'BackgroundColor',[0.899999976158142 0.899999976158142 0.899999976158142],...
'Call',{@main_uiCB,h1,'slider_projOrigPitch_CB'},...
'Max',90,...
'Min',-90,...
'Position',[500 21 101 14],...
'Style','slider',...
'SliderStep',[0.00555555555555556 0.0277777777777778],...
'Tag','slider_projOrigPitch');

uicontrol('Parent',h1,...
'HorizontalAlignment','left',...
'Position',[605 76 51 15],...
'String','Lon  0',...
'Style','text',...
'Tag','text_projOrigLon');

uicontrol('Parent',h1,...
'HorizontalAlignment','left',...
'Position',[605 47 51 15],...
'String','Lat  0',...
'Style','text',...
'Tag','text_projOrigLat');

uicontrol('Parent',h1,...
'HorizontalAlignment','left',...
'Position',[605 21 51 15],...
'String','Pitch 0',...
'Style','text',...
'Tag','text_projOrigPitch');

uicontrol('Parent',h1,...
'Call',{@main_uiCB,h1,'radio_proj_CB'},...
'Position',[400 64 79 15],...
'String','Linear',...
'Style','radiobutton',...
'Value',1,...
'Tag','radio_projLinear');

uicontrol('Parent',h1,...
'Call',{@main_uiCB,h1,'radio_proj_CB'},...
'Position',[400 34 96 15],...
'String','Orthographic',...
'Style','radiobutton',...
'Tag','radio_projOrtho');

uicontrol('Parent',h1, 'Position',[517 100 111 17],...
'FontSize',10,...
'String','Projection center',...
'Tag', 'txt_Pc',...
'Style','text');

uicontrol('Parent',h1,...
'Call',{@main_uiCB,h1,'radio_animForward_CB'},...
'Position',[700 53 79 15],...
'String','Forward',...
'Style','radiobutton',...
'Tooltip','Animations run from past to present',...
'Value',1,...
'Tag','radio_animForward');

uicontrol('Parent',h1,...
'Call',{@main_uiCB,h1,'radio_animBackward_CB'},...
'Position',[810 53 79 15],...
'String','Backward',...
'Style','radiobutton',...
'Tooltip','Animations run from present to past',...
'Tag','radio_animBackward');

uicontrol('Parent',h1,...
'FontSize',10,...
'Position',[240 93 81 16],...
'String','Stage Poles',...
'Style','text');

uicontrol('Parent',h1, 'Position',[409 99 71 17],...
'FontSize',10,...
'String','Projection',...
'Tag','txt_Prj',...
'Style','text');

function main_uiCB(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
	feval(callback_name,hObject,guidata(h1));

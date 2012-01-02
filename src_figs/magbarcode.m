function magbarcode(chron_file)
% Creates an figure with the Geomagnetic Bar code
%
% CHRON_FILE is the full name of the isochrons file.
% If not given, the default file is used

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

	mir_dirs = getappdata(0,'MIRONE_DIRS');
	if (~isempty(mir_dirs))
		d_path = [mir_dirs.home_dir filesep 'data' filesep];
	else
		d_path = [cd filesep 'data' filesep];
	end

	if (nargin == 0)
		chron_file = [d_path 'Cande_Kent_95.dat'];
	end

	% Define height of the bar code
	tmax = 165;             % Max time (in Ma) represented
	width = 8.0;            % Figure width in cm
	height = 17;            % Figure height in cm
	barHeight_1M = 0.75;    % Scale factor for the height a 1 Ma bar will show in the display
	tscal = tmax * barHeight_1M;
	dy = 1;             % dy is used to shift down all bar code. It represents 1 Ma in the time scale

	% Create figure for bar code display
	F = figure('Units','centimeters',...
		'Position',[1 1 [width height]],...
		'Toolbar','none',...
		'Menubar','none',...
		'Numbertitle','off',...
		'Name','Geomagnetic Bar Code',...
		'doublebuffer','on',...
		'Visible','off',...
		'Color','k',...
		'renderer','Zbuffer',...	% Otherwise the patch command below would set it to OpenGL and
		'Tag','figure1');			% despite what TMW says, in R13 it is awfully bugy.

	% Create axes for bar code display
	axes('Units','normalized',...
		'Position',[0 0 1 1],...
		'XLim',[0 width],...
		'YLim',[0 height],...
		'Color','w',...
		'XTick',[],...
		'YTick',[],...
		'YDir','reverse',...
		'DataAspectRatio',[1 1 1],...
		'DefaultPatchEdgecolor', 'none', ...
		'Tag','axes1');

	handles = guihandles(F);
	handles.scale_vert = barHeight_1M;
	handles.dy = dy;
	guidata(F, handles);

	pos=[0.96 0 .04 1];
	S=['set(gca,''ylim'',[' num2str(tscal-height) ' ' num2str(tscal) ']-get(gcbo,''value''))'];
	uicontrol('style','slider','units','normalized','position',pos,...
		'callback',S,'min',0,'max',tscal-height,'Value',tscal-height);

	fid = fopen(chron_file,'r');
	todos = fread(fid,'*char');     [chron age_start age_end age_txt] = strread(todos,'%s %f %f %s');
	fclose(fid);    clear todos

	y = barHeight_1M * ([age_start'; age_end'] + dy);
	y = y(:);
	x = [repmat(0.1,length(y),1); repmat(0.4,length(y),1)]*width; 
	y = [y; y(end:-1:1)];
	vert = [x y];

	n_ages = length(age_start);
	n2 = 2 * n_ages;
	c1 = (1:n2-1)';     c3 = n2*2 - c1;
	c2 = c3 + 1;        c4 = c1 + 1;
	faces = [c1 c2 c3 c4];

	cor = repmat([0 0 0; 1 1 1],n_ages-1,1);    cor = [cor; [0 0 0]];

	hp = patch('Faces',faces,'Vertices',vert,'FaceVertexCData',cor,'FaceColor','flat');
	set(hp,'ButtonDownFcn',{@bdn_MagBar,handles})

	for (i=1:n_ages)
		if (~strcmp(age_txt(i),'a'))   % Plot anomaly names
			text(width*.42, barHeight_1M*(age_start(i)+dy),age_txt(i),'VerticalAlignment','cap')
		end
	end

	% x = width*[.1 .4 .4 .1 .1];
	% for i=1:length(chron)
	%     patch(x,barHeight_1M * ([age_start(i) age_start(i) age_end(i) age_end(i) age_start(i)] + dy), 'k')
	%     if (~strcmp(age_txt(i),'a'))   % Plot anomaly names
	%         text(width*.42, barHeight_1M*(age_start(i)+dy),age_txt(i),'VerticalAlignment','cap')
	%     end
	% end

	% Draw a vertical line for time ruler
	line('XData',width*[.55 .55],'YData',barHeight_1M*([0+dy tmax+dy]),'LineWidth',2)  

	for (i=0:5:tmax)
		line('XData',width*[.55 .58],'YData',barHeight_1M * ([i i]+dy),'LineWidth',.5)   % Time tick marks
		text(width*.6, barHeight_1M * (i + dy), [num2str(i) ' Ma'])
	end
	for (i=0:tmax)
		line('XData',width*[.55 .565],'YData',barHeight_1M * ([i i]+dy),'LineWidth',.5)   % 1 Ma ticks
	end

	% Draw two vertical lines for known gemagnetic periods
	line('XData',width*[.015 .015],'YData',barHeight_1M*([0+dy 5.4+dy]),'LineWidth',1)  
	line('XData',width*[.08 .08],'YData',barHeight_1M*([0+dy 5.4+dy]),'LineWidth',1)  
	text(width*.04, barHeight_1M*(.35 + dy), 'Bru','rotation',90,'HorizontalAlignment','center')
	line('XData',width*[.015 .08],'YData',barHeight_1M*([.73 .73]+dy),'LineWidth',.5)	% period separator
	text(width*.04, barHeight_1M*(1.62 + dy), 'Mathu','rotation',90,'HorizontalAlignment','center')
	line('XData',width*[.015 .08],'YData',barHeight_1M*([2.5 2.5]+dy),'LineWidth',.5)	% period separator
	text(width*.04, barHeight_1M*(2.95 + dy), 'Gau','rotation',90,'HorizontalAlignment','center')
	line('XData',width*[.015 .08],'YData',barHeight_1M*([3.4 3.4]+dy),'LineWidth',.5)	% period separator
	text(width*.04, barHeight_1M*(4.4 + dy), 'Gilbert','rotation',90,'HorizontalAlignment','center')
	line('XData',width*[.015 .08],'YData',barHeight_1M*([5.4 5.4]+dy),'LineWidth',.5)	% period separator

	%--------------------------
	% Draw two vertical lines for geological periods
	line('XData',width*[.79 .79],'YData',barHeight_1M*([0+dy tmax+dy]),'LineWidth',1)  
	line('XData',width*[.88 .88],'YData',barHeight_1M*([0+dy tmax+dy]),'LineWidth',1)  

	% Plistocene
	text(width*.84, barHeight_1M*(1 + dy), ' Plistocene','rotation',90,'HorizontalAlignment','center')
	line('XData',width*[.79 .88],'YData',barHeight_1M*([2 2]+dy),'LineWidth',.5)   % period separator

	% Pliocene
	text(width*.84, barHeight_1M*(3.5 + dy), 'Pliocene','rotation',90,'HorizontalAlignment','center')
	line('XData',width*[.79 .88],'YData',barHeight_1M*([5 5]+dy),'LineWidth',.5)   % period separator

	% Miocene
	text(width*.84, barHeight_1M*(14.75 + dy), 'Miocene','rotation',90,'HorizontalAlignment','center')
	line('XData',width*[.79 .88],'YData',barHeight_1M*([24.5 24.5]+dy),'LineWidth',.5)   % period separator

	% Oligocene
	text(width*.84, barHeight_1M*(31.25 + dy), 'Oligocene','rotation',90,'HorizontalAlignment','center')
	line('XData',width*[.79 .88],'YData',barHeight_1M*([38.0 38.0]+dy),'LineWidth',.5)   % period separator

	% Eocene
	text(width*.84, barHeight_1M*(46.5 + dy), 'Eocene','rotation',90,'HorizontalAlignment','center')
	line('XData',width*[.79 .88],'YData',barHeight_1M*([55.0 55.0]+dy),'LineWidth',.5)   % period separator

	% Paleocene
	text(width*.84, barHeight_1M*(60.0 + dy), 'Eocene','rotation',90,'HorizontalAlignment','center')
	line('XData',width*[.79 .88],'YData',barHeight_1M*([65.0 65.0]+dy),'LineWidth',.5)   % period separator

	% Maastricthian
	text(width*.84, barHeight_1M*(69.0 + dy), 'Maastricthian','rotation',90,'HorizontalAlignment','center')
	line('XData',width*[.79 .88],'YData',barHeight_1M*([73.0 73.0]+dy),'LineWidth',.5)   % period separator

	% Campanian
	text(width*.84, barHeight_1M*(78.0 + dy), 'Campanian','rotation',90,'HorizontalAlignment','center')
	line('XData',width*[.79 .88],'YData',barHeight_1M*([83.0 83.0]+dy),'LineWidth',.5)   % period separator

	% Santonian
	text(width*.84, barHeight_1M*(85.2 + dy), 'Santonian','rotation',90,'HorizontalAlignment','center')
	line('XData',width*[.79 .88],'YData',barHeight_1M*([87.4 87.4]+dy),'LineWidth',.5)   % period separator

	% Coniacian
	text(width*.84, barHeight_1M*(87.95 + dy), 'Con','rotation',90,'HorizontalAlignment','center')
	line('XData',width*[.79 .88],'YData',barHeight_1M*([88.5 88.5]+dy),'LineWidth',.5)   % period separator

	% Turonian
	text(width*.84, barHeight_1M*(89.75 + dy), 'Turonian','rotation',90,'HorizontalAlignment','center')
	line('XData',width*[.79 .88],'YData',barHeight_1M*([91.0 91.0]+dy),'LineWidth',.5)   % period separator

	% Cenomanian
	text(width*.84, barHeight_1M*(94.25 + dy), 'Cenomanian','rotation',90,'HorizontalAlignment','center')
	line('XData',width*[.79 .88],'YData',barHeight_1M*([97.5 97.5]+dy),'LineWidth',.5)   % period separator

	% Albian
	text(width*.84, barHeight_1M*(100.25 + dy), 'Albian','rotation',90,'HorizontalAlignment','center')
	line('XData',width*[.79 .88],'YData',barHeight_1M*([103.0 103.0]+dy),'LineWidth',.5)   % period separator

	% Aptian
	text(width*.84, barHeight_1M*(110.0 + dy), 'Aptian','rotation',90,'HorizontalAlignment','center')
	line('XData',width*[.79 .88],'YData',barHeight_1M*([119.0 119.0]+dy),'LineWidth',.5)   % period separator

	% Barremian
	text(width*.84, barHeight_1M*(122.0 + dy), 'Barremian','rotation',90,'HorizontalAlignment','center')
	line('XData',width*[.79 .88],'YData',barHeight_1M*([125.0 125.0]+dy),'LineWidth',.5)   % period separator

	% Hauterivian
	text(width*.84, barHeight_1M*(128.0 + dy), 'Hauterivian','rotation',90,'HorizontalAlignment','center')
	line('XData',width*[.79 .88],'YData',barHeight_1M*([131.0 131.0]+dy),'LineWidth',.5)   % period separator

	% Valanginian
	text(width*.84, barHeight_1M*(134.5 + dy), 'Valanginian','rotation',90,'HorizontalAlignment','center')
	line('XData',width*[.79 .88],'YData',barHeight_1M*([138.0 138.0]+dy),'LineWidth',.5)   % period separator

	% Berriasian
	text(width*.84, barHeight_1M*(141.0 + dy), 'Berriasian','rotation',90,'HorizontalAlignment','center')
	line('XData',width*[.79 .88],'YData',barHeight_1M*([144.0 144.0]+dy),'LineWidth',.5)   % period separator

	% Tithonian
	text(width*.84, barHeight_1M*(147.0 + dy), 'Tithonian','rotation',90,'HorizontalAlignment','center')
	line('XData',width*[.79 .88],'YData',barHeight_1M*([150.0 150.0]+dy),'LineWidth',.5)   % period separator

	% Kimmeridgian
	text(width*.84, barHeight_1M*(153.0 + dy), 'Kimmeridgian','rotation',90,'HorizontalAlignment','center')
	line('XData',width*[.79 .88],'YData',barHeight_1M*([156.0 156.0]+dy),'LineWidth',.5)   % period separator

	% Oxfordian
	text(width*.84, barHeight_1M*(159.5 + dy), 'Oxfordian','rotation',90,'HorizontalAlignment','center')
	line('XData',width*[.79 .88],'YData',barHeight_1M*([163.0 163.0]+dy),'LineWidth',.5)   % period separator

	set(F,'Visible','on')

% -----------------------------------------------------------------------------------------
function bdn_MagBar(obj,eventdata,handles)
	stype = get(handles.figure1,'selectiontype');
	if (~strcmp(stype,'open')),		return,		end

	% Create a pico marker
	p = get(handles.axes1,'currentpoint');
	age = p(1,2) / handles.scale_vert - handles.dy;
	yh = 0.24;					% Stick height (in cm)
	%xw = 0.4;					% Stick width (falta dar-lhe uso)
	x = [.8 .6 .4 .4 .55 .8];
	y = [0 -yh/2 -yh/2 yh/2 yh/2 0] + p(1,2);
	h = patch(x,y,'r');

	hp = findobj(handles.axes1,'Tag','Picos');
	n = length(hp);
	ages = ones(n,1) * 1000;    num = zeros(n,1);
	for (k=1:n)
		ages(k) = getappdata(hp(k),'Age');
		num(k) = get(hp(k),'UserData');
	end
	[num,ind] = sort(num); 
	hp = hp(ind);
	ages = sort(ages);
	ind = find(ages > age);
	if (~isempty(ind))			% New marker was inserted between existing ones
		n_pico = ind(1);
		ind = ind + 1;			% Move indices numbers to account for the new inserted one
		opt = 'ins';
	else
		n_pico = n + 1;			% New marker was at the end
		opt = 'add';
	end
	if (~isempty(ind))
		m = 1;
		for (k=ind(1)-1:ind(end)-1)
			set(hp(k),'Userdata',ind(m))
			m = m + 1;
		end
	end

	set(h,'Tag','Picos','UserData',n_pico,'ButtonDownFcn',{@bdn_pico,handles,h,yh})
	setappdata(h,'Age',age)
	listbox_message(num2str(age),n_pico,opt)

% -----------------------------------------------------------------------------------------
function bdn_pico(obj,eventdata,handles,h,yh)
	handles = guidata(handles.figure1);       % We may need to have an updated version of handles
	stype = get(handles.figure1,'selectiontype');
	n_pico = get(h,'UserData');
	hp = findobj(handles.axes1,'Tag','Picos');
	n = length(hp);
	num = get(hp,'UserData');
	if (iscell(num))   num = cat(1,num{:});    end
	[num,ind] = sort(num); 
	hp = hp(ind);
	if strcmp(stype,'open')
		% Remove the pico marker, but we have also to update the 'Userdata' of the remaining picos
		ind = find(num > n_pico);
		if (~isempty(ind))
			num(ind) = num(ind) - 1;
			for (k=1:n)
				set(hp(k),'Userdata',num(k))
			end
		end
		delete(h)
		listbox_message([],n_pico,'del')
		set(handles.figure1,'WindowButtonMotionFcn','','WindowButtonUpFcn','')
	else
		% Prepare things to the ButtonMtion function
		ind = find(n_pico == num);
		if (ind > 1)
				ind_b = ind - 1;        age_b = getappdata(hp(ind_b),'Age');
		else    age_b = 0;
		end
		if (ind ~= num(end))
				ind_a = ind + 1;        age_a = getappdata(hp(ind_a),'Age');
		else    age_a = 1000;
		end
		set(handles.figure1,'WindowButtonMotionFcn',{@wbm_pico,handles,h,yh,n_pico,age_b,age_a}, ...
			'WindowButtonUpFcn',{@wbu_pico,handles});
	end

% -----------------------------------------------------------------------------------------
function wbm_pico(obj,eventdata,handles,h,yh,n_pico,age_b,age_a)
% age_b is the age of the precedent marker (or 0 if it doesn't exist)
% age_a is the age of the following marker (or max_age if it doesn't exist)
	p = get(handles.axes1,'currentpoint');
	age = p(1,2) / handles.scale_vert - handles.dy;
	if (age < 0),	return,		end								% Don't get out of the Bar Code
	if ((age < age_b) || (age > age_a)),	return,		end		% Do not let markers cross each other
	% %disp([num2str(ind_c) ' ' num2str(ind_b) ' ' num2str(ind_a) ' ' num2str(nc)])

	y = [0 -yh/2 -yh/2 yh/2 yh/2 0] + p(1,2);
	set(h,'YData',y)
	setappdata(h,'Age',age)
	listbox_message(num2str(age),n_pico,'rep')

% -----------------------------------------------------------------------------------------
function wbu_pico(obj,eventdata,handles)
	handles = guidata(handles.figure1);       % We may need to have an updated version of handles
	guidata(handles.figure1,handles)
	set(handles.figure1,'WindowButtonMotionFcn','','WindowButtonUpFcn','')

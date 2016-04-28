function varargout = demets_od(varargin)
% Helper window to adjust isochrones for the DeMets & Wilson 2008 Outward Displacement

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

% $Id: demets_od.m 3656 2012-08-09 13:20:31Z j $

	if (isempty(varargin))
		errordlg('DEMETS_OD called with empty arguments','Error'),	return
	end

	hObject = figure('Tag','figure1','Visible','off');
	demets_od_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hObject,'right')

	% Initialize those
	handles.hCallingFig = [];
	handles.isoca1 = [];
	handles.isoca2 = [];
	handles.do_graphic = 0;
	handles.path_continent = [pwd filesep 'continents' filesep];
	handles.OD = [];			% To hold the Outward Displacement data

	handles.hCallingFig = varargin{1};        % This the Mirone's fig handle
	handMir = guidata(handles.hCallingFig);
	if (handMir.no_file)
		errordlg('You didn''t even load a file. What are you expecting then?','ERROR')
		delete(hObject);    return
	end
	if (~handMir.geog)
		errordlg('This operation is currently possible only for geographic type data','ERROR')
		delete(hObject);    return
	end
	handles.last_dir = handMir.last_dir;
	handles.work_dir = handMir.work_dir;
	handles.home_dir = handMir.home_dir;

	%------------ Give a Pro look (3D) to the frame boxes  -------------------------------
	new_frame3D(hObject, [handles.textLines handles.textPole handles.textODf])
	%------------- END Pro look (3D) -----------------------------------------------------

	% Add this figure handle to the carraças list
	plugedWin = getappdata(handles.hCallingFig,'dependentFigs');
	plugedWin = [plugedWin hObject];
	setappdata(handles.hCallingFig,'dependentFigs',plugedWin);

	set(hObject,'Visible','on');
	guidata(hObject, handles);
	if (nargout),	varargout{1} = hObject;		end

% -------------------------------------------------------------------------------------
function edit_pLon_ini_CB(hObject, handles)
	handles.pLon_ini = str2double(get(hObject,'String'));
	guidata(handles.figure1, handles);

% -------------------------------------------------------------------------------------
function edit_pLat_ini_CB(hObject, handles)
	handles.pLat_ini = str2double(get(hObject,'String'));
	guidata(handles.figure1, handles);

% -------------------------------------------------------------------------------------
function edit_pAng_ini_CB(hObject, handles)
	handles.pAng_ini = str2double(get(hObject,'String'));
	guidata(handles.figure1, handles);

% -------------------------------------------------------------------------------------
function push_polesList_CB(hObject, handles)
	fid = fopen([handles.path_continent 'lista_polos.dat'],'rt');
	c = fread(fid,'*char').';
	fclose(fid);
	s = strread(c,'%s','delimiter','\n');

	[s,v] = choosebox('Name','One Euler list',...
						'PromptString','List of poles:',...
						'SelectString','Selected poles:',...
						'ListSize',[450 300],'ListString',s);

	if (v == 1)         % Finite pole
		handles.pLon_ini = s(1);
		handles.pLat_ini = s(2);
		handles.pAng_ini = s(3);
		set(handles.edit_pLon_ini, 'String', num2str(s(1)))
		set(handles.edit_pLat_ini, 'String', num2str(s(2)))
		set(handles.edit_pAng_ini, 'String', num2str(s(3)))
		guidata(handles.figure1,handles)
	elseif (v == 2)     % Stage poles
		%set(handles.edit_polesFile,'String',s)
	end

% -------------------------------------------------------------------------------------
function toggle_pickLines_CB(hObject, handles)
if (get(hObject,'Value'))
    % Test if we have potential target lines and their type
    h_mir_lines = findobj(handles.hCallingFig,'Type','line');		% Fish all objects of type line in Mirone figure
    if (isempty(h_mir_lines))										% We don't have any lines
        str = ['If you hited this button on purpose, than you deserve the following insult.',...
                'You #!|"*!%!?~^)--$&.',... 
                'THERE ARE NO LINES IN THAT FIGURE.'];
        errordlg(str,'Chico Clever');	set(hObject,'Value',0);		return;
    end
    if (length(h_mir_lines) == 1)									% We don't have at least two lines
        str = ['If you hited this button on purpose, than you deserve the following insult.',...
                'You -$&#!*!%!?~^)-.|"/',... 
                'THERE IS ONLY ONE LINE IN THAT FIGURE.'];
        errordlg(str,'Chico Clever');   set(hObject,'Value',0);     return;
    end
    
    % The above test is not enough. For exemple, coastlines are not eligible neither,
    % but is very cumbersome to test all the possibilities of pure non-eligible lines.
    set(handles.hCallingFig,'pointer','crosshair')
    h_line = get_polygon(handles.hCallingFig);		% Get first line handle
    if (~isempty(h_line))
        x = get(h_line,'XData');		y = get(h_line,'YData');
        handles.isoca1 = [x(:) y(:)];
        set(handles.text_lineA,'String','Got Line A', 'ForegroundColor',[0 0.75 0])
    else
        handles.isoca1 = [];
        set(handles.text_lineA,'String','Line A (not yet)', 'ForegroundColor','r')
    end
    h_line = get_polygon(handles.hCallingFig);		% Get second line handle
    if (~isempty(h_line))
        x = get(h_line,'XData');		y = get(h_line,'YData');
        handles.isoca2 = [x(:) y(:)];
        set(handles.text_lineB,'String','Got Line B', 'ForegroundColor',[0 0.75 0])
    else
        handles.isoca2 = [];
        set(handles.text_lineB,'String','Line B (not yet)', 'ForegroundColor','r')
    end
    set(handles.hCallingFig,'pointer','arrow')
    if (isempty(handles.isoca1) || isempty(handles.isoca2))
        set(hObject,'Value',0)
        handles.do_graphic = false;
    else
        handles.do_graphic = true;
    end
    set(hObject,'Value',0)
    figure(handles.figure1)			% Bring this figure to front again
else        % What should I do?
    %handles.do_graphic = 0;
end
guidata(handles.figure1, handles);

% -------------------------------------------------------------------------------------
function push_OK_CB(hObject, handles)
% ...
	if (isempty(handles.isoca1) || isempty(handles.isoca2))
		errordlg('Displace what? It would help if you provide me TWO lines.','Chico Clever')
		return
	end
	if (isempty(handles.pLon_ini) || isempty(handles.pLat_ini) || isempty(handles.pAng_ini))
		errordlg('I need to know a Euler pole for this exercise.','Error'),		return
	end
	if (isempty(handles.OD))
		errordlg('Need file with the Outward displacements','Chico Clever'),	return
	end

	apply_OD(handles)

% -------------------------------------------------------------------------------
function apply_OD(handles)
% Apply the Outward Displacement (in)correction request

	D2R = pi / 180;		R2D = 180 / pi;
	isoca1 = handles.isoca1 * D2R;		% A maluca
	isoca2 = handles.isoca2 * D2R;		% A fixa
	OD = handles.OD;

	[rlon1,rlat1] = rot_euler(isoca1(:,1), isoca1(:,2),handles.pLon_ini*D2R,handles.pLat_ini*D2R,handles.pAng_ini*D2R,'rad',-1);
	[rlon2,rlat2] = rot_euler(isoca2(:,1), isoca2(:,2),handles.pLon_ini*D2R,handles.pLat_ini*D2R,-handles.pAng_ini*D2R,'rad',-1);
	shifts_1 = interp1(OD(:,1)*D2R, OD(:,2), rlat1, 'linear', 0);		% Shifts in km
	shifts_2 = interp1(OD(:,1)*D2R, OD(:,2), rlat2, 'linear', 0);

	isoca1 = isoca1 * R2D;			isoca2 = isoca2 * R2D;
	rlat1  = rlat1 * R2D;			rlon1  = rlon1 * R2D;
	rlat2  = rlat2 * R2D;			rlon2  = rlon2 * R2D;
	[s,a12] = vdist(isoca1(:,2),isoca1(:,1),rlat1,rlon1);	% coords must be in degrees here
	for (k = 1:numel(rlon1))
		[rlat1(k), rlon1(k)] = vreckon(rlat1(k), rlon1(k), shifts_1(k)*1000, a12(k), 1);
	end
	[s,a12] = vdist(isoca2(:,2),isoca2(:,1),rlat2,rlon2);
	for (k = 1:numel(rlon2))
		[rlat2(k), rlon2(k)] = vreckon(rlat2(k), rlon2(k), shifts_2(k)*1000, a12(k), 1);
	end

	isoca1 = isoca1 * D2R;			isoca2 = isoca2 * D2R;	% Put them back in radians
	rlat1  = rlat1 * D2R;			rlon1  = rlon1 * D2R;
	rlat2  = rlat2 * D2R;			rlon2  = rlon2 * D2R;
	% Now invert rotations to get the "original" corrected positions
	[isoca1(:,1),isoca1(:,2)] = rot_euler(rlon1, rlat1, handles.pLon_ini*D2R,handles.pLat_ini*D2R,-handles.pAng_ini*D2R,'rad',-1);
	[isoca2(:,1),isoca2(:,2)] = rot_euler(rlon2, rlat2, handles.pLon_ini*D2R,handles.pLat_ini*D2R,handles.pAng_ini*D2R,'rad',-1);

	h1 = line('parent',get(handles.hCallingFig,'CurrentAxes'),'XData',isoca1(:,1)/D2R,'YData',isoca1(:,2)/D2R,'Color','r');
	h2 = line('parent',get(handles.hCallingFig,'CurrentAxes'),'XData',isoca2(:,1)/D2R,'YData',isoca2(:,2)/D2R,'Color','y');
	setappdata(h1,'LineInfo','OD adjusted Line');
	setappdata(h2,'LineInfo','OD adjusted Line');
	draw_funs(h1,'isochron',{'OD adjusted Line'},'Userdata',1)
	draw_funs(h2,'isochron',{'OD adjusted Line'},'Userdata',1)

% --------------------------------------------------------------------------------
function edit_ODfile_CB(hObject, handles)
	fname = get(hObject,'String');
	if isempty(fname),		return,		end
	% Let the push_ODfile_CB do all the work
	push_ODfile_CB(handles.push_ODfile, handles, fname)

% --------------------------------------------------------------------------------
function push_ODfile_CB(hObject, handles, opt)
% ...
	if (nargin == 3)    fname = opt;
	else                opt = [];
	end

	if (isempty(opt))		% Otherwise we already know fname from the 4th input argument
		[FileName,PathName] = put_or_get_file(handles, ...
			{'*.dat;*.DAT', 'Data files (*.dat)';'*.*', 'All Files (*.*)'},'Enter OD file','get');
		if (isequal(FileName,0)),		return,		end
		fname = [PathName FileName];
	end

	if (exist(fname,'file') ~= 2)
		errordlg(['Outward displacement file ' fname ' does not exist. Ignoring OD request'],'Error')
		fname = [];
	end
	out = text_read(fname);
	if (isempty(out) || size(out,2) == 1)
		errordlg(['File ' fname ' doesn''t have any recognized numeric data, or one column only.'],'Error');
		fname = [];
	end

	if (isempty(fname))
		handles.OD = [];
		set(handles.edit_ODfile, 'Str', '')
	else
		handles.OD = out;		% Store the OD data
		set(handles.edit_ODfile, 'Str', fname)
	end
	guidata(handles.figure1, handles);
	
% -----------------------------------------------------------------------------------
function figure1_KeyPressFcn(hObject, eventdata)
	handles = guidata(hObject);
	if isequal(get(hObject,'CurrentKey'),'escape')
		delete(handles.figure1);
	end

% -----------------------------------------------------------------------
function demets_od_LayoutFcn(h1)

set(h1, 'Pos',[520 540 440 260],...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'KeyPressFcn',@figure1_KeyPressFcn,...
'MenuBar','none',...
'Name','Outward Displacements',...
'NumberTitle','off',...
'RendererMode','manual',...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1');

uicontrol('Parent',h1, 'Pos',[10 39 421 61],'Style','frame', 'Tag','frame3');
uicontrol('Parent',h1, 'Pos',[10 118 421 62],'Style','frame', 'Tag','frame2');
uicontrol('Parent',h1, 'Pos',[10 199 421 54],'Style','frame', 'Tag','frame1');

uicontrol('Parent',h1, 'Pos',[150 214 141 23],...
'Call',{@demets_od_uiCB,h1,'toggle_pickLines_CB'},...
'String','Pick lines from Figure',...
'Style','togglebutton',...
'TooltipString','Select the two lines from a Mirone figure',...
'Tag','toggle_pickLines');

uicontrol('Parent',h1, 'Pos',[19 131 61 21],...
'BackgroundColor',[1 1 1],...
'Call',{@demets_od_uiCB,h1,'edit_pLon_ini_CB'},...
'String','',...
'Style','edit',...
'TooltipString','Euler pole Longitude',...
'Tag','edit_pLon_ini');

uicontrol('Parent',h1, 'Pos',[109 131 61 21],...
'BackgroundColor',[1 1 1],...
'Call',{@demets_od_uiCB,h1,'edit_pLat_ini_CB'},...
'String','',...
'Style','edit',...
'TooltipString','Euler pole Latitude',...
'Tag','edit_pLat_ini');

uicontrol('Parent',h1, 'Pos',[199 131 61 21],...
'BackgroundColor',[1 1 1],...
'Call',{@demets_od_uiCB,h1,'edit_pAng_ini_CB'},...
'String','',...
'Style','edit',...
'TooltipString','Euler pole Angle',...
'Tag','edit_pAng_ini');

uicontrol('Parent',h1, 'Pos',[288 130 121 23],...
'Call',{@demets_od_uiCB,h1,'push_polesList_CB'},...
'String','Poles selector',...
'TooltipString','Select a pole from the default list',...
'Tag','push_polesList');

uicontrol('Parent',h1, 'Pos',[190 171 60 18],...
'FontSize',10,...
'String','Pole',...
'Style','text',...
'Tag','textPole');

uicontrol('Parent',h1, 'Pos',[189 243 60 16],...
'FontSize',10,...
'String','Lines',...
'Style','text',...
'Tag','textLines');

uicontrol('Parent',h1, 'Pos',[23 155 51 15], 'String','Longitude', 'Style','text');
uicontrol('Parent',h1, 'Pos',[115 155 51 15],'String','Latitude', 'Style','text');
uicontrol('Parent',h1, 'Pos',[204 155 51 15],'String','Angle', 'Style','text');

uicontrol('Parent',h1, 'Pos',[134 90 180 18],...
'FontSize',10,...
'String','Outward Displacement file',...
'Style','text',...
'Tag','textODf');

uicontrol('Parent',h1, 'Pos',[18 50 371 21],...
'BackgroundColor',[1 1 1],...
'Call',{@demets_od_uiCB,h1,'edit_ODfile_CB'},...
'HorizontalAlignment','left',...
'Style','edit',...
'TooltipString','Name of a file with the DeMets and Wilson outward displacements',...
'Tag','edit_ODfile');

uicontrol('Parent',h1, 'Pos',[389 50 21 21],...
'Call',{@demets_od_uiCB,h1,'push_ODfile_CB'},...
'FontSize',10,...
'FontWeight','bold',...
'String','...',...
'Tag','push_ODfile');

uicontrol('Parent',h1, 'Pos',[20 71 80 18],...
'FontSize',10,...
'String','File name',...
'Style','text');

uicontrol('Parent',h1, 'Pos',[19 215 120 21],...
'FontSize',12,...
'ForegroundColor',[1 0 0],...
'HorizontalAlignment','right',...
'String','Line A (not yet)',...
'Style','text',...
'Tag','text_lineA');

uicontrol('Parent',h1, 'Pos',[304 215 120 21],...
'FontSize',12,...
'ForegroundColor',[1 0 0],...
'HorizontalAlignment','left',...
'String','Line B (not yet)',...
'Style','text',...
'Tag','text_lineB');

uicontrol('Parent',h1, 'Pos',[355 7 76 23],...
'Call',{@demets_od_uiCB,h1,'push_OK_CB'},...
'FontSize',10,...
'FontWeight','bold',...
'String','OK',...
'Tag','push_OK');

function demets_od_uiCB(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
	feval(callback_name,hObject,guidata(h1));

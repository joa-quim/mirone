function varargout = interp_chimoce(varargin)
% Helper figure to ...

%	Copyright (c) 2004-2018 by J. Luis
%
%             DO WHAT THE FUCK YOU WANT TO PUBLIC LICENSE
%                     Version 2, December 2004
% 
%  Everyone is permitted to copy and distribute verbatim or modified
%  copies of this license document, and changing it is allowed as long
%  as the name is changed.
% 
%             DO WHAT THE FUCK YOU WANT TO PUBLIC LICENSE
%    TERMS AND CONDITIONS FOR COPYING, DISTRIBUTION AND MODIFICATION
% 
%   0. You just DO WHAT THE FUCK YOU WANT TO.
%
%	Contact info: w3.ualg.pt/~jluis/mirone
% --------------------------------------------------------------------

% $Id: interp_chimoce.m 10353 2018-03-30 22:33:50Z j $

	if (nargin > 1 && ischar(varargin{1}))
		gui_CB = str2func(varargin{1});
		[varargout{1:nargout}] = feval(gui_CB,varargin{2:end});
	else
		h = interp_chimoce_OF(varargin{:});
		if (nargout),	varargout{1} = h;   end
	end

% ---------------------------------------------------------------------------------
function hObject = interp_chimoce_OF(varargin)
% varargin is a 3 element array because this function is called as a callback in load_xyz
% Element of interest is the 3rd element that holds the line handle where necessary info is stored.

	hObject = interp_chimoce_LayoutFcn;
	handles = guihandles(hObject);

	hMirFig = get(get(varargin{3}, 'Parent'), 'Parent');
	move2side(hMirFig, hObject)

	set(hObject,'Vis','on');

	data = getappdata(varargin{3}, 'all_data');
	ind = strcmp(data.var_names, 'dist');
	ind = ind | strcmp(data.var_names, 'station');
	ind = ind | strcmp(data.var_names, 'depth');
	ind = ind | strcmp(data.var_names, 'st_num');
	ind = ind | strcmp(data.var_names, 'lon');
	ind = ind | strcmp(data.var_names, 'lat');
	set(handles.listbox_vars, 'Str', data.var_names(~ind))
	set(handles.popup1, 'Str', data.st_num)
	
	% Add this figure handle to the carraças list
	plugedWin = getappdata(hMirFig,'dependentFigs');
	plugedWin = [plugedWin hObject];
	setappdata(hMirFig,'dependentFigs',plugedWin);

	handles.hMirFig = hMirFig;
	handles.hLine   = varargin{3};
	guidata(hObject, handles);

	if (nargin > 3),	external_drive(handles, 'interp_chimoce', varargin{4:end}),	end

% ---------------------------------------------------------------------------
function radio_profile_CB(hObject, handles)
	if (~get(hObject,'Value')),		set(hObject,'Value',1),		return,		end
	set(handles.radio_interpolate, 'Val', 0)
	set([handles.check_allStations handles.popup1], 'Enable', 'on')

% ---------------------------------------------------------------------------
function radio_interpolate_CB(hObject, handles)
	if (~get(hObject,'Value')),		set(hObject,'Value',1),		return,		end
	set(handles.radio_profile, 'Val', 0)
	set([handles.check_allStations handles.popup1], 'Enable', 'off')

% ---------------------------------------------------------------------------
function push_OK_CB(hObject, handles)
% ...
	val = get(handles.listbox_vars, 'Val');		str = get(handles.listbox_vars, 'Str');
	var = str{val};
	all = getappdata(handles.hLine, 'all_data');

	if (get(handles.radio_profile, 'Val'))
		val = get(handles.popup1, 'Val');
		st = all.st_num(val);
		stations = all.station;
		if (~get(handles.check_allStations, 'Val'))
			ind = find(stations == st);
			y = all.depth(ind);
			x = all.(var);
			x = x(ind);
			ecran(x,y, [var ' vs depth'])
		else
			tmp = cell(numel(all.st_num), 2);		np = 1;
			for (k = 1:numel(all.st_num))
				ind = find(stations == all.st_num(k));
				tmp{k, 1} = all.depth(ind);
				x = all.(var);
				tmp{k, 2} = x(ind);
				n = numel(tmp{k, 1});
				if (n > np),	np = n;		end
			end
			% Now that we know the size of each profile we can copy them all into a matrix of YYs
			XX = zeros(np, numel(all.st_num)) * NaN;
			YY = zeros(np, numel(all.st_num)) * NaN;
			for (k = 1:numel(all.st_num))
				YY(1:numel(tmp{k}), k) = tmp{k,1}(:);
				XX(1:numel(tmp{k}), k) = tmp{k,2}(:);
			end
			ecran(XX,YY, [var ' vs depth'])
		end
	else
		x = all.dist;
		y = all.depth;
		z = all.(var);
		griding_mir(handles.hMirFig, 'surface', [min(x) max(x) min(y) max(y)], [x y z])
	end

% ---------------------------------------------------------------------------
function h1 = interp_chimoce_LayoutFcn()

h1 = figure('Position',[748 968 270 170],...
'Color', get(0,'factoryUicontrolBackgroundColor'),...
'MenuBar','none',...
'Name','Interp ChimOce',...
'NumberTitle','off',...
'Resize','off',...
'HandleVisibility','callback',...
'Vis','off',...
'Tag','figure1');

uicontrol('Parent',h1, 'Position',[10 10 131 140.5],...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'String',{'Listbox'},...
'Style','listbox',...
'Value',1,...
'BackgroundColor',[1 1 1],...
'TooltipString','Select the variable do analize',...
'Tag','listbox_vars');

uicontrol('Parent',h1, 'Position',[26 151 90 16],...
'String','Variables in file',...
'Style','text');

uicontrol('Parent',h1, 'Position',[166 103 45 16],...
'HorizontalAlignment','right',...
'String','Stations',...
'Style','text');

uicontrol('Parent',h1, 'Position',[211.5 104.5 50 16],...
'String',' ',...
'Style','popupmenu',...
'Value',1,...
'BackgroundColor',[1 1 1],...
'TooltipString','Select station to plot',...
'Tag','popup1');

uicontrol('Parent',h1, 'Position',[151 123 75 17],...
'Callback',@interpChimOce_uiCB,...
'String','Plot profile',...
'Style','radiobutton',...
'TooltipString','Plot profile(s) of selected variable',...
'Tag','radio_profile');

uicontrol('Parent',h1, 'Position',[171 83 80 17],...
'String','All stations',...
'Style','checkbox',...
'TooltipString','Plot data from all stations',...
'Tag','check_allStations');

uicontrol('Parent',h1, 'Position',[151 51 110 17],...
'Callback',@interpChimOce_uiCB,...
'String','Interpolate section',...
'Style','radiobutton',...
'Value',1,...
'TooltipString','Call tool to interpolate selected var along the section.',...
'Tag','radio_interpolate');

uicontrol('Parent',h1, 'Position',[190 10 70 21],...
'Callback',@interpChimOce_uiCB,...
'String','OK',...
'Tag','push_OK',...
'FontSize',9,...
'FontWeight','bold');

function interpChimOce_uiCB(hObject, evt)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));
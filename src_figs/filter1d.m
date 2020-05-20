function varargout = filter1d(varargin)
% Helper window to GMT filter1d.

%	Copyright (c) 2004-2020 by J. Luis
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

	% Here varargin{1} is a line handle
	if (isempty(varargin)),		return,		end

	hObject = figure('Vis','off');
	filter1d_LayoutFcn(hObject);
	handles = guihandles(hObject);

	handXY = guidata(get(get(varargin{1}, 'Parent'), 'Parent'));	% handles of the Ecran fig
	handles.hXYFig = handXY.figure1;
	move2side(handXY.figure1, hObject,'bottom')
	handles.filt_type = '-Fb';	% Default to first one

	handles.xx = get(handXY.hLine(1),'XData');		handles.yy = get(handXY.hLine(1),'YData');
	handles.hFiltLine = [];		% Will hold the handle of the filtered line
	handles.cmd = '';			% To hold the filter1d command

	set(hObject,'Visible','on');
	guidata(hObject, handles);
	inc = mean(diff(double(handles.xx)));
	set(handles.edit_inc, 'Tooltip', sprintf('Data average separation is %.12g But leave blank if data is equispaced', inc))
	if (nargout),	varargout{1} = hObject;		end

%----------------------------------------------------------------------
function popup_ftype_CB(hObject, handles)
	val = get(hObject,'Value');		str = get(hObject, 'String');
	switch str{val}
		case 'boxcar',				handles.filt_type = '-Fb';
		case 'cosine arch',			handles.filt_type = '-Fc';
		case 'gaussian',			handles.filt_type = '-Fg';
		case 'median',				handles.filt_type = '-Fm';
		case 'maximum likelihood',	handles.filt_type = '-Fp';
		case 'lower',				handles.filt_type = '-Fl';
		case 'upper',				handles.filt_type = '-Fu';
	end
	guidata(hObject, handles);

%----------------------------------------------------------------------
function edit_width_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx) || xx < 0),	set(hObject, 'String', ''),	return,	end

%----------------------------------------------------------------------
function edit_inc_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx) || xx < 0),	set(hObject, 'String', ''),	return,	end
	handles.data_inc = xx;
	guidata(hObject,handles)

%----------------------------------------------------------------------
function radio_low_CB(hObject, handles)
	if (~get(hObject,'Val')),	set(hObject,'Val',1),	return,		end
	set(handles.radio_high,'Val',0)

%----------------------------------------------------------------------
function radio_high_CB(hObject, handles)
	if (~get(hObject,'Val')),	set(hObject,'Val',1),	return,		end
	set(handles.radio_low,'Val',0)

%----------------------------------------------------------------------
function push_Apply_CB(hObject, handles)
% ...
	filt_W = get(handles.edit_width,'String');
	if (isempty(filt_W))
		errordlg('Filter width is mandatory', 'Error'),		return
	end

	if (get(handles.check_robust, 'Val'))
		handles.filt_type(3) = upper(handles.filt_type(3));
	else
		handles.filt_type(3) = lower(handles.filt_type(3));
	end
	xx = str2double(get(handles.edit_inc,'String'));
	opt_D = '';
	if (~isnan(xx)),	opt_D = sprintf('-D%.12g', xx);		end
	opt_F = [handles.filt_type ddewhite(filt_W)];
	if (get(handles.radio_high, 'Val')),	opt_F = [opt_F '+h'];	end

	cmd = sprintf('filter1d %s %s -E', opt_F, opt_D);
	if (strcmp(cmd, handles.cmd)),		return,		end		% If nothing changed, nothing new to do

	D = gmtmex(cmd, [handles.xx(:) handles.yy(:)]);

	if (isempty(handles.hFiltLine))		% First run
		handXY = guidata(handles.hXYFig);
		handles.hFiltLine = line('XData',D.data(:,1),'YData',D.data(:,2), 'Parent', handXY.axes1, ...
			'LineStyle', '--', 'Color','r', 'LineWidth', get(handXY.hLine(1), 'LineWidth'));
	else								% Just update filtered line
		set(handles.hFiltLine, 'YData',D.data(:,2))
	end
	handles.cmd = cmd;
	guidata(handles.figure1, handles)

%----------------------------------------------------------------------
function push_applyNreturn_CB(hObject, handles)
	push_Apply_CB([], handles)
	handXX = guidata(handles.hXYFig);
	set(handles.hFiltLine, 'LineStyle', '-')
	draw_funs([], 'set_line_uicontext_XY', handles.hFiltLine)	% Set lines's uicontextmenu
	handXX.hLine(end+1) = handles.hFiltLine;		% Updated this one
	guidata(handXX.figure1, handXX)
	delete(handles.figure1)

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, evt)
	handles = guidata(hObject);
	delete(handles.hFiltLine)		% Delete the filter line
	delete(handles.figure1);

%----------------------------------------------------------------------
function h1 = filter1d_LayoutFcn(h1)

	set(h1, 'Position',[707 995 291 121],...
		'CloseRequestFcn',@figure1_CloseRequestFcn,...
		'Color',get(0,'factoryUicontrolBackgroundColor'),...
		'MenuBar','none',...
		'Name','filter1d',...
		'NumberTitle','off',...
		'Resize','off',...
		'HandleVisibility','callback',...
		'Tag','figure1');

	uicontrol('Parent',h1, 'Position',[10 81 110 21],...
		'Callback',@filter1d_uiCB,...
		'String',{'boxcar'; 'cosine arch'; 'gaussian'; 'median'; 'maximum likelihood'; 'lower'; 'upper'},...
		'Style','popupmenu',...
		'Value',1,...
		'BackgroundColor',[1 1 1],...
		'Tag','popup_ftype');

	uicontrol('Parent',h1, 'Position',[120 81 70 21],...
		'Callback',@filter1d_uiCB,...
		'Style','edit',...
		'BackgroundColor',[1 1 1],...
		'Tooltip','Specify filter full width in xx data units',...
		'Tag','edit_width');

	uicontrol('Parent',h1, 'Position',[24.6 102 68 14],...
		'String','Filter type',...
		'Style','text');

	uicontrol('Parent',h1, 'Position',[122.2 102 68.4 14],...
		'String','Filter width',...
		'Style','text');

	uicontrol('Parent',h1,...
		'Style','edit',...
		'Position',[210 81 70 21],...
		'BackgroundColor',[1 1 1],...
		'Tooltip','increment is used when series is NOT equidistantly sampled. Then increment will be the abscissae resolution',...
		'Tag','edit_inc');

	uicontrol('Parent',h1, 'Position',[10 43 70 18],...
		'Callback',@filter1d_uiCB,...
		'String','Low pass',...
		'Style','radiobutton',...
		'Value',1,...
		'Tooltip','Perform low-pass filtering',...
		'Tag','radio_low');

	uicontrol('Parent',h1, 'Position',[81 43 70 18],...
		'Callback',@filter1d_uiCB,...
		'String','High pass',...
		'Style','radiobutton',...
		'Tooltip','Perform high-pass filtering',...
		'Tag','radio_high');

	uicontrol('Parent',h1, 'Position',[211 102 70 14],...
		'String','xx resolution',...
		'Style','text');

	uicontrol('Parent',h1, 'Position',[171.6 43 60 18],...
		'String','Robust',...
		'Style','checkbox',...
		'Tooltip','Use robust filter',...
		'Tag','check_robust');

	uicontrol('Parent',h1, 'Position',[190 11 91 21],...
		'Callback',@filter1d_uiCB,...
		'String','Apply n return',...
		'FontSize',9,...
		'FontWeight','bold',...
		'Tag','push_applyNreturn');

	uicontrol('Parent',h1, 'Position',[110 11 70 21],...
		'Callback',@filter1d_uiCB,...
		'String','Apply',...
		'FontSize',9,...
		'FontWeight','bold',...
		'Tag','push_Apply');

function filter1d_uiCB(hObject, eventdata)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));

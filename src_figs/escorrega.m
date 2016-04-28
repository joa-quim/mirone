function varargout = escorrega(varargin)
% Figure with only a slider+edit box. Used for showing VE and in fututre to set transaparency and alike

%	Copyright (c) 2004-2013 by J. Luis
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

% $Id: escorrega.m 3818 2013-01-07 16:23:39Z j $

	hObject = figure('Vis','off');
	escorrega_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hObject,'center')

	opt = varargin{1};
	if (strncmp(opt,'vertical_exaggeration',3))
		handles.hCallingFig = varargin{2};
		handEcran = guidata(handles.hCallingFig);
		posAx1 = get(handEcran.axes1,'Pos');			% Assume the units is 'normalized'
		DAR = get(handEcran.axes1, 'DataAspectRatio');
		VertExa = DAR(1)/DAR(2) * posAx1(4)/posAx1(3);
		screen = get(0,'ScreenSize');
		posFigInPix = get(handEcran.figure1,'Pos');		% Assume the units are 'pixels'
		maxHeightInPix = posAx1(4) * (screen(4)-100);	% Maximum height in pixels that axes1 may have
		curHeightInPix = posAx1(4) * posFigInPix(4);
		maxVE = VertExa * maxHeightInPix / curHeightInPix;
		set(handles.edit, 'String', VertExa)
		set(handles.slider, 'Min',VertExa, 'Max',maxVE, 'Val',VertExa)	% Is 'off' anyway
		set(handles.text,'Str','Vertical exaggeration')
		set(handles.figure1,'Name','VE')
	end

	guidata(hObject, handles);
	set(hObject,'Visible','on');
	if (nargout),	varargout{1} = hObject;		end

% ---------------------------------------------------------------------
function slider_CB(hObject, handles)
% Not used yet
	val = get(hObject,'Value');
	set(handles.edit,'String', sprintf('%.2f',val))

% ---------------------------------------------------------------------
function edit_CB(hObject, handles)
% ...
	val = str2double(get(hObject,'String'));
	if (isnan(val) || val < get(handles.slider,'Min') || val > get(handles.slider,'Max'))
		set(hObject,'String',handles.val),		return
	end
	set(handles.slider,'Val',val)


% --- Creates and returns a handle to the GUI figure. 
function escorrega_LayoutFcn(h1)

set(h1,...
'Units','centimeters',...
'Position',[13.7214409375 20.1194923958333 7.164760104166666 1.08396739583333],...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'MenuBar','none',...
'Name','escorrega',...
'NumberTitle','off',...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1');

uicontrol('Parent',h1, 'Position',[5 5 201 17],...
'BackgroundColor',[0.9 0.9 0.9],...
'Call',@escorrega_uiCB,...
'Style','slider',...
'Enable','off',...
'Tag','slider');

uicontrol('Parent',h1, 'Position',[206 4 61 20],...
'BackgroundColor',[1 1 1],...
'Call',@escorrega_uiCB,...
'String','',...
'Style','edit',...
'Tag','edit');

uicontrol('Parent',h1, 'Position',[21 25 121 14],...
'HorizontalAlignment','left',...
'String','Text',...
'Style','text',...
'Tag','text');

function escorrega_uiCB(hObject, evt)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));

function varargout = about_box(varargin)
% Display the "About" thing

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
	about_box_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hObject,'center')

    handMir = varargin{1};

	% Load background images
	logo  = imread([handMir.path_data 'logo.png']);
	fundo = imread([handMir.path_data 'cafe_ico.png']);
	image(logo,'Parent',handles.axes1)
	set(handles.axes1,'Visible', 'off');
	image(fundo,'Parent',handles.axes2)
	set(handles.axes2,'Visible', 'off');
	bg_color = double( squeeze(logo(1,1,:)) ) / 255;
	set(hObject, 'Color', bg_color)

    set(handles.JL,'BackgroundColor', bg_color)
    set(handles.text_prog,'String',varargin{2}, 'BackgroundColor', bg_color)
	if (numel(varargin) == 3) 
        str_prog = 'Mirone. The ultimate indiscreet grid viewer';
        str_analpha = ['Mirone, Version ' varargin{3}];
        str_url = 'w3.ualg.pt/~jluis/mirone';
        set(handles.text_ProgName,'String',str_prog, 'BackgroundColor', bg_color)
        set(handles.text_AnalphaBeta,'String',str_analpha, 'BackgroundColor', bg_color)
        set(handles.text_url,'String',str_url, 'BackgroundColor', bg_color)
	end
	
	set(hObject,'Vis','on');
	if (strncmp(computer,'PC',2))
		WindowAPI(hObject, 'TopMost')
	end
	guidata(hObject,handles)
	if (nargout),    varargout{1} = hObject;    end


% --- Executes on key press over figure1 with no controls selected.
function figure1_KeyPressFcn(hObject, eventdata)
% Check for "escape"
if isequal(get(hObject,'CurrentKey'),'escape')
	handles = guidata(hObject);
	delete(handles.figure1)
end

% --- Creates and returns a handle to the GUI figure. 
function about_box_LayoutFcn(h1)

set(h1, 'Position',[303 228 269 380],...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'KeyPressFcn',@figure1_KeyPressFcn,...
'MenuBar','none',...
'Toolbar', 'none',...
'Name','about_box',...
'NumberTitle','off',...
'Resize','off',...
'Tag','figure1');

uicontrol('Parent',h1,'Position',[13 99 241 30],...
'Max',10,...
'String',{'M-GMT is a MATLAB GUI to the Generic Mapping'; 'Tools (GMT) software'; '' },...
'Style','text',...
'Tag','text_ProgName');

axes('Parent',h1,'Units','pixels', 'Position',[15 139 240 240],...
'Color',get(0,'defaultaxesColor'),...
'Tag','axes1');

uicontrol('Parent',h1,'Position',[11 69 293 20],...
'FontSize',11,...
'HorizontalAlignment','left',...
'String','M-GMT, Version Analpha-Beta 0.5 (Built 4)',...
'Style','text',...
'Tag','text_AnalphaBeta');

axes('Parent',h1,'Units','pixels',...
'Color',get(0,'defaultaxesColor'),...
'Position',[10 49 75 18],'Tag','axes2');

uicontrol('Parent',h1,'Position',[86 45 45 21],...
'FontSize',12, 'HorizontalAlignment','left',...
'String','J Luis',...
'Style','text',...
'Tag','JL');

uicontrol('Parent',h1, 'Position',[138 47 113 17],...
'ForegroundColor',[0 0.501960784313725 0.250980392156863],...
'String','w3.ualg.pt/~jluis/m_gmt',...
'Style','text',...
'Tag','text_url');

uicontrol('Parent',h1,'Position',[10 91 251 1],'Style','frame');
uicontrol('Parent',h1,'Position',[10 41 251 1],'Style','frame');

uicontrol('Parent',h1, 'Position',[10 16 201 15],...
'HorizontalAlignment','left',...
'Style','text','Tag','text_prog');

uicontrol('Parent',h1, 'Position',[216 11 46 21],...
'Call', 'delete(gcf)', ...
'String','OK');

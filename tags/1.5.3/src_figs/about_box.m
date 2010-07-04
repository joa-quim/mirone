function varargout = about_box(varargin)
% M-File changed by desGUIDE  
% varargin   command line arguments to about_box (see VARARGIN)

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
	about_box_LayoutFcn(hObject);
	handles = guihandles(hObject);
	movegui(hObject,'northeast')

    handMir = varargin{1};

	% Load background images
	logo  = imread([handMir.path_data 'logo.png']);
	fundo = imread([handMir.path_data 'cafe_ico.png']);
	image(logo,'Parent',handles.axes1)
	set(handles.axes1,'Visible', 'off');
	image(fundo,'Parent',handles.axes2)
	set(handles.axes2,'Visible', 'off');

    set(handles.text_prog,'String',varargin{2})
	if (numel(varargin) == 3) 
        str_prog = 'Mirone. The ultimate indescrete grid viewer';
        str_analpha = ['Mirone, Version ' varargin{3}];
        str_url = 'w3.ualg.pt/~jluis/mirone';
        set(handles.text_ProgName,'String',str_prog)
        set(handles.text_AnalphaBeta,'String',str_analpha)
        set(handles.text_url,'String',str_url)
	end
	
	set(hObject,'Visible','on');
	guidata(hObject,handles)
	
	% Choose default command line output for about_box
	if (nargout),    varargout{1} = hObject;    end

% -------------------------------------------------------------------
function pushbutton_OK_Callback(hObject, eventdata, handles)
	delete(handles.figure1);

% --- Executes on key press over figure1 with no controls selected.
function figure1_KeyPressFcn(hObject, eventdata)
% Check for "escape"
handles = guidata(hObject);
if isequal(get(hObject,'CurrentKey'),'escape')
    delete(handles.figure1);
end

% --- Creates and returns a handle to the GUI figure. 
function about_box_LayoutFcn(h1)

set(h1,...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'KeyPressFcn',@figure1_KeyPressFcn,...
'MenuBar','none',...
'Name','about_box',...
'NumberTitle','off',...
'Position',[303 228 269 380],...
'Resize','off',...
'Tag','figure1');

uicontrol('Parent',h1,'Position',[13 99 241 30],...
'Max',10,...
'String',{  'M-GMT is a MATLAB GUI to the Generic Mapping'; 'Tools (GMT) software'; '' },...
'Style','text',...
'Tag','text_ProgName');

axes('Parent',h1,'Units','pixels',...
'Color',get(0,'defaultaxesColor'),...
'Position',[15 139 240 240],...
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
'FontSize',12,...
'HorizontalAlignment','left',...
'String','J Luis',...
'Style','text');

uicontrol('Parent',h1, 'Position',[138 47 113 17],...
'ForegroundColor',[0 0.501960784313725 0.250980392156863],...
'String','w3.ualg.pt/~jluis/m_gmt',...
'Style','text',...
'Tag','text_url');

uicontrol('Parent',h1,'Position',[10 91 251 1],'Style','frame','Tag','frame3');
uicontrol('Parent',h1,'Position',[10 41 251 1],'Style','frame','Tag','frame4');

uicontrol('Parent',h1, 'Position',[10 16 201 15],...
'HorizontalAlignment','left',...
'Style','text','Tag','text_prog');

uicontrol('Parent',h1, 'Position',[216 11 46 21],...
'Callback',{@about_box_uicallback,h1,'pushbutton_OK_Callback'},...
'String','OK',...
'Tag','pushbutton_OK');

function about_box_uicallback(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
feval(callback_name,hObject,[],guidata(h1));

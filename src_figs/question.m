function out = question(varargin)
% Helper window to ask for a numeric value

%	Copyright (c) 2004-2014 by J. Luis
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

% $Id: question.m 4555 2014-09-14 17:40:32Z j $

	hObject = figure('Vis','off');
	question_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hObject,'center');

	if (isempty(varargin))
		set(handles.text_desc,'Str', 'Nikles')
	else
		desc = varargin{1};
		if (isa(desc, 'cell')),		desc = desc{1};		end
		set(handles.text_desc,'Str', desc)
		if (numel(varargin) >= 2)
			set(handles.figure1,'Name', varargin{2})
			if (numel(varargin) >= 3)
				% Not used yet. Should be the edit box width
				if (numel(varargin) >= 4)
					set(handles.edit_in,'Str', varargin{4})
					if (numel(varargin) == 5)
						set(handles.check_invert, 'Vis', 'on')
					end
				end
			end
		end
	end

	handles.out = [];
	guidata(hObject, handles);

	% Make the GUI modal
	set(hObject,'Visible','on','WindowStyle','modal')

	% UIWAIT makes yes_or_no wait for user response (see UIRESUME)
	uiwait(handles.figure1);
	handles = guidata(hObject);
	out = handles.out;
	delete(hObject)

% ------------------------------------------------------------------------------
function push_OK_CB(hObject, handles)
% ...
	handles.out{1} = get(handles.edit_in, 'Str');
	handles.out{2} = get(handles.check_invert, 'Val');
	guidata(handles.figure1, handles)
	uiresume(handles.figure1);

% --- Executes on key press over figure1 with no controls selected.%
function figure1_KeyPressFcn(hObject, evt)
	if isequal(get(hObject,'CurrentKey'),'escape')
		handles = guidata(hObject);
		uiresume(handles.figure1);
	end

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, evt)
	handles = guidata(hObject);
	uiresume(handles.figure1);

% --- Creates and returns a handle to the GUI figure. 
function question_LayoutFcn(h1)

set(h1, 'Position',[520 722 200 79],...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'KeyPressFcn',@figure1_KeyPressFcn,...
'CloseRequestFcn',@figure1_CloseRequestFcn,...
'MenuBar','none',...
'Name','-?-',...
'NumberTitle','off',...
'Resize','off',...
'Tag','figure1');

uicontrol('Parent',h1, 'Position',[10 60 120 14],...
'HorizontalAlignment','left',...
'String','Text',...
'Style','text',...
'Tag','text_desc');

uicontrol('Parent',h1, 'Pos',[9 37 182 22],...
'BackgroundColor',[1 1 1],...
'HorizontalAlignment','left',...
'String','',...
'Style','edit',...
'Tag','edit_in');

uicontrol('Parent',h1, 'Position',[10 6 55 23],...
'String','Invert',...
'Style','checkbox',...
'Visible','off',...
'TooltipString','Invert the sense. That is, set the vaue OUTSIDE the polygon',...
'Tag','check_invert');

uicontrol('Parent',h1, 'Pos',[141 8 50 21],...
'Call',@question_uiCB,...
'FontSize',9,...
'FontWeight','bold',...
'String','OK',...
'Tag','push_OK');

function question_uiCB(hObject, eventdata)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));
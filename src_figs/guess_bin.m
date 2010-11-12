function varargout = guess_bin(varargin)
% Helper window to confirm/modify the guessing job of guess_file about binary files
%

%	Copyright (c) 2004-2010 by J. Luis
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

	if (nargin == 0)
		errordlg('Guess in bin ERROR: wrong number of input arguments','Error'),	return
	end

	hObject = figure('Tag','figure1','Visible','off');
	guess_bin_LayoutFcn(hObject);
	handles = guihandles(hObject);

	if (nargin == 2)			% We have a confirmation request
		set(handles.popup_nCols,'Val', varargin{1}-1)
		n = 0;
		switch varargin{2}(1:3)
			case 'sin',		n = 1;
			case 'dou',		n = 2;
			case 'int',		n = 3;
			case 'sho',		n = 4;
		end
		if (n == 0)
			errordlg('Guess in bin ERROR: wrong data type on input.','Error')
			delete(hObject),	return
		end
		set(handles.popup_type,'Val', n)
	end

	handles.output = 0;
	guidata(hObject, handles);

	set(hObject,'Visible','on');
	% UIWAIT makes guess_bin wait for user response
	uiwait(handles.figure1);

	handles = guidata(hObject);
	if (~isempty(handles.output))
		out.nCols = get(handles.popup_nCols,'Val');
		types = {'single' 'double' 'int' 'short'};
		out.type = types{get(handles.popup_type,'Val')};
		varargout{1} = out;
	else
		varargout{1} = handles.output;
	end
	delete(handles.figure1);

% ---------------------------------------------------------------------
function push_OK_CB(hObject, eventdata, handles)
	uiresume(handles.figure1);

% ---------------------------------------------------------------------
function push_Cancel_CB(hObject, eventdata, handles)
	handles.output = [];        % User gave up, return nothing
	guidata(hObject, handles);  uiresume(handles.figure1);

% --------------------------------------------------------------------------------------------------
% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata)
	handles = guidata(hObject);
	handles.output = [];        % User gave up, return nothing
	guidata(hObject, handles);    uiresume(handles.figure1);

% --- Executes on key press over figure1 with no controls selected.
function figure1_KeyPressFcn(hObject, eventdata)
	if isequal(get(hObject,'CurrentKey'),'escape')
		handles = guidata(hObject);
		handles.output = [];    % User said no by hitting escape
		guidata(hObject, handles);    uiresume(handles.figure1);
	end

% --- Creates and returns a handle to the GUI figure. 
function guess_bin_LayoutFcn(h1)

set(h1, 'Position',[520 659 228 141],...
'CloseRequestFcn',@figure1_CloseRequestFcn,...
'KeyPressFcn',@figure1_KeyPressFcn,...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'MenuBar','none',...
'Name','Guess in Bin',...
'NumberTitle','off',...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1');

uicontrol('Parent',h1, 'Position',[10 42 61 20],...
'BackgroundColor',[1 1 1],...
'String',{'2'; '3'; '4'; '5'; '6'},...
'Style','popupmenu',...
'TooltipString','Confirm number of columns in file.',...
'Value',1,...
'Tag','popup_nCols');

uicontrol('Parent',h1, 'Position',[98 42 123 21],...
'BackgroundColor',[1 1 1],...
'String',{'single (4 bytes)'; 'double (8 bytes)'; 'int (4 bytes)'; 'short int (2 bytes)'},...
'Style','popupmenu',...
'TooltipString','Confirm the data type',...
'Value',1,...
'Tag','popup_type');

uicontrol('Parent',h1, 'Position',[10 88 211 51],...
'FontSize',9,...
'String',{'This is a binary file.'; 'My guess of its contents is shown.'; 'Please confirm or correct.' },...
'Style','text');

uicontrol('Parent',h1, 'Position',[152 9 69 22],...
'Call',{@guess_bin_uiCB,h1,'push_OK_CB'},...
'FontWeight','bold',...
'String','OK',...
'Tag','push_OK');

uicontrol('Parent',h1, 'Position',[70 9 69 22],...
'Call',{@guess_bin_uiCB,h1,'push_Cancel_CB'},...
'FontWeight','bold',...
'String','Cancel',...
'Tag','push_Cancel');

uicontrol('Parent',h1, 'Position',[10 62 60 15],...
'String','N columns',...
'Style','text');

uicontrol('Parent',h1, 'Position',[129 64 60 15],...
'String','Data type',...
'Style','text');

function guess_bin_uiCB(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
	feval(callback_name,hObject,guidata(h1));

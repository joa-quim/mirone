function varargout = guess_bin(varargin)
% Helper window to confirm/modify the guessing job of guess_file about binary files
%
% Usage:
%	out = guess_bin(n_columns, data_type);
%		n_columns		Guessed number of columns
%		data_type		Guessed data type ('single' 'double' 'int' 'short int')
%
%	out = guess_bin(0);
%		Offer the default options of n_columns = 2 & data_type = 'single'
%
%	Output is a structure with
%		out.nCols		number of columns
%		out.type		Data type
%		out.twoD		Boolean telling if Z (in case it exists) is to be discarded
%	OR
%		out = [] if user hit cancel

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

% $Id$

	if (nargin == 0)
		errordlg('Guess in bin ERROR: wrong number of input arguments','Error'),	return
	end

	hObject = figure('Vis','off');
	guess_bin_LayoutFcn(hObject);
	handles = guihandles(hObject);

	if (nargin == 2)			% We have a confirmation request
		set(handles.popup_nCols,'Val', varargin{1}-1)
		if (varargin{1} > 2)
			set(handles.check_2D,'Vis','on')
		end
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
		out.nCols = get(handles.popup_nCols,'Val') + 1;
		types = {'single' 'double' 'int' 'short'};
		out.type = types{get(handles.popup_type,'Val')};
		out.twoD = ((get(handles.popup_nCols,'Val') > 1) && get(handles.check_2D,'Val'));
		out.multiseg = get(handles.check_multiseg,'Val');
		varargout{1} = out;
	else
		varargout{1} = handles.output;
	end
	delete(handles.figure1);

% ---------------------------------------------------------------------
function popup_nCols(hObject, handles)
	if (get(hObject,'Val') > 1),	set(handles.check_2D,'Vis','on')
	else							set(handles.check_2D,'Vis','off')
	end

% ---------------------------------------------------------------------
function push_OK_CB(hObject, handles)
	uiresume(handles.figure1);

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

set(h1, 'Position',[500 420 228 151],...
'CloseRequestFcn',@figure1_CloseRequestFcn,...
'KeyPressFcn',@figure1_KeyPressFcn,...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'MenuBar','none',...
'Name','Guess in Bin',...
'NumberTitle','off',...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1');

uicontrol('Parent',h1, 'Position',[10 52 61 21],...
'Call',@guess_bin_uiCB,...
'BackgroundColor',[1 1 1],...
'String',{'2'; '3'; '4'; '5'; '6'},...
'Style','popupmenu',...
'TooltipString','Confirm number of columns in file.',...
'Value',1,...
'Tag','popup_nCols');

uicontrol('Parent',h1, 'Position',[98 52 123 21],...
'BackgroundColor',[1 1 1],...
'String',{'single (4 bytes)'; 'double (8 bytes)'; 'int (4 bytes)'; 'short int (2 bytes)'},...
'Style','popupmenu',...
'Tooltip','Confirm the data type',...
'Value',1,...
'Tag','popup_type');

uicontrol('Parent',h1, 'Position',[10 95 211 51],...
'FontSize',9,...
'String',{'This is a binary file.'; 'My guess of its contents is shown.'; 'Please confirm or correct.' },...
'Style','text');

uicontrol('Parent',h1, 'Position',[10  74 60 15],'String','N columns','Style','text');
uicontrol('Parent',h1, 'Position',[129 74 60 15],'String','Data type','Style','text');

uicontrol('Parent',h1, 'Position',[10 28 90 21],...
'Style','checkbox',...
'String','Discard Z',...
'Tooltip','Do not store the Z. Lighter to render but will not do 3D in Fledermaus',...
'Value',1,...
'Visible','off',...
'Tag','check_2D');

uicontrol('Parent',h1, 'Position',[10 8 90 21],...
'Style','checkbox',...
'String','Split multiseg',...
'Tooltip','If file is multisegment plot separate lines, one for each segment.',...
'Value',0,...
'Visible','on',...
'Tag','check_multiseg');

uicontrol('Parent',h1, 'Position',[161 9 60 22],...
'Call',@guess_bin_uiCB,...
'FontWeight','bold',...
'String','OK',...
'Tag','push_OK');

function guess_bin_uiCB(hObject, eventdata)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));

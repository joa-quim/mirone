function varargout = multiscale(varargin) 
% Auxiliary figure to interface with the mirblock.c MEX

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
	multiscale_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hObject,'right')
 
	handles.nWin = 3;		% Default to a 3x3 window
	handles.method = 0;		% "Terrain Ruggedness Index"
	handles.method_name = 'Terrain Ruggedness Index';

	% Update handles structure
	guidata(hObject, handles);

	% Make the GUI modal
	set(handles.figure1,'WindowStyle','modal')
	set(hObject,'Visible','on');

	% UIWAIT makes yes_or_no wait for user response (see UIRESUME)
	uiwait(handles.figure1);
	handles = guidata(hObject);
	varargout{1} = [];
	if (handles.ok)
		varargout{1} = struct('method',handles.method, 'size',handles.nWin, 'name',handles.method_name);
	end
	
	delete(handles.figure1)

% ----------------------------------------------------------------------------
function popup_algo_CB(hObject, handles)
% VERY IMPORTANT. The Method's order in popup must agree strictly with the mirblock.c MEX
	val = get(hObject,'Value');		str = get(hObject,'String');
	handles.method = val - 1;
	handles.method_name = str{val};
	guidata(handles.figure1,handles)

% ----------------------------------------------------------------------------
function edit_nWin_CB(hObject, handles)
	xx = round( str2double(get(hObject,'String')) );
	if (isnan(xx)),		set(hObject,'String',3),	return,		end
	if (rem(xx,2) == 0)
		errordlg('Window dimension must be an odd number. Fixing it','Error')
		xx = xx + 1;
	end
	handles.nWin = xx;
	if (handles.nWin < 3)       % Minimum allowed is 3
		set(hObject,'String',3)
		handles.nWin = 3;
	end
	guidata(handles.figure1,handles)

% ----------------------------------------------------------------------------
function push_OK_CB(hObject, handles)
	handles.ok = true;
	guidata(hObject, handles);
	% Use UIRESUME instead of delete because the OutputFcn needs to get the updated handles structure.
	uiresume(handles.figure1);	

% ----------------------------------------------------------------------------
function push_cancel_CB(hObject, handles)
	handles.ok = false;
	guidata(hObject, handles);
	% Use UIRESUME instead of delete because the OutputFcn needs to get the updated handles structure.
	uiresume(handles.figure1);

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata)
	handles = guidata(hObject);
	if (exist('OCTAVE_VERSION','builtin'))		% To know if we are running under Octave
		do_uiresume = ( isprop(hObject, '__uiwait_state__') && strcmp(get(hObject, '__uiwait_state__'), 'active') );
	else
		do_uiresume = strcmp(get(handles.figure1, 'waitstatus'), 'waiting');
	end
	if (do_uiresume)		% The GUI is still in UIWAIT, us UIRESUME
		handles.ok = false;
		guidata(hObject, handles);	uiresume(handles.figure1);
	else					% The GUI is no longer waiting, just close it
		delete(handles.figure1);
	end

% --- Executes on key press over figure1 with no controls selected.
function figure1_KeyPressFcn(hObject, eventdata)
	handles = guidata(hObject);
	if isequal(get(hObject,'CurrentKey'),'escape')
		handles.ok = false;
        guidata(hObject, handles);
        uiresume(handles.figure1);
	end    


% --- Creates and returns a handle to the GUI figure. 
function multiscale_LayoutFcn(h1)

set(h1, 'Position',[520 670 380 120],...
'CloseRequestFcn',@figure1_CloseRequestFcn,...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'KeyPressFcn',@figure1_KeyPressFcn,...
'MenuBar','none',...
'Name','Multiscale',...
'NumberTitle','off',...
'RendererMode','manual',...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@multiscale_uiCB,...
'Position',[60 49 181 22],...
'String',{'Terrain Ruggedness Index'; 'Topographic Position Index'; 'Roughness (aka Range)'; ...
	'Mean'; 'Min'; 'Max'; 'Slope'; 'Aspect'; 'RMS'; 'Trend'; 'Residue'; 'RMS of Residue'; ...
	'AGC (Full Amp)'; 'AGC (Local Amp)' },...
'Style','popupmenu',...
'TooltipString','Select one method',...
'Value',1,...
'Tag','popup_algo');

uicontrol('Parent',h1, 'Position',[339 50 31 21],...
'BackgroundColor',[1 1 1],...
'Call',@multiscale_uiCB,...
'String','3',...
'Style','edit',...
'Tooltip','Width of the rectangular neighborhood (MUST bo an odd number)',...
'Tag','edit_nWin');

uicontrol('Parent',h1, 'Position',[7 54 51 15],...
'FontName','Helvetica', 'FontSize',9,...
'HorizontalAlignment','right', 'String','Method',...
'Style','text');

uicontrol('Parent',h1, 'Position',[260 52 78 17],...
'FontName','Helvetica', 'FontSize',9,...
'HorizontalAlignment','right', 'String','Window size',...
'Style','text',...
'TooltipString','Width of the rectangular neighborhood (MUST bo an odd number)');

uicontrol('Parent',h1, 'Position',[22 87 340 22],...
'FontAngle','italic',...
'FontName','Helvetica', 'FontSize',12, 'FontWeight','bold',...
'String','Multi-scale methods for Terrain Analysis',...
'Style','text');

uicontrol('Parent',h1, 'Position',[225 9 66 21],...
'Call',@multiscale_uiCB,...
'FontName','Helvetica', 'FontSize',9,...
'String','OK', 'Tag','push_OK');

uicontrol('Parent',h1, 'Position',[305 9 66 21],...
'Call',@multiscale_uiCB,...
'FontName','Helvetica',...
'FontSize',9,...
'String','Cancel',...
'Tag','push_cancel');

function multiscale_uiCB(hObject, eventdata)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));

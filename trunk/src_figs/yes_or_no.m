function varargout = yes_or_no(varargin)
%   Just ask a yes or no question.
%   By default (h=yes_or_no('title','Warning')) the question is the 'grid_max_size ... sand ... bla bla'
%   The form h=yes_or_no('string',str,'title','Warning') returns the answer
%   to the question contained in the cellstring STR

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

if (nargin == 4 && isempty(varargin{3}))
    % Come here when called by one of the buttons callback
    gui_Callback = str2func(varargin{1});
    feval(gui_Callback,varargin{2:end})
else
    h = yes_or_no_OpeningFcn(varargin{:});
    if (nargout)    varargout{1} = h;   end
end


% --- Executes just before yes_or_no is made visible.
function output = yes_or_no_OpeningFcn(varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to yes_or_no (see VARARGIN)

hObject = yes_or_no_LayoutFcn;
handles = guihandles(hObject);

% Choose default command line output for yes_or_no
handles.output = 'Yes';

% Insert custom Title and Text if specified by the user
if(nargin > 1)
    for (index = 1:2:(nargin-1))
        switch lower(varargin{index})
            case 'title'
                set(hObject, 'Name', varargin{index+1});
            case 'string'
                set(handles.text1, 'String', varargin{index+1});
        end
    end
end

% Determine the position of the dialog - centered on the callback figure
% if available, else, centered on the screen
FigPos=get(0,'DefaultFigurePosition');
OldUnits = get(hObject, 'Units');
set(hObject, 'Units', 'pixels');
OldPos = get(hObject,'Position');
FigWidth = OldPos(3);
FigHeight = OldPos(4);
ScreenUnits=get(0,'Units');
set(0,'Units','pixels');
ScreenSize=get(0,'ScreenSize');
set(0,'Units',ScreenUnits);

if isempty(gcbf)
    FigPos(1)=1/2*(ScreenSize(3)-FigWidth);
    FigPos(2)=2/3*(ScreenSize(4)-FigHeight);
else
    GCBFOldUnits = get(gcbf,'Units');
    set(gcbf,'Units','pixels');
    GCBFPos = get(gcbf,'Position');
    set(gcbf,'Units',GCBFOldUnits);
    FigPos(1:2) = [(GCBFPos(1) + GCBFPos(3) / 2) - FigWidth / 2, ...
                   (GCBFPos(2) + GCBFPos(4) / 2) - FigHeight / 2];
end
FigPos(3:4)=[FigWidth FigHeight];
if (FigPos(2)+FigPos(4) > ScreenSize(4))            % Figure is actually partially outside screen
    FigPos(2) = ScreenSize(4) - FigPos(4) - 30;     % The 30 is to account for figure's blue bar
end
set(hObject, 'Position', FigPos, 'Units', OldUnits, 'Visible', 'on');
guidata(hObject, handles);

% Make the GUI modal
set(handles.figure1,'WindowStyle','modal')

% UIWAIT makes yes_or_no wait for user response (see UIRESUME)
uiwait(handles.figure1);
handles = guidata(hObject);
output = yes_or_no_OutputFcn(hObject, [], handles);

% --- Outputs from this function are returned to the command line.
function varargout = yes_or_no_OutputFcn(hObject, eventdata, handles)
% Get default command line output from handles structure
varargout{1} = handles.output;
% The figure can be deleted now
delete(handles.figure1);

% ----------------------------------------------------------------------
function pushbutton1_Callback(hObject, eventdata, handles)
	handles.output = get(hObject,'String');
	guidata(hObject, handles);
	% Use UIRESUME instead of delete because the OutputFcn needs to get the updated handles structure.
	uiresume(handles.figure1);

% ----------------------------------------------------------------------
function pushbutton2_Callback(hObject, eventdata, handles)
	handles.output = get(hObject,'String');
	guidata(hObject, handles);
	% Use UIRESUME instead of delete because the OutputFcn needs to get the updated handles structure.
	uiresume(handles.figure1);

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
	handles = guidata(hObject);
	if (exist('OCTAVE_VERSION','builtin'))		% To know if we are running under Octave
		do_uiresume = ( isprop(hObject, '__uiwait_state__') && strcmp(get(hObject, '__uiwait_state__'), 'active') );
	else
		do_uiresume = strcmp(get(handles.figure1, 'waitstatus'), 'waiting');
	end
	if (do_uiresume)		% The GUI is still in UIWAIT, us UIRESUME
		uiresume(handles.figure1);
	else					% The GUI is no longer waiting, just close it
		delete(handles.figure1);
	end

% --- Executes on key press over figure1 with no controls selected.
function figure1_KeyPressFcn(hObject, eventdata, handles)
	% Check for "enter" or "escape"
	if isequal(get(hObject,'CurrentKey'),'escape')
        % User said Yes (the default answer) by hitting escape
        handles.output = 'Yes';
        guidata(hObject, handles);
        uiresume(handles.figure1);
	end    
	if isequal(get(hObject,'CurrentKey'),'return')
        uiresume(handles.figure1);
	end    


% --- Creates and returns a handle to the GUI figure. 
function h1 = yes_or_no_LayoutFcn()
h1 = figure(...
'CloseRequestFcn','yes_or_no(''figure1_CloseRequestFcn'',gcf,[],guidata(gcf))',...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'IntegerHandle','off',...
'KeyPressFcn','yes_or_no(''figure1_KeyPressFcn'',gcbo,[],guidata(gcbo))',...
'MenuBar','none',...
'Name','yes_or_no',...
'NumberTitle','off',...
'Position',[657 484 360 165],...
'Resize','off',...
'Visible', 'off',...
'HandleVisibility','callback',...
'Tag','figure1');

uicontrol('Parent',h1,...
'Callback','yes_or_no(''pushbutton1_Callback'',gcbo,[],guidata(gcbo))',...
'FontWeight','bold',...
'Position',[93 8 81 25],...
'String','Yes','Tag','pushbutton1');

uicontrol('Parent',h1,...
'Callback','yes_or_no(''pushbutton2_Callback'',gcbo,[],guidata(gcbo))',...
'Position',[234 8 81 25],...
'String','No','Tag','pushbutton2');

uicontrol('Parent',h1,'FontSize',10,...
'HorizontalAlignment','left',...
'Position',[11 45 341 112],...
'FontUnits','points',...
'FontSize',10.0,...
'String',{'According to your own declaration of grid_max_size '; '(Yes you did declare it, even if you don''t know it because'; ' "I don''t like Manuals"),'; 'this grid is too big for your teeth. In these cases it is'; 'VERY WISE to STOP here and use the "Overview tool".'; ''; 'Shell I stop now?' },...
'Style','text','Tag','text1');

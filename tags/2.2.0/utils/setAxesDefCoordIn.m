function handles = setAxesDefCoordIn(handles, opt)
% Sets the value of the axes uicontextmenu that selects what will be donne
% when Loading a file in terms of needing, or not, to project it.

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

    if (nargin == 1),   handles = guidata(handles.figure1);     end     % Otherwise uses input handles

    if (handles.is_projected)
        if (~handles.geog)          % Is projected
            set(handles.hAxMenuLF, 'Label', 'Load files in geogs', ...
                'Call', {@CoordMode_CB,handles.figure1}, 'Vis', 'on', 'Separator','on');
            handles.defCoordsIn = 1;
        else                        % Hum, probably a GDAL file in geogs
            set(handles.hAxMenuLF, 'Vis', 'off', 'Separator','off');
            handles.defCoordsIn = -1;
        end
    else                                % Otherwise remove eventual previous one Load in projected coords Load files in geogs
        set(handles.hAxMenuLF, 'Vis', 'off', 'Separator','off');    % No projection known, no choices offered
        handles.defCoordsIn = 0;        % We have no idea about the coordinates type
    end
    
    if (handles.defCoordsIn == 1)
        set(handles.hAxMenuLF, 'Vis', 'on', 'Separator','on');
        set(handles.hAxMenuDM, 'Call', {@DisplayMode_CB,handles.figure1}, 'Vis', 'on', 'Separator','on');
    else
        set(handles.hAxMenuLF, 'Vis', 'off', 'Separator','off');    % No projection known, no choices offered
        set(handles.hAxMenuDM, 'Vis', 'off', 'Separator','off');
    end
    
    if (nargout == 0),    guidata(handles.figure1, handles);    end

% --------------------------------------------------------------------
function CoordMode_CB(hObject, event, hFig)
	% Swapp Labels which inform of the current default mode of interpreting the coords of loaded files
    % Also uppdates the handles.defCoords_type field.
    handles = guidata(hFig);        % We need the most updated version
	if (strncmp(get(hObject,'Label'),'Load in',7))          % Was projected, now geogs
        set(hObject,'Label','Load files in geogs');         handles.defCoordsIn = 1;
	else
        set(hObject,'Label','Load in projected coords');    handles.defCoordsIn = -1;
	end
    guidata(handles.figure1, handles)

% --------------------------------------------------------------------
function DisplayMode_CB(hObject, event, hFig)
	if (strncmp(get(hObject,'Label'),'Display projected',17))   % Was projected, now geogs
        set(hObject,'Label','Display coords in geogs');
        setappdata(hFig,'DispInGeogs',1)
	else
        set(hObject,'Label','Display projected coords');
        setappdata(hFig,'DispInGeogs',0)
	end

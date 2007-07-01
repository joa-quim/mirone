function handles = setAxesDefCoordIn(handles)
    % Sets the value of the axes uicontextmenu that selects what will be donne
    % when Loading a file in terms of needing, or not, to project it.

    % Fish eventual proj strings
    projGMT = getappdata(handles.figure1,'ProjGMT');
    projWKT = getappdata(handles.axes1,'ProjWKT');

    if (~isempty(projWKT) || ~isempty(projGMT))
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

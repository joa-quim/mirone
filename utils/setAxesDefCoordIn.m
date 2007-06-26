function setAxesDefCoordIn(handles)
    % Sets the value of the axes uicontextmenu that selects what will be donne
    % when Loading a file in terms of needing, or not, to project it.

    % Fish eventual proj strings
    projGMT = getappdata(handles.figure1,'ProjGMT');
    projWKT = getappdata(handles.axes1,'ProjWKT');

    if (~isempty(projWKT) || ~isempty(projGMT))
        if (~handles.geog)
            set(handles.hAxMenuLF, 'Label', 'Load in projected coords', ...
                'Call', {@CoordMode_CB,handles}, 'Vis', 'on', 'Separator','on');
            handles.defCoordsIn = -1;
        else
            set(handles.hAxMenuLF, 'Vis', 'off', 'Separator','off');
            handles.defCoordsIn = 1;
        end
    else                                % Otherwise remove eventual previous one Load in projected coords Load files in geogs
        set(handles.hAxMenuLF, 'Vis', 'off', 'Separator','off');    % No projection known, no choices offered
        handles.defCoordsIn = 0;        % We have no idea about the coordinates type
    end
    guidata(handles.figure1, handles)

% --------------------------------------------------------------------
function CoordMode_CB(hObject, event, handles)
	% Swapp Labels which inform of the current default mode of interpreting the coords of loaded files
    % Also uppdates the handles.defCoords_type field.
	if (strncmp(get(hObject,'Label'),'Load in',7))          % Was projected, now geogs
        set(hObject,'Label','Load files in geogs');         handles.defCoordsIn = 1;
	else
        set(hObject,'Label','Load in projected coords');    handles.defCoordsIn = -1;
	end
    guidata(handles.figure1, handles)

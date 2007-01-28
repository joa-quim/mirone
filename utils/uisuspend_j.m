function uistate = uisuspend_j(fig, setdefaults)
%UISUSPEND suspends all interactive properties of the figure.
%
%   UISTATE=UISUSPEND(FIG) suspends the interactive properties of a 
%   figure window and returns the previous state in the structure
%   UISTATE.  This structure contains information about the figure's
%   WindowButton* functions and the pointer.  It also contains the 
%   ButtonDownFcn's for all children of the figure.
%
%   UISTATE=UISUSPEND(FIG,FALSE) returns the structure as above but leaves
%   the current settings unchanged.
%   
%   See also UIRESTORE.

%   Copyright 1984-2004 The MathWorks, Inc.
%   $Revision: 1.10.4.4 $ $Date: 2004/08/16 01:54:55 $
%
%   UISUSPEND of R13 is very bugy. This is a replacement from R14 but with
%   'plotedit' desabled (We don't use that in Mirone, well for the time beeing ...)

if (nargin < 2)
    setdefaults = logical(1);
end

%chi = findobj(fig);
chi = findobj(fig,'Type','axes');               % Megas saving (and speed)
chi = [chi; findobj(fig,'Type','image')];
chi = [chi; findobj(fig,'Type','line')];
chi = [chi; findobj(fig,'Type','patch')];
chi = [chi; findobj(fig,'Type','text')];

uistate = struct(...
        'figureHandle',          fig, ...
        'Children',              chi, ...
        'WindowButtonMotionFcn', Lwrap(get(fig, 'WindowButtonMotionFcn')), ...
        'WindowButtonDownFcn',   Lwrap(get(fig, 'WindowButtonDownFcn')), ...
        'WindowButtonUpFcn',     Lwrap(get(fig, 'WindowButtonUpFcn')), ...
        'KeyPressFcn',           Lwrap(get(fig, 'KeyPressFcn')), ...
        'Pointer',               get(fig, 'Pointer'), ...
        'PointerShapeCData',     get(fig, 'PointerShapeCData'), ...
        'PointerShapeHotSpot',   get(fig, 'PointerShapeHotSpot'), ...
        'ButtonDownFcns',        Lwrap(get(chi, {'ButtonDownFcn'})), ...
        'Interruptible',         Lwrap(get(chi, {'Interruptible'})), ...
        'BusyAction',            Lwrap(get(chi, {'BusyAction'})), ...
        'UIContextMenu',         Lwrap(get(chi, {'UIContextMenu'})) );

if (setdefaults)
    % disable plot editing and annotation buttons
    %plotedit(fig,'setenabletools','off'); % ploteditEnable      COMENTED (JL)
    % do nothing figureHandle
    % do nothing for Children
    set(fig, 'WindowButtonMotionFcn', get(0, 'DefaultFigureWindowButtonMotionFcn'))
    set(fig, 'WindowButtonDownFcn',   get(0, 'DefaultFigureWindowButtonDownFcn'))
    set(fig, 'WindowButtonUpFcn',     get(0, 'DefaultFigureWindowButtonUpFcn'))
    set(fig, 'KeyPressFcn',           get(0, 'DefaultFigureKeyPressFcn'))
    set(fig, 'Pointer',               get(0, 'DefaultFigurePointer'))
    set(fig, 'PointerShapeCData',     get(0, 'DefaultFigurePointerShapeCData'))
    set(fig, 'PointerShapeHotSpot',   get(0, 'DefaultFigurePointerShapeHotSpot'))
    set(chi, 'ButtonDownFcn',         '')
    set(chi, 'Interruptible',         'on');
    set(chi, 'BusyAction',            'Queue')
    % do nothing for UIContextMenu
end

% wrap cell arrays in another cell array for passing to the struct command
function x = Lwrap(x)
if (iscell(x)),  x = {x};    end

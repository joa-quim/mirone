function uistate = uisuspend_fig(fig, setdefaults)
%UISUSPEND_FIG very minimalist version of UISUSPEND that suspends only fig properties.
%
%   UISTATE=UISUSPEND_FIG(FIG) suspends the interactive properties of a 
%   figure window and returns the previous state in the structure
%   UISTATE.  This structure contains information about the figure's
%   WindowButton* functions and the pointer.  
%
%   UISTATE=UISUSPEND_FIG(FIG,FALSE) returns the structure as above but leaves
%   the current settings unchanged.
%
%  Joaquim Luis 25-01-07

if (nargin < 2)
    setdefaults = logical(1);
end

uistate = struct(...
        'figureHandle',          fig, ...
        'WindowButtonMotionFcn', Lwrap(get(fig, 'WindowButtonMotionFcn')), ...
        'WindowButtonDownFcn',   Lwrap(get(fig, 'WindowButtonDownFcn')), ...
        'WindowButtonUpFcn',     Lwrap(get(fig, 'WindowButtonUpFcn')), ...
        'KeyPressFcn',           Lwrap(get(fig, 'KeyPressFcn')), ...
        'Pointer',               get(fig, 'Pointer'), ...
        'PointerShapeCData',     get(fig, 'PointerShapeCData'), ...
        'PointerShapeHotSpot',   get(fig, 'PointerShapeHotSpot') );

if (setdefaults)
    set(fig, 'WindowButtonMotionFcn', get(0, 'DefaultFigureWindowButtonMotionFcn'))
    set(fig, 'WindowButtonDownFcn',   get(0, 'DefaultFigureWindowButtonDownFcn'))
    set(fig, 'WindowButtonUpFcn',     get(0, 'DefaultFigureWindowButtonUpFcn'))
    set(fig, 'KeyPressFcn',           get(0, 'DefaultFigureKeyPressFcn'))
    set(fig, 'Pointer',               get(0, 'DefaultFigurePointer'))
    set(fig, 'PointerShapeCData',     get(0, 'DefaultFigurePointerShapeCData'))
    set(fig, 'PointerShapeHotSpot',   get(0, 'DefaultFigurePointerShapeHotSpot'))
end

% wrap cell arrays in another cell array for passing to the struct command
function x = Lwrap(x)
if (iscell(x)),  x = {x};    end

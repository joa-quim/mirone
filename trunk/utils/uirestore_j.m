function uirestore_j(uistate,kidsOnly)
%UIRESTORE Restores the interactive functionality figure window.
%  UIRESTORE(UISTATE) restores the state of a figure window to
%  what it was before it was suspended by a call to UISUSPEND.
%  The input UISTATE is the structure returned by UISUSPEND which
%  contains information about the window properties and the button
%  down functions for all of the objects in the figure.
%
%  UIRESTORE(UISTATE, 'children') updates ONLY the children of
%  the figure.
%
%  UIRESTORE(UISTATE, 'nochildren') updates ONLY the figure.
%
%  UIRESTORE(UISTATE, 'uicontrols') updates ONLY the uicontrol children
%  of the figure
%
%  UIRESTORE(UISTATE, 'nouicontrols') updates the figure and all non-uicontrol
%  children of the figure
%
%  See also UISUSPEND.

%   Copyright 1984-2005 The MathWorks, Inc.
%   $Revision: 1.16.4.4 $ $Date: 2005/06/21 19:43:40 $
%
%   UISUSPEND/RESTORE of R13 is very bugy. This is a replacement from R14 but with
%   'plotedit' desabled (We don't use that in Mirone, well for the time beeing ...)

fig = uistate.figureHandle;

% No need to restore anything if the figure isn't there
if (~ishandle(fig)),    return;     end

% No need to restore if the figure is being deleted
if strcmp('on',get(fig,'BeingDeleted'))
    return
end

updateFigure   = 1;
updateChildren = 1;
updateUICtrl   = 1;

if (nargin == 2)
    if strcmpi(kidsOnly, 'children');
        updateFigure   = 0;
    elseif strcmpi(kidsOnly, 'nochildren');
        updateChildren = 0;
        updateUICtrl   = 0;
    elseif strcmpi(kidsOnly,'uicontrols')
        updateChildren = 0;
    elseif strcmpi(kidsOnly,'nouicontrols')
        updateUICtrl   = 0;
    end
end

% if ~isempty(uistate.ploteditEnable)
%     plotedit(fig, 'setenabletools',uistate.ploteditEnable);
% end

if (updateFigure)
    set(fig, 'WindowButtonMotionFcn', uistate.WindowButtonMotionFcn)
    set(fig, 'WindowButtonDownFcn',   uistate.WindowButtonDownFcn)
    set(fig, 'WindowButtonUpFcn',     uistate.WindowButtonUpFcn);
    set(fig, 'KeyPressFcn',           uistate.KeyPressFcn)
    set(fig, 'Pointer',               uistate.Pointer)
    set(fig, 'PointerShapeCData',     uistate.PointerShapeCData)
    set(fig, 'PointerShapeHotSpot',   uistate.PointerShapeHotSpot)

    if (isfield(uistate,'docontext') && uistate.docontext)
        % Do not set an invalid handle or a handle that belongs to another figure. g239172
        if isempty(uistate.UIContextMenu{1}) || ( ishandle(uistate.UIContextMenu{1}) && ...
              isequal(fig,ancestor(uistate.UIContextMenu{1},'figure') ))
            set(fig,'UIContextMenu', uistate.UIContextMenu{1});
        end
    end
end

if (updateChildren && updateUICtrl)     % updating children including uicontrols
    LupdateChildren(uistate, [], []);
elseif (updateChildren)                 % updating non-uicontol children only
    LupdateChildren(uistate, 'uicontrol', false);
elseif (updateUICtrl)                   % updating only uicontrol children
    LupdateChildren(uistate, 'uicontrol', true);
end

%----------------------------------------------------
function LupdateChildren(uistate, childType, include)
    fig = [];
    for i=1:length(uistate.Children)
        chi = uistate.Children(i);
        if (~ishandle(chi)),    continue;   end
        if isempty(fig) || ~ishandle(fig)
            fig = ancestor(chi,'figure');
        end
        if ~isempty(childType)
            type = get(chi,'type');
            same = strcmp(type, childType);
            if (include && ~same),   continue;   end
            if (~include && same),   continue;   end
        end
        set(chi, {'ButtonDownFcn'},  uistate.ButtonDownFcns(i))
        set(chi, {'Interruptible'},  uistate.Interruptible(i))
        set(chi, {'BusyAction'},     uistate.BusyAction(i))
        if isfield(uistate,'docontext') && uistate.docontext
            % Do not set an invalid handle or a handle that belongs to another figure. g239172
            if isempty(uistate.UIContextMenu{i}) || ( ishandle(uistate.UIContextMenu{i}) && ...
                 isequal(fig,ancestor(uistate.UIContextMenu{i},'figure')))
                    set(chi, {'UIContextMenu'}, uistate.UIContextMenu(i));
            end
        end
    end

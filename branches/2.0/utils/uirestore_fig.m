function uirestore_fig(uistate)
%UIRESTORE_FIG Restores the interactive functionality figure window.
%
%  UIRESTORE_FIG(UISTATE) restores the state of a figure window to
%  what it was before it was suspended by a call to UISUSPEND_FIG.
%  The input UISTATE is the structure returned by UISUSPEND_FIG which
%  contains information about the figure's WindowButton* functions and the pointer
%
%  Joaquim Luis 25-01-07

fig = uistate.figureHandle;

% No need to restore anything if the figure isn't there
if (~ishandle(fig)),    return;     end

% No need to restore if the figure is being deleted
if strcmp('on',get(fig,'BeingDeleted'))
    return
end

set(fig, 'WindowButtonMotionFcn', uistate.WindowButtonMotionFcn)
set(fig, 'WindowButtonDownFcn',   uistate.WindowButtonDownFcn)
set(fig, 'WindowButtonUpFcn',     uistate.WindowButtonUpFcn);
set(fig, 'KeyPressFcn',           uistate.KeyPressFcn)
set(fig, 'Pointer',               uistate.Pointer)
set(fig, 'PointerShapeCData',     uistate.PointerShapeCData)
set(fig, 'PointerShapeHotSpot',   uistate.PointerShapeHotSpot)

%     if (isfield(uistate,'docontext') && uistate.docontext)
%         % Do not set an invalid handle or a handle that belongs to another figure. g239172
%         if isempty(uistate.UIContextMenu{1}) || ( ishandle(uistate.UIContextMenu{1}) && ...
%               isequal(fig,ancestor(uistate.UIContextMenu{1},'figure') ))
%             set(fig,'UIContextMenu', uistate.UIContextMenu{1});
%         end
%     end

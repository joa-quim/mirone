function uirestore_fig(uistate)
%UIRESTORE_FIG Restores the interactive functionality figure window.
%
%  UIRESTORE_FIG(UISTATE) restores the state of a figure window to
%  what it was before it was suspended by a call to UISUSPEND_FIG.
%  The input UISTATE is the structure returned by UISUSPEND_FIG which
%  contains information about the figure's WindowButton* functions and the pointer

%	Copyright (c) 2004-2011 by J. Luis
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

	fig = uistate.figureHandle;

	% No need to restore anything if the figure isn't there
	if (~ishandle(fig)),	return,		end

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

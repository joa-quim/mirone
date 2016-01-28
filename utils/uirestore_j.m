function uirestore_j(uistate,kidsOnly)
%UIRESTORE_J Restores the interactive functionality figure window.
%  UIRESTORE_J(UISTATE) restores the state of a figure window to
%  what it was before it was suspended by a call to UISUSPEND_J.
%  The input UISTATE is the structure returned by UISUSPEND_J which
%  contains information about the window properties and the button
%  down functions for all of the objects in the figure.
%
%  UIRESTORE_J(UISTATE, 'nochildren') updates ONLY the figure.
%  UIRESTORE_J(UISTATE, 'uicontrols') updates ONLY the uicontrol children of the figure
%  UIRESTORE_J(UISTATE, 'nouicontrols') updates the figure and all non-uicontrol figure's children

%   UISUSPEND/RESTORE of R13 is very bugy. This is a minimalist replacement which only
%	resets properties that matter for Mirone usage.
  
%   Coffeeright J. Luis 2004-2012

% $Id: uirestore_j.m 3624 2012-07-22 20:17:57Z j $

	hFig = uistate.figureHandle;

	if (~ishandle(hFig)),	return,		end

	updateChildren = true;
	updateUICtrl   = true;

	if (nargin == 2)
		if strcmpi(kidsOnly, 'nochildren');
			updateChildren = false;
			updateUICtrl   = false;
		elseif strcmpi(kidsOnly,'uicontrols')
			updateChildren = false;
		elseif strcmpi(kidsOnly,'nouicontrols')
			updateUICtrl   = false;
		end
	end

	set(hFig, 'WindowButtonMotionFcn', uistate.WindowButtonMotionFcn)
	set(hFig, 'WindowButtonDownFcn',   uistate.WindowButtonDownFcn)
	set(hFig, 'WindowButtonUpFcn',     uistate.WindowButtonUpFcn)
	set(hFig, 'KeyPressFcn',           uistate.KeyPressFcn)
	set(hFig, 'Pointer',               uistate.Pointer)
	set(hFig, 'PointerShapeCData',     uistate.PointerShapeCData)
	set(hFig, 'PointerShapeHotSpot',   uistate.PointerShapeHotSpot)

	if (updateChildren && updateUICtrl)		% updating children including uicontrols
		restoreChildren(uistate, []);
	elseif (updateChildren)					% updating non-uicontol children only
		restoreChildren(uistate, false);
	elseif (updateUICtrl)					% updating only uicontrol children
		restoreChildren(uistate, true);
	end

%-------------------------------------------------------------------------------
function restoreChildren(uistate, include)
	for i = 1:numel(uistate.Children)
		if (~ishandle(uistate.Children(i))),	continue,	end
		if ~isempty(include)
			same = strcmp(get(uistate.Children(i),'type'), 'uicontrol');
			if (include && ~same),		continue,	end		% TRUE => This child IS NOT a 'uicontrol'
			if (~include && same),		continue,	end
		end
		set(uistate.Children(i), {'ButtonDownFcn'},  uistate.ButtonDownFcns(i))
		set(uistate.Children(i), {'Interruptible'},  uistate.Interruptible(i))
		set(uistate.Children(i), {'BusyAction'},     uistate.BusyAction(i))
	end

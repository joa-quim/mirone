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

updateFigure   = true;
updateChildren = true;
updateUICtrl   = true;

if (nargin == 2)
    if strcmpi(kidsOnly, 'children');
        updateFigure   = false;
    elseif strcmpi(kidsOnly, 'nochildren');
        updateChildren = false;
        updateUICtrl   = false;
    elseif strcmpi(kidsOnly,'uicontrols')
        updateChildren = false;
    elseif strcmpi(kidsOnly,'nouicontrols')
        updateUICtrl   = false;
    end
end

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
              isequal(fig,ancestor_m(uistate.UIContextMenu{1},'figure') ))
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
			fig = ancestor_m(chi,'figure');
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
				 isequal(fig,ancestor_m(uistate.UIContextMenu{i},'figure')))
					set(chi, {'UIContextMenu'}, uistate.UIContextMenu(i));
			end
		end
	end

% ---------------------------------------------------------------------------
function p = ancestor_m(p,type,varargin)
%ANCESTOR  Get object ancestor.
%    P = ANCESTOR(H,TYPE) returns the handle of the closest ancestor
%    of h of one of the types in TYPE, or empty if none exists. TYPE
%    may be a single string (single type) or cell array of strings
%    (types). If H is a vector of handles the P is a cell array of the
%    same length as H and P{n} is the ancestor of H(n). If H has one
%    of the specified types then the ancestor of H is H itself.
%    P = ANCESTOR(H,TYPE,'TOPLEVEL') finds the highest level ancestor of
%    one of the types in TYPE
%
%    If H is not an Handle Graphics object, ANCESTOR returns empty.
%
%  Examples:
%    p = ancestor(gca,'figure');
%    p = ancestor(gco,{'hgtransform','hggroup','axes'},'toplevel');

%   Copyright 1984-2004 The MathWorks, Inc.
%   $Revision: 1.1.6.9 $  $Date: 2004/08/16 01:47:07 $

if (numel(p) > 1)
    n = numel(p);
    pv = cell(n,1);
    for k=1:n
        pv{k} = ancestor_m(p(k),type,varargin{:});
    end
    p = pv;
    return
end

if (nargin == 2)        % ancestor(h,type)
    while ~isempty(p) && ~isatype(p,type)
        p = get(handle_j(p),'parent');
    end
elseif nargin==3 % ancestor(h,type,'toplevel')
    P=[];
    if isatype(p,type)
        P = p;
    end
    while ~isempty(p)
        p = get(handle_j(p),'parent');
        if isatype(p,type)
            P = p;
        end
    end
    p=P;
end

%-------------------------------------------------------------%
function istype = isatype(h,type)

istype = false;
if ischar(type) % ancestor(h,'type'..)
    if strcmpi(get(h,'type'),type) || isa(handle_j(h),type)
        istype=true;
    end
else % ancestor(h,{'type1','type2',..}...)
    % % make sure it's a cell array
    % type = cellstr(type);
    if any(strcmpi(get(h,'type'),type))
        istype = true;
    else
        % check each cell
        for k=1:length(type)
            if isa(handle_j(h),type{k})
                istype=true;
            end
        end
    end
end

%-------------------------------------------------------------%
function type = handle_j(h)
% Attempt to short-circuit the fact that the built in function HANDLE
% is not compilable and is not documented. Apparently it returns the TYPE
% of the handles H as an object class. Here I just return the handle in input
% because the ANCESTOR function seams to be called within Mirone only with H
% as a figure handle. In this case, things seams to work
% and I hope to not f... the ancestor functionality
type = h;

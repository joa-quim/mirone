function uistate = uisuspend_j(fig, setdefaults)
%UISUSPEND_J suspends interactive properties of the figure.
%
%   UISTATE = UISUSPEND_J(FIG) suspends some interactive properties of a 
%   figure window and returns the previous state in the structure
%   UISTATE.  This structure contains information about the figure's
%   WindowButton* functions and the pointer.
%
%   UISTATE = UISUSPEND_J(FIG,TRUE) returns the structure as above but it also
%	contains the ButtonDownFcn's for all children of the figure.
%
%   UISTATE = UISUSPEND_J(FIG,FALSE) returns the structure as above but leaves
%   the current settings unchanged.

%	This is a minimalist replacement of UISUSPEND which only saves properties that
%	matter for Mirone usage.

%   Coffeeright J. Luis 2004-2012

% $Id: uisuspend_j.m 3624 2012-07-22 20:17:57Z j $

	if (nargin < 2)
		setdefaults = true;
		onlyFig = true;
	else
		onlyFig = false;
	end

	if (onlyFig)
		uistate = struct(...
				'figureHandle',          fig, ...
				'WindowButtonMotionFcn', Lwrap(get(fig, 'WindowButtonMotionFcn')), ...
				'WindowButtonDownFcn',   Lwrap(get(fig, 'WindowButtonDownFcn')), ...
				'WindowButtonUpFcn',     Lwrap(get(fig, 'WindowButtonUpFcn')), ...
				'KeyPressFcn',           Lwrap(get(fig, 'KeyPressFcn')), ...
				'Pointer',               get(fig, 'Pointer'), ...
				'PointerShapeCData',     get(fig, 'PointerShapeCData'), ...
				'PointerShapeHotSpot',   get(fig, 'PointerShapeHotSpot') );

		hIm = findobj(fig,'Type','image');  % For now we want only image obj to interrupt the pixval_stsbar ButtonDownFcn

		set(fig, 'WindowButtonMotionFcn', '')
		set(fig, 'WindowButtonDownFcn',   '')
		set(fig, 'WindowButtonUpFcn',     '')
		set(fig, 'KeyPressFcn',           '')
		set(fig, 'Pointer',               'arrow')
		set(fig, 'PointerShapeCData',     get(0, 'DefaultFigurePointerShapeCData'))
		set(fig, 'PointerShapeHotSpot',   [1 1])
		set(hIm, 'ButtonDownFcn',         '')

	else

		% Try to minimize as possible the highly inefficient Octave's findobj
		h1 = findobj(fig,'-depth',1,'Type','axes');
		h2 = findobj(h1,'-depth',1,'Type','image');
		h3 = findobj(h1,'-depth',1,'Type','line');
		h4 = findobj(h1,'-depth',1,'Type','patch');
		h5 = findobj(h1,'-depth',1,'Type','text');
		chi = [h1; h2; h3; h4; h5];

		uistate = struct(...
				'figureHandle',          fig, ...
				'WindowButtonMotionFcn', Lwrap(get(fig, 'WindowButtonMotionFcn')), ...
				'WindowButtonDownFcn',   Lwrap(get(fig, 'WindowButtonDownFcn')), ...
				'WindowButtonUpFcn',     Lwrap(get(fig, 'WindowButtonUpFcn')), ...
				'KeyPressFcn',           Lwrap(get(fig, 'KeyPressFcn')), ...
				'Pointer',               get(fig, 'Pointer'), ...
				'PointerShapeCData',     get(fig, 'PointerShapeCData'), ...
				'PointerShapeHotSpot',   get(fig, 'PointerShapeHotSpot'), ...
				'Children',              chi, ...
				'ButtonDownFcns',        Lwrap(get(chi, {'ButtonDownFcn'})), ...
				'Interruptible',         Lwrap(get(chi, {'Interruptible'})), ...
				'BusyAction',            Lwrap(get(chi, {'BusyAction'})) );

		if (setdefaults)
			set(fig, 'WindowButtonMotionFcn', '')
			set(fig, 'WindowButtonDownFcn',   '')
			set(fig, 'WindowButtonUpFcn',     '')
			set(fig, 'KeyPressFcn',           '')
			set(fig, 'Pointer',               'arrow')
			set(fig, 'PointerShapeCData',     get(0, 'DefaultFigurePointerShapeCData'))
			set(fig, 'PointerShapeHotSpot',   [1 1])
			set(chi, 'ButtonDownFcn',         '')
			set(chi, 'Interruptible',         'on');
			set(chi, 'BusyAction',            'Queue')
		end

	end

% wrap cell arrays in another cell array for passing to the struct command
function x = Lwrap(x)
	if (iscell(x)),  x = {x};    end

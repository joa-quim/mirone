function uistate = uisuspend_j(fig, setdefaults)
%UISUSPENDJ suspends all interactive properties of the figure.
%
%   UISTATE=UISUSPEND_J(FIG) suspends the interactive properties of a 
%   figure window and returns the previous state in the structure
%   UISTATE.  This structure contains information about the figure's
%   WindowButton* functions and the pointer.
%
%   UISTATE=UISUSPEND_J(FIG,TRUE) returns the structure as above but it also
%	contains the ButtonDownFcn's for all children of the figure.
%
%   UISTATE=UISUSPEND_J(FIG,FALSE) returns the structure as above but leaves
%   the current settings unchanged.
  
%   Coffeeright J. Luis 2004-2012


	onlyFig = true;
	if (nargin < 2)
		setdefaults = true;
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
			
			set(fig, 'WindowButtonMotionFcn', get(0, 'DefaultFigureWindowButtonMotionFcn'))
			set(fig, 'WindowButtonDownFcn',   get(0, 'DefaultFigureWindowButtonDownFcn'))
			set(fig, 'WindowButtonUpFcn',     get(0, 'DefaultFigureWindowButtonUpFcn'))
			set(fig, 'KeyPressFcn',           get(0, 'DefaultFigureKeyPressFcn'))
			set(fig, 'Pointer',               get(0, 'DefaultFigurePointer'))
			set(fig, 'PointerShapeCData',     get(0, 'DefaultFigurePointerShapeCData'))
			set(fig, 'PointerShapeHotSpot',   get(0, 'DefaultFigurePointerShapeHotSpot'))
			set(hIm, 'ButtonDownFcn',         '')

	else

		% Try to minimize as possible the highly inefficient Octave's findobj
		h1 = findobj(fig,'-depth',1,'Type','axes');
		h2 = findobj(h1,'-depth',1,'Type','image');
		h3 = findobj(h2,'-depth',1,'Type','line');
		h4 = findobj(h2,'-depth',1,'Type','patch');
		h5 = findobj(h2,'-depth',1,'Type','text');
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
				'BusyAction',            Lwrap(get(chi, {'BusyAction'})), ...
				'UIContextMenu',         Lwrap(get(chi, {'UIContextMenu'})) );

		if (setdefaults)
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

	end

% wrap cell arrays in another cell array for passing to the struct command
function x = Lwrap(x)
	if (iscell(x)),  x = {x};    end

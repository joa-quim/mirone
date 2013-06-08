function [out1,out2,out3] = click_e_point(arg1,arg2,arg3)
%CLICK_E_PUT Graphical input from mouse.
%
%   [X,Y] = CLICK_E_PUT(N) gets N points from the current axes and returns 
%   the X- and Y-coordinates in vectors X and Y. 
%	Hitting the return key always terminates the input.
%
%   [X,Y] = CLICK_E_PUT gathers an unlimited number of points until the
%   return key is pressed.
% 
%   [X,Y,BUTTON] = CLICK_E_PUT(N) returns a third result, BUTTON, that 
%   contains a vector of integers specifying which mouse button was
%   used (1,2,3 from left) or ASCII numbers if a key on the keyboard
%   was used.
%
%   [...] = CLICK_E_PUT(N,POINTER)     or  [...] = CLICK_E_PUT(POINTER)
%   where POINTER is either a valid pointer name or a 16x16 matrix deffining a custom pointer
%
%   [...] = CLICK_E_PUT(N,POINTER,PSHS)
%   where PSHS is a two element vector, sets the 'PointerShapeHotSpot' property

%   Coffeeright J. Luis 2012

	PSHS = [];
	if (nargin >= 2)
		if ~(ischar(arg2) || size(arg2,1) == 16)
			error('click_e_input error: second argument must be a string containing the pointer name');
		elseif size(arg2,1) == 16
			new_pointer = 'custom';
			nargs = 1;
			if (nargin == 3 && numel(arg3) == 2)    % Got a 'PointerShapeHotSpot'
				PSHS = arg3;
			end
		else
			new_pointer = arg2;
			nargs = 1;				% So that the rest of the code may work as if this change has not been made
		end
	elseif (nargin == 1)
		if ischar(arg1)
			new_pointer = arg1;
			nargs = 0;				% So that the rest of the code may work as if this change has not been made
		elseif length(arg1) == 16
			new_pointer = 'custom';
			nargs = 0;
		else
			new_pointer = 'fullcrosshair';
			nargs = 1;
		end
	else
		new_pointer = 'fullcrosshair';
		nargs = 0;
	end

	out1 = []; out2 = []; out3 = []; y = [];

	hFig = gcf;
	figure(hFig);
	b = [];
	if nargs == 0
		nWishPts = -1;
	else
		nWishPts = arg1;
		if  (ischar(nWishPts) || numel(nWishPts) ~= 1 || ...
				~(fix(nWishPts) == nWishPts) || nWishPts <= 0)
			error('click_e_input error: Requires a positive integer.')
		end
	end

	% Remove figure button functions
	state = uisuspend_j(hFig);
	if strcmp(new_pointer,'custom')
		set(hFig,'pointer',new_pointer,'PointerShapeCData',arg2)
		if (~isempty(PSHS)),    set(hFig,'PointerShapeHotSpot',PSHS);    end
	else
		set(hFig,'pointer',new_pointer);
	end
	fig_units = get(hFig,'units');
	current_char = 0;

	while (nWishPts ~= 0)
		try
			keydown = waitforbuttonpress;
		catch
			if (~ishandle(hFig))
				error('click_e_input: figure deleted');
			end
			set(hFig,'units',fig_units);
			uirestore_j(state, 'nochildren');
			error(['Interrupted: ' lasterr]);
		end
		drawnow

		ptr_fig = get(0,'CurrentFigure');
		if (ptr_fig == hFig)
			if keydown
				current_char = get(hFig, 'CurrentCharacter');
				button = abs(current_char);
				scrn_pt = get(0, 'PointerLocation');
				set(hFig,'units','pixels')
				figPos = get(hFig, 'Position');
				pt = [scrn_pt(1) - figPos(1), scrn_pt(2) - figPos(2)];
				set(hFig,'CurrentPoint',pt);
			else
				button = get(hFig, 'SelectionType');
				if strcmp(button,'open')
					if (isempty(b))		% A double click when the b vector is still empty would cause an error
						button = 1;		% Pretend it was a normal click
					else
						button = b(numel(b));
					end
				elseif strcmp(button,'normal')
					button = 1;
				elseif strcmp(button,'extend')
					button = 2;
				elseif strcmp(button,'alt')
					button = 3;
				else
					error('Invalid mouse selection.')
				end
			end
			pt = get(gca, 'CurrentPoint');

			nWishPts = nWishPts - 1;

			if (current_char == 13)		% Return key
				break
			end

			out1 = [out1; pt(1,1)];
			y = [y; pt(1,2)];
			b = [b; button];
		end
	end

	uirestore_j(state, 'nochildren');
	set(hFig,'units',fig_units);

	if (nargout > 1)
		out2 = y;
		if (nargout > 2),   out3 = b;   end
	else
		out1 = [out1 y];
	end

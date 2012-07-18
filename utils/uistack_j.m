function uistack_j(hand, opt, step)
%UISTACK_J Restack objects children of AXES.
%   UISTACK_J(H) raises the stacking order of the objects, H. 
%   UISTACK_J(H,OPT) where OPT is 'up','down','top','bottom'
%   UISTACK_J(H,OPT,STEP) where STEP is the distance to move 'up' and 'down'

%   UISTACK replacement that should only be used with handles children of
%   Image and Surface objects (e.g. NO UICONTROLS).
%   It also fixes a bug on the ML R13 function that doesn't allow
%   restacking between different object types.

% $Id: $

	if (~all(ishandle(hand)))
		error('UISTACK_JUISTACK_J: invalid handle in input.');
	end
	if (nargin == 1)
		opt = 'up';		step = 1;
	elseif (nargin == 2)
		step = 1;
	end

	Parent = get(hand,{'Parent'});
	Parent = [Parent{:}];

	children = allchild(Parent);
	handPos = find(children == hand);

	if (opt(1) == 'u')			% 'up'
		restacked = [-ones(step,1); children];
		handPos = handPos + step;
		for (n = 1:numel(hand))
			ind = [(1:handPos(n)-step-1) handPos(n) (handPos(n)-step:handPos(n)-1) (handPos(n)+1:numel(restacked))];
			restacked = restacked(ind);
		end
		restacked(restacked == -1) = [];

	elseif (opt(1) == 'd')		% 'down'
		restacked = [children; -ones(step,1)];
		for (n = numel(hand):-1:1)
			ind = [(1:handPos(n)-1) (handPos(n)+1:handPos(n)+step) handPos(n) (handPos(n)+step+1:numel(restacked))];
			restacked = restacked(ind);
		end
		restacked(restacked == -1) = [];

	elseif (opt(1) == 'b')		% 'bottom'
		children(handPos) = [];
		restacked = [children; hand];

	elseif (opt(1) == 't')		% 'top'
		children(handPos) = [];
		restacked = [hand; children];

	else
		error('UISTACK_J: invalid stack option');      
	end

	% If image, light, or surface object exists, put them at the bottom of the stack
	hImg = findobj(Parent,'type','image');
	hLight = findobj(Parent,'type','light');
	hSurf = findobj(Parent,'type','surface');
	hAx_bot = [hImg; hLight; hSurf];
	if (~isempty(hAx_bot))
		restacked(restacked == hAx_bot) = [];
		restacked = [restacked; hAx_bot];
	end

	set(Parent,'children',restacked);

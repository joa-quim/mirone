function [xt, yt] = make_arrow(lh , hscale, vscale, ah, vFac, full)
% returns x and y coord vectors for an arrow or just of the arrow head
%
% LH is either a line handle or a [x1 x2; y1 y2] matrix with the arrow tail coordinates
% HSCALE & VSCALE are the core of this business. If only HSCALE is provided, VSCALE = HSCALE
% The following snipet illustrates how to compute them
%  	axLims = [get(axesHandle,'XLim') get(axesHandle,'YLim')];
% 	oldUnits = get(axesHandle,'Units');			set(axesHandle,'Units','points');
% 	Pos = get(axesHandle,'Position');			set(axesHandle,'Units',oldUnits);
% 	vscale = 1/Pos(4) * diff(axLims(1:2));	hscale = 1/Pos(3) * diff(axLims(3:4));
% 	vscale = (vscale + hscale) / 2;			hscale = vscale;	% For not having a head direction dependency
% 
% AH is the arrow head height in points
% VFAC represents the arrow's head V shape factor. VFAC = 1 makes a triangular head. Default is 1.3
% FULL = 'yes' (default) returns the coord of the full arrow. Otherwise only
% the coords of the arrow head are returned
%
% Reworked version of the ML buildArrow

% Joaquim Luis	27-Aug-2008

% $Id$

	n_args = nargin;
	if (n_args == 2)
		vscale = hscale;		ah = 12;	vFac = 1.3;		full = 'yes';
	elseif (n_args == 3)
		ah = 12;	vFac = 1.3;		full = 'yes';
	elseif (n_args == 4)
		vFac = 1.3;		full = 'yes';
	elseif (n_args == 5)
		full = 'yes';
	elseif (n_args == 0 || n_args > 6)
		error('MAKE_ARROW: Wrong number of arguments')
	end

	if (ishandle(lh))		% get the line width
		try
			lw = get(lh, 'LineWidth');		lw_2 = lw / 2;
		catch
			error('MAKE_ARROW: First argument was a valid line handle')
		end
		X = get(lh, 'XData');	Y = get(lh, 'YData');
	else
		lw = 1;					lw_2 = lw / 2;		% pretend the line width is 1 pt
		X = lh(1,:);			Y = lh(2,:);
	end

	% calculate the x, and y lengths (in points/map units)
	dx = (X(end) - X(end - 1)) / hscale;
	dy = (Y(end) - Y(end - 1)) / vscale;

	% calculate the cosine and sine
	hy = (dx^2 + dy^2)^.5;
	co =  dx / hy;		si =  dy / hy;

	aw = ah * .66 + lw_2;	% 3 : 2 aspect ratio
	ahV = ah * vFac;		% Sets V shape factor on arrow head

	% determine arrowhead style
	xt = [ -ah,  -ahV,  0,  -ahV, -ah ];
	yt = [ lw_2,  aw/2, 0, -aw/2, -lw_2 ];

	% Offset of the size of arrow head to avoid increasing the total length of that extra size
	xt_off = [ ah ah ];
	yt_off = [ 0  0 ];

	% rotate the head based on the slope of the last line
	foo = [co -si; si  co] * [xt; yt];

	% Rotate the offset
	foo_off = [co -si; si  co] * [xt_off; yt_off];

	% convert points back to to data units and add in the offset
	xt = foo(1,:) * hscale + X(end);
	yt = foo(2,:) * vscale + Y(end);

	if (full(1) == 'y')		% Join the line coords as the tail of the full arrow
		foot_x = ([foo(1,end) foo(1,1)] + foo_off(1,:)) * hscale + X(end-1);
		foot_y = ([foo(2,end) foo(2,1)] + foo_off(2,:)) * vscale + Y(end-1);
		xt = [xt foot_x xt(1)];		yt = [yt foot_y yt(1)];

		if ( hy < ah )		% When head is still longer than dist from initial to final pt draw the header only
			xt(end-1) = xt(end);		xt(end-2) = xt(end-3);
			yt(end-1) = yt(end);		yt(end-2) = yt(end-3);
		end
	end

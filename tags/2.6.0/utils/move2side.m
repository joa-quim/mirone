function move2side(hFigStatic, hFigMov, opt)
% Put a second figure right to the side of a reference (standing) figure
%
% Try to position figure "hFigMov" glued to the right, left or bottom of "hFigStatic" figure
% OPT is an optional char string with 'left', 'right' or 'bottom' (if not provided, defaults to 'right')
% that selects at which side to position the moving figure.
%
% move2side(hFigMov, opt) works nearly as above but without a reference StaticFigure
%	OPT can now have any of the MOVEGUI position strings and this function works like it.
%
% Some cases intercept the movegui ML command, but I had to because TMW simply doesn't
% seam to Grok the concept of stability between versions and OSs.
% movegui incredibly f. keeps behaving slightly different across versions.
%
%	EXAMPLE:
%		load clown
%		imagesc(X);		hS = gcf;	set(hS,'colormap',map)
%		hL = warndlg('I am a clown');
%		hR = warndlg('I am a clown');
%		move2side(hS, hL, 'left')
%		move2side(hS, hR)
%
%   AUTHOR
%		Joaquim Luis (jluis@ualg.pt)    13-Sep-2007
%		University of Algarve

%	Copyright (c) 2004-2012 by J. Luis
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

	n_args = nargin;
	if (n_args < 2)
		error('move2side: need at least two inputs')
	elseif (n_args == 2)
		if (isa(hFigMov,'char'))
			opt = hFigMov;
			hFigMov = hFigStatic;
			hFigStatic = false;
		else
			opt = 'right';
		end
	elseif (n_args == 3 && ~ischar(opt))
		error('move2side: third argument must be a character string')
	end

	if ( ~(ishandle(hFigMov) && strcmp(get(hFigMov,'Type'), 'figure')) )
		error('move2side: first argument must be a valid figure handle')
	end
	if (hFigStatic ~= 0 && ~(ishandle(hFigStatic) && strcmp(get(hFigStatic,'Type'), 'figure')) )
		error('move2side: second argument must be a valid figure handle')
	end

	% get the screen size
	units_bak = get(0, 'units');		set(0, 'units', 'pixels');
	ecran = get(0,'ScreenSize');		set(0, 'units', units_bak);
	
	% save original figures units and temp set them to pixels
	if (hFigStatic ~= 0)
		FigStaticUnit = get(hFigStatic,'Units');
		set(hFigStatic,'Units','Pixels')
	end
	FigMovUnit = get(hFigMov,'Units');	set(hFigMov,'Units','Pixels')
	
	% Get figures dimensions
	posFigMov = get(hFigMov,'Pos');
	outPosFigMov = get(hFigMov,'outerposition');
	if (hFigStatic ~= 0)
		posFigStatic = get(hFigStatic,'Pos');
	else
		if (lower(opt(1)) == 'r')
			posFigStatic = [ecran(3) ecran(2:4)];
		else
			posFigStatic = [-1 0 ecran(3:4)];
		end
	end
	
	refine_yLL = true;				% For the horizontal cases
	if (lower(opt(1)) == 'r')		% Put moving figure on the RIGHT side of reference figure
		xLL = posFigStatic(1) + posFigStatic(3) + 6;	% + 6 is empiric
		xLR = xLL + posFigMov(3);
		if (xLR > ecran(3))         % If figure is partially out, bring it totally into screen
			xLL = ecran(3) - posFigMov(3);
		end
	elseif (lower(opt(1)) == 'l')	% Put moving figure on the LEFT side of reference figure
		xLL = posFigStatic(1) - posFigMov(3) - 6;		% + 6 is empiric
		if (xLL < 0)				% If figure is partially out, bring it totally into screen
			xLL = 4;				% 4 is nicier than 0
		end
	elseif (lower(opt(1)) == 'b')	% Put moving figure on the BOTOM side of reference figure
		[posFigMov, bars_height] = move_to(hFigMov, 'south', ecran);
		xLL = posFigStatic(1) + posFigStatic(3)/2  - posFigMov(3)/2;	% More or less centered
		yLL = posFigStatic(2) - posFigMov(4) - bars_height - 3;		% But this can leak through the screen bottom ...
		yLL = max(yLL, posFigMov(2)+35);		% Take the highest of the two estimations. 35 is for bottom bar
		refine_yLL = false;
	else							% Use any of the movegui postion options (no StaticFig)
		posFigMov = move_to(hFigMov, lower(opt), ecran);
		xLL = posFigMov(1);
		yLL = posFigMov(2);
		refine_yLL = false;
	end

	if (refine_yLL),	yLL = (posFigStatic(2) + posFigStatic(4)/2+12) - (posFigMov(4) / 2 - 22);	end
	if ( (yLL + outPosFigMov(4) + 5) > ecran(4) )		% Figure is partially out from top (5 is au-pif)
		yLL = ecran(4) - outPosFigMov(4) + 4;
	end
	set(hFigMov,'Pos',[xLL yLL posFigMov(3:4)])
	
	% Reset original figures units
	if (hFigStatic ~= 0),	set(hFigStatic,'Units',FigStaticUnit),		end
	set(hFigMov,'Units',FigMovUnit)


% ----------------------------------------------------------------------
function [new_pos, bars_height] = move_to(fig, position, ecran)
% Bit of code extracted and adapted from ML's movegui.
	oldpos  = get(fig, 'position');
	wfudge =  6;
	hfudge = 24;

	if ~isempty(findall(fig,'-depth',1,'type','uimenu'))
		hfudge = hfudge + 32;
	end

	n_toolbars = numel(findall(fig,'-depth',1,'type','uitoolbar'));
	if (n_toolbars > 0)
		hfudge = hfudge + 24 * n_toolbars;
	end

	oldpos(3) = oldpos(3) + wfudge;
	oldpos(4) = oldpos(4) + hfudge;

	fleft   = oldpos(1);	fwidth  = oldpos(3);
	fbottom = oldpos(2);	fheight = oldpos(4);

	% make sure the figure is not bigger than the screen size
	fwidth = min(fwidth, ecran(3));		fheight = min(fheight, ecran(4));

	rwidth  = ecran(3) - fwidth;		% remaining width
	rheight = ecran(4) - fheight;		% remaining height

	switch position
		case 'north',		new_pos = [rwidth/2,   rheight];
		case 'south',		new_pos = [rwidth/2,         0];
		case 'east',		new_pos = [  rwidth, rheight/2];
		case 'west',		new_pos = [       0, rheight/2];
		case 'northeast',	new_pos = [  rwidth,   rheight];
		case 'southeast',	new_pos = [  rwidth,         0];
		case 'northwest',	new_pos = [       0,   rheight];
		case 'southwest',	new_pos = [       0,         0];
		case 'center',		new_pos = [rwidth/2, rheight/2];
		case 'onscreen'
			if (fleft < 0),			fleft = 0;		end
			if (fbottom < 0),		fbottom = 0;	end
			if (fleft > rwidth),	fleft = rwidth;	end
			if (fbottom > rheight),	fbottom = rheight;	end
			new_pos = [fleft, fbottom];
	end

	new_pos(3:4) = [fwidth - wfudge, fheight - hfudge];
	if (nargout == 2),		bars_height = hfudge;	end

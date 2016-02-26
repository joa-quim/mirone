function cayula_cornillon(varargin)
% ...

%	Copyright (c) 2004-2016 by J. Luis
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

% $Id: $

	if (isempty(varargin)),		return,		end
	hMirHand = varargin{1};
	hPolyg = varargin{2};

	minPopProp = 0.20;				% 
	minPopMeanDiff = 0.4;
	minTheta = 0.70;
	minSinglePopCohesion = 0.90;
	minGlobalPopCohesion = 0.92;

	[X,Y,Z,head] = load_grd(hMirHand);
	if isempty(Z),		return,		end		% An error message was already issued

	x = get(hPolyg,'XData');	y = get(hPolyg,'YData');
	xp(1) = min(x);		xp(2) = max(x);
	yp(1) = min(y);		yp(2) = max(y);
	if (diff(xp) < 1)		% TEMPORARIO. Para fazer so uma janela
		rect_crop = [xp(1) yp(1) (xp(2) - xp(1)) (yp(2) - yp(1))];		% Insure a known order
		[w, r_c] = cropimg(head(1:2),head(3:4),Z,rect_crop,'out_grid');
		head(2) = head(1) + (r_c(4)-1)*head(8);			head(1) = head(1) + (r_c(3)-1)*head(8);
		head(4) = head(3) + (r_c(2)-1)*head(9);			head(3) = head(3) + (r_c(1)-1)*head(9);

		[x, y, z, exitType] = getFrontInWindow(w, head, minTheta, minPopProp, minPopMeanDiff, ...
											 minSinglePopCohesion, minGlobalPopCohesion);
		if (~exitType)
			line('XData',x, 'YData',y, 'Parent',hMirHand.axes1, 'LineWidth',hMirHand.DefLineThick, ...
				'Color',hMirHand.DefLineColor, 'userdata',z,'Tag','SSTfront');
		end
		return
	end

	% Here we follow the algorithm described in Nieto et al 2012
	[n_rows, n_cols] = size(Z);
	winW16 = 16;	winW32 = 2*16;		winW48 = 3*16;
	s = 0;				% s = 1 means subwindows do NOT share a common border. With s = 0 they do.
	xSide16 = winW16 * head(8);
	ySide16 = winW16 * head(9);
	xSide32 = (winW32 - s) * head(8);			% X and Y sizes of the 32x32 window in real coordinates
	ySide32 = (winW32 - s) * head(9);			% -1 because coords are in grid registration
	nWinRows = floor(n_rows / winW16);		nWinCols = floor(n_cols / winW16);
	for (wRow = 1:nWinRows-2)
		r1 = (wRow - 1) * winW16 + 1;		r2 = r1 + winW48 - s;	% Start and stop indices and coords of current window
		y0 = head(3) + (wRow-1) * ySide16;				%y1 = y0 + (winW48-1) * head(9);
		for (wCol = 1:nWinCols-2)
			c1 = (wCol - 1) * winW16 + 1;	c2 = c1 + winW48 - s;
			x0 = head(1) + (wCol-1) * xSide16;			%x1 = x0 + (winW48-1) * head(8);
			wPad = Z(r1:r2, c1:c2);						% 48x49 (or 48x48 if s == 1) Window
			% |4|3|
			% |1|2|
			% Rows and columns of the 4 sub-windows (avoids a double loop). We loop along the order of the
			% above sketch. There will be a common region when we slide the 33x33 window with shifts of 16
			% When we at 1 (a 33x33 window) that common zone is depicted at position 1 form the sketch below,
			% and so forth for the other three positions.
			rr = [1 1 2 2];		cc = [1 2 2 1];

			% Set the index of the four corner windows of 17x17 representing each the same physical region
			% when we loop 4 times the sliding window of 33x33 in the larger 49x49 window (above sketch).
			% |2|1|
			% |3|4|
			if (s == 1),	corners = [17 32 17 32; 17 32 1 16; 1 16 1 16; 1 16 17 32];		% Less good
			else			corners = [17 33 17 33; 17 33 1 17; 1 17 1 17; 1 17 17 33];
			end

			for (k = 1:4)					% Loop over the 4 sliding 32x32 sub-windows of the larger 48x48 one.
				m1 = (rr(k) - 1) * winW16 + 1;			m2 = m1 + 2*winW16 - s;		% indices of the sliding 33x33 window
				n1 = (cc(k) - 1) * winW16 + 1;			n2 = n1 + 2*winW16 - s;
				w = wPad(m1:m2, n1:n2);					% Sub-window with size 33x33
				subWinX0 = x0 + (cc(k)-1) * xSide16;	subWinX1 = subWinX0 + xSide32;	% Corner coordinates
				subWinY0 = y0 + (rr(k)-1) * ySide16;	subWinY1 = subWinY0 + ySide32;
				R = [subWinX0 subWinX1 subWinY0 subWinY1];
				[x, y, z, exitType] = getFrontInWindow(w, R, minTheta, minPopProp, minPopMeanDiff, ...
													   minSinglePopCohesion, minGlobalPopCohesion, corners(k,:));
				if (~exitType)
					line('XData',x, 'YData',y, 'Parent',hMirHand.axes1, 'LineWidth',hMirHand.DefLineThick, ...
						'Color',hMirHand.DefLineColor, 'userdata',z, 'Tag','SSTfront');
					%rectX = [subWinX0 subWinX0 subWinX1 subWinX1 subWinX0];
					%rectY = [subWinY0 subWinY1 subWinY1 subWinY0 subWinY0];
					%hl = line('XData',rectX, 'YData',rectY, 'Parent',hMirHand.axes1);
					%hl = line('XData',[x0 x0 x1 x1 x0], 'YData',[y0 y1 y1 y0 y0], 'Parent',hMirHand.axes1);
					%draw_funs(hl,'line_uicontext')		% Set lines's uicontextmenu
				end
			end

		end
		drawnow
	end

% ---------------------------------------------------------------------------------------------------
function [xdata, ydata, z, exitType] = getFrontInWindow(w, head, minTheta, minPopProp, minPopMeanDiff, ...
                                       minSinglePopCohesion, minGlobalPopCohesion, corners)
% ...
% If CORNERS exist than it holds the start and stop index of 1/4 sub-window of W inside which
% we will retain the fronts found. Detections outside this sub-window are ignored
%
% Here we follow, patially, what is donne in MGET (FrontUtils.cpp by Jason Roberts) that in turn,
% according to program's notes, followed an original Fortran code by Dave Ulman.

	if (nargin < 8),	corners = [];	end
	xdata = [];		ydata = [];		z = [];		exitType = 0;

	mask = isnan(w);	
	haveNaNs = any(mask(:));
	n_NaNs = 0;
	if (haveNaNs)
		n_NaNs = sum(mask(:));
		if (n_NaNs / numel(w) > 0.5)			% Window must be at least > half filled
			exitType = -1;
			return
		end
	end

	mi_ma = double([min(w(:)) max(w(:))]);			% <=== MUDAR PARA GRDUTILS
	n = ceil(diff(mi_ma) / 0.02);
	[y,xout] = histo_m('hist',w(:), n, mi_ma);

	thresValue = xout(1);	totalCount = numel(w) - n_NaNs;
	threshPopACount = 0;	threshSeparation = -1;
	threshPopAMean = 0;		threshPopBMean = 0;

	w(mask) = 0;				% Replace NaNs with 0's
	totalSum = sum(w(:));
	totalSumSquares = sum(w(:) .* w(:));

	for (k = 2:n-1)				% Ignore first and last as candidates
		popASum = sum(y(1:k) .* xout(1:k));
		popBSum = sum(y(k+1:end) .* xout(k+1:end));
		popACount = sum(y(1:k));
		popBCount = sum(y(k+1:end));

		popAMean = popASum / popACount;
		popBMean = popBSum / popBCount;
		separation = popACount * popBCount * (popAMean - popBMean) * (popAMean - popBMean);
		if (separation > threshSeparation)
			threshSeparation = separation;
			thresValue = xout(k);
			threshPopACount = popACount;
			threshPopAMean = popAMean;
			threshPopBMean = popBMean;
		end
	end

	% Continue only if the proportional size of the smaller population exceeds the minimum allowed value.
	% This test corresponds to equation 14 in Cayula-Cornillon
	if (threshPopACount / totalCount < minPopProp)
		exitType = 1;
		return
	end
	if (1.0 - threshPopACount / totalCount < minPopProp)
		exitType = 1;
		return
	end
	% Abort this window if the difference in the populations' means is less than a minimum value.
	if (threshPopBMean - threshPopAMean < minPopMeanDiff)
		exitType = 2;
		return
	end

	% Calculate the criterion function THETA(TAUopt) discussed on page 72 of the paper.
	totalMean = totalSum / totalCount;
	variance = totalSumSquares - (totalMean * totalMean * totalCount);
	theta = threshSeparation / (variance * totalCount);
	if (theta < minTheta)
		exitType = 3;
		return
	end

	% Cohesion
	% Count the number of times a population A cell is immediately adjacent to another population A cell,
	% and the same for population B. A cell can be adjacent on four sides. Count only two of them
	% (bottom and right side) because doing all four would be redundant. Do not count diagonal neighbors.
	countANextToA = 0;		countBNextToB = 0;
	countANextToAOrB = 0;	countBNextToAOrB = 0;
	[n_rows, n_cols] = size(w);
	for (col = 1:n_cols-1)
		for (row = 1:n_rows-1)
			if (haveNaNs && (mask(row,col) || mask(row+1,col) || mask(row, col + 1)))
				continue
			end
			% Examine the bottom neighbor
			if (w(row,col) <= thresValue)
				countANextToAOrB = countANextToAOrB + 1;
				if (w(row + 1, col) <= thresValue)
					countANextToA = countANextToA + 1;
				end
			else
				countBNextToAOrB = countBNextToAOrB + 1;
				if (w(row + 1, col) > thresValue)
					countBNextToB = countBNextToB + 1;
				end
			end
			% Examine the right neighbor
			if (w(row,col) <= thresValue)
				countANextToAOrB = countANextToAOrB + 1;
				if (w(row, col + 1) <= thresValue)
					countANextToA = countANextToA + 1;
				end
			else
				countBNextToAOrB = countBNextToAOrB + 1;
				if (w(row, col + 1) > thresValue)
					countBNextToB = countBNextToB + 1;
				end
			end

		end
	end
	
	popACohesion = countANextToA / countANextToAOrB;
	popBCohesion = countBNextToB / countBNextToAOrB;
	globalCohesion = (countANextToA + countBNextToB) / (countANextToAOrB + countBNextToAOrB);
	if (popACohesion < minSinglePopCohesion)
		exitType = 4;
		return
	end
	if (popBCohesion < minSinglePopCohesion)
		exitType = 4;
		return
	end
	if (globalCohesion < minGlobalPopCohesion)
		exitType = 4;
		return
	end

	% OK, if we reach here we have a front. Compute its contour. This part is not discussed in any
	% of the papers I saw and is different from the aproach followed in MGET.
	X = linspace(head(1),head(2), n_cols);
	Y = linspace(head(3),head(4), n_rows);
	if (isempty(corners))		% Use the full window
		w = double(w);			% It has to be because of the contourc
		if (haveNaNs)
			w(w == 0) = NaN;	% Need to restore the NaNs to not invent new contours around zeros
		end
		c = contourc(X,Y,w,[thresValue thresValue]);
	else
		% The four corners have these indices [17 32 17 32; 17 32 1 16; 1 16 1 16; 1 16 17 32]
		% and the variable corners has one of its rows (the current to be retained sub-window)
		X = X(corners(3):corners(4));
		Y = Y(corners(1):corners(2));
		w = double(w(corners(1):corners(2), corners(3):corners(4)));
		if (haveNaNs)
			w(w == 0) = NaN;
		end
		c = contourc(X,Y,w,[thresValue thresValue]);
	end
	limit = size(c,2);
	i = 1;
	while(i < limit)
		npoints = c(2,i);
		x = c(1,i+1:i+npoints);		y = c(2,i+1:i+npoints);
		if (x(1) == x(end) || npoints < 7)		% We don't want closed lines nor very short ones
			i = i+npoints+1;
			continue
		end
		if (isempty(xdata))			% First contour
			xdata = x;	ydata = y;
		else
			xdata = [xdata NaN x];
			ydata = [ydata NaN y];
		end
		i = i+npoints+1;
		%xdata = c(1,2:npoints+1);		ydata = c(2,2:npoints+1);
	end
	z = thresValue;
	
	if (isempty(xdata))
		exitType = 5;
	end


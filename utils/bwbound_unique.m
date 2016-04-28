function B = bwbound_unique(B)
% Remove duplicate points from the output of BWBOUNDARIES
%
% When bwboundaries vectorize raster lines it returns vectors that go arround the raster line.
% This results in duplicated points which screw up further computations relying on polyline uiqueness.
% To complicate things further, the vector lines do not always start at one of the line ends.
% This function tries to remove the duplicates and re-order vertex so that we have a decent line.

%	Copyright (c) 2004-2013 by J. Luis
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

% $Id: bwbound_unique.m 4348 2014-03-14 17:34:54Z j $

	c = false(numel(B),1);

	for (k = 1:numel(B))
		y = B{k}(:,1);		x = B{k}(:,2);
		n_pts = numel(x);
		if (n_pts <= 2)
			if ((n_pts == 2 && x(1) == x(2) && y(1) == y(2)) || n_pts == 1)		% May happen that it's a duplicate point
				c(k) = true;		% Mark for deletion
			end
			continue
		end
		n = fix(n_pts/2) + 1; 
		if ( (y(n-1) == y(n+1)) && (x(n-1) == x(n+1)) )
			x = x(1:n);  	y = y(1:n);
		else
			dfx = diff(x(1:2:end));
			dfy = diff(y(1:2:end));
			ind = find(dfx == 0 & dfy == 0);
			if (isempty(ind))
				dfx = diff(x(2:2:end));
				dfy = diff(y(2:2:end));
				ind = find(dfx == 0 & dfy == 0);
			end

			ind = 2 * ind + 1;			% The differences start a 2nd pt

			if (isempty(ind)),		continue,		end

			% The above test is sometimes not good enough. Do one more to detect 'excursions'
			if (numel(ind) == 1)
				try
					confirm = (x(ind-3) == x(ind+2) && y(ind-3) == y(ind+2));
					if (~confirm),	continue,		end
				end
			end
			
			if (ind(1) >= n)					% Turning point is at or ahead of mid point
				x = x((ind-n+1):ind);		y = y((ind-n+1):ind);
			else 							% (ind < n)	 Before the mid point
				x = x(ind:ind+n-1);			y = y(ind:ind+n-1);
			end
		end
		B{k} = [y x];
	end

	if (any(c)),	B(c) = [];		end

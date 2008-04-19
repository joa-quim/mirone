function B = bwbound_unique(B)
% Remove duplicate points from the output of BWBOUNDARIES
%
% When bwboundaries vectorize raster lines it returns vectors that go arround the raster line.
% This results in duplicated points which screw up further computations relying on polyline uiqueness.
% To complicate things further, the vector lines do not always start at one of the line ends.
% This function tries to remove the duplicates and re-order vertex so that we have a decent line.

	for (k = 1:numel(B))
		got_it = true;
		y = B{k}(:,1);		x = B{k}(:,2);
		n_pts = numel(x);
		if ( n_pts <= 2 ),	continue,		end
		n = fix(n_pts/2) + 1; 
		if ( (y(n-1) == y(n+1)) && (x(n-1) == x(n+1)) )
			x = x(1:n);  	y = y(1:n);
		else
			dfx = diff(x(1:2:end));
			dfy = diff(y(1:2:end));
			ind = find(dfx == 0 & dfy == 0);
			ind = 2 * ind + 1;
			if (isempty(ind))
				dfx = diff(x(2:2:end));
				dfy = diff(y(2:2:end));
				ind = find(dfx == 0 & dfy == 0);
				ind = 2 * ind + 1;			% The differences start a 2nd pt
			end
			if (~isempty(ind)),		ind = ind(1);	end
			
			if (ind >= n)					% Turning point is at or ahead of mid point
				x = x((ind-n+1):ind);		y = y((ind-n+1):ind);
			else 							% (ind < n)	 Before the mid point
				x = x(ind:ind+n-1);			y = y(ind:ind+n-1);
			end
		end
		if (got_it),	B{k} = [y x];	end
	end

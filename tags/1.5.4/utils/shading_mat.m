function    img = shading_mat(img,R,scale)
% New version that calls a mex function named mex_illuminate
% With this we gained in performance by a factor > 2

	if (nargin == 2),    scale = 'scale';    end

	rows = size(img,1);
	cols = size(img,2);
	if (strcmp(scale,'scale'))      % Scale R into [-1 1] interval
		if (~isa(R,'double')),      R = double(R);  end
		Rmax = max(R(:));   Rmin = min(R(:));
		esc = 1 / (Rmax -Rmin);
		R = (-1 + 2*((R-Rmin)*esc)) * 0.95;
	end

	% If R is not of the same size of 'img' resize it to be. May "legally" occur with drappings
	if ( ~isequal(size(R), [rows cols]) )
		R = cvlib_mex('resize',R,[rows cols]);
	end

	[r,g,b] = mex_illuminate(img,R);        % mex_illuminate should output img
	img = reshape([r g b],rows,cols,3);

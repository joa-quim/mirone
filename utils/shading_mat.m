function    img = shading_mat(img,R,scale)
% New version that calls a mex function named mex_illuminate
% With this we gained in performance by a factor > 2

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

% $Id: shading_mat.m 3499 2012-04-22 23:39:22Z j $

	if (nargin == 2),    scale = 'scale';    end

	rows = size(img,1);
	cols = size(img,2);
	if (strcmp(scale,'scale'))		% Scale R into [-1 1] interval
		if (~isa(R,'double')),		R = double(R);  end
		Rmax = max(R(:));			Rmin = min(R(:));
		esc = 1 / (Rmax -Rmin);
		R = (-1 + 2*((R-Rmin)*esc)) * 0.95;
	end

	% If R is not of the same size of 'img' resize it to be. May "legally" occur with drappings
	if ( ~isequal(size(R), [rows cols]) )
		R = cvlib_mex('resize',R,[rows cols]);
	end

	img = mex_illuminate(img,R);		% It now uses OpenMP

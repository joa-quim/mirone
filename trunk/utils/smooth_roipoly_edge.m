function Z_filt = smooth_roipoly_edge(head, have_local_nans, Z, Z_rect, Z_filt, r_c, mask, thick)
% Apply a ROI mask and smooth the transitional.
%
% When doing a ROI filtering the transition between the filtered and non-filt region
% can result in very sharp transition, makein an uggly effect. This function smooths
% out the zone over the edges of the mask.
% have_local_nans	True if Z_rect has NaNs somewhere. If not known, set it to []
% Z			The full grid array
% Z_rect	The rectangular zone that holds the ROI (not filtered)
% Z_filt	The rectangular zone that holds the ROI (already filtered)
% r_c		The row colum vector issued by cropimg
% mask		The mask as calculated by poipoly_j
% thick		The thickness of the smooth transitional zone in number of grid cells

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

	if (nargin == 7)	thick = 3;		end

	mask_b = mask_border(mask,thick);		mask_b = mask_b(:);

	X = (head(1) + (r_c(3)-1)*head(8)):head(8):(head(1) + (r_c(4)-1)*head(8));
	Y = (head(3) + (r_c(1)-1)*head(9)):head(9):(head(3) + (r_c(2)-1)*head(9));
	opt_R = sprintf('-R%.10f/%.10f/%.10f/%.10f', X(1), X(end), Y(1), Y(end));
	opt_I = sprintf('-I%.10f/%.10f',head(8),head(9));
	Z_filt(~mask) = Z_rect(~mask);		% Outside mask zone Z_filt is set to original value
	
	Z_filt = Z_filt(:);     Z_filt(mask_b) = [];
	[X,Y] = meshgrid(X,Y);
	X = X(:);   X(mask_b) = [];
	Y = Y(:);   Y(mask_b) = [];
	Z_filt = gmtmbgrid_m(X, Y, Z_filt, opt_R, opt_I, '-Mz');    clear X Y;

	if (have_local_nans)
		rect_nans = isnan(Z_rect);
		if (any(rect_nans(:)))				% Check if we have NaNs inside the Z_rect rectangle
			Z_filt(rect_nans) = NaN;
		end
	end

% -----------------------------------------------------------------------
function mask_b = mask_border(mask,w)
% Takes the MASK matrix which has ones defining the ROI and computes
% MASK_B that is composed of ones at the border of the ROI and W nodes
% to the inside

	[m,n] = size(mask);
	mask_b = [false(m,w) mask(:,1:end-w)];		% shift mask rigth
	mask_in = mask & mask_b;
	mask_b = [mask(:,w+1:end) false(m,w)];		% shift mask left
	mask_in = mask_in & mask_b;
	mask_b = [mask(w+1:end,:); false(w,n)];		% shift mask up
	mask_in = mask_in & mask_b;
	mask_b = [false(w,n); mask(1:end-w,:)];		% shift mask down
	mask_in = mask_in & mask_b;
	mask_b = xor(mask_in,mask);

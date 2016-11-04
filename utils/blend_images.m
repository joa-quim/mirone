function [img, inc, msg]  = blend_images(img_a, img_b, alpha_a, alpha_b, BB_a, BB_b, master, transparency)
% Blend two images of potential different resolutions and convering different areas BUT with an area in common.
% MASTER  -> A text string with either '1' (or 'A') or '2' (or 'B'). '1' IMG_A is the master image and it is the
%            IMG_B image that will be reinterpolated (if needed) to fir IMG_A resolution. '2' is obviously the vice-versa
%
% ALPHA_A -> An alpha mask with same size of IMG_A or empty ([]) if IMG_A has no transparency mask zone.
% ALPHA_B -> Same ofr IMG_B.
%
%            When master image is IMG_A and ALHPHA_B is not empty, the holes is IMG_B defined in ALHPHA_B
%            are filled with the data from IMG_A.
%
% BB_A|B  -> Bounding box of the A and B images ([x_min x_max y_min y_max]).
%
% TRANSPARENCY -> A value in the [0-1] interval to decide how blending is done on the common zone of the
%                 two images that is common to slave image alpha mask. Default = 0, means one smashes the other. 
%
% IMG    -> Returns the slave image reinterpolated to master image resolution and its transparent zones
%           filled from the data of master image.
% INC    -> New resolution ([inc_x inc_y]) of the resized master image.
% MSG    -> An error message in case things when wrong. In that case IMG = INC = []

	if (nargin < 8),	transparency = 0;	end
	msg = '';	img = [];	inc = [];
	if ((ndims(img_a) == 2 && ndims(img_b) == 3) || (ndims(img_a) == 3 && ndims(img_b) == 2))
		msg = 'The two images must both be indexed OR true color but NOT mixed types';
		return
	end

	if (master == 'A' || master == '1')
		[img, inc]  = blend_it(img_a, img_b, alpha_b, BB_a, BB_b, transparency);
	else
		[img, inc]  = blend_it(img_b, img_a, alpha_a, BB_b, BB_a, transparency);
	end

function [img_b, inc, msg] = blend_it(img_a, img_b, alpha_b, BB_a, BB_b, transparency)
% Resize IMG_B and implant it in the common area that it shares with IMG_A
% In areas where ALPHA_B is transparent put the IMG_A data.
	msg = '';	inc = [];
	P1.x = [BB_a(1) BB_a(1) BB_a(2) BB_a(2) BB_a(1)];	P1.hole = 0;
	P1.y = [BB_a(3) BB_a(4) BB_a(4) BB_a(3) BB_a(3)];
	P2.x = [BB_b(1) BB_b(1) BB_b(2) BB_b(2) BB_b(1)];	P2.hole = 0;
	P2.y = [BB_b(3) BB_b(4) BB_b(4) BB_b(3) BB_b(3)];
	P3 = PolygonClip(P1, P2, 1);				% Intersection of the two rectangles
	if (isempty(P3)),	msg = 'The two images do not overlap.';		return,		end
	rx_min = min(P3.x);			rx_max = max(P3.x);		ry_min = min(P3.y);			ry_max = max(P3.y);
	rect_crop = [rx_min ry_min rx_max-rx_min ry_max-ry_min];

	% Resize img_b to fit the img_a resolution
	inc_x = diff(BB_b(1:2)) / (size(img_a,2) - 1);
	inc_y = diff(BB_b(3:4)) / (size(img_a,1) - 1);
	inc = double([inc_x inc_y]);
	X_b = round(diff(BB_b(1:2)) / inc_x) + 1;
	Y_b = round(diff(BB_b(3:4)) / inc_y) + 1;
	img_b = cvlib_mex('resize', img_b, [Y_b X_b], 'bicubic');			% Resize IMG_B to the resolution of IMG_A

	[r_c] = cropimg(BB_b(1:2), BB_b(3:4), img_b, rect_crop, 'out_ind');		% Get new IMG_B indices of the intersection rectangle
	[I,P1.hole] = cropimg(BB_a(1:2),BB_a(3:4), img_a, rect_crop,'out_grid');% Crop IMG_A to common rectangle. P1.hole to shut up MLint
	img_a = cvlib_mex('resize',I,[diff(r_c(1:2)) diff(r_c(3:4))]+1,'bicubic');	% Don't remember anymore why, but it's right.
	if (~isempty(alpha_b) && numel(alpha_b) > 1)
		alpha_b = logical(cvlib_mex('resize', alpha_b, [Y_b X_b], 'bilinear'));	% Resize also the alpha mask
		%alpha_b = cvlib_mex('resize',alpha_b,[diff(r_c(1:2)) diff(r_c(3:4))]+1,'bilinear');
		alpha_b = alpha_b(r_c(1):r_c(2),r_c(3):r_c(4));
	end

	% Actually never tested this case.
	if (transparency)
		tmp = img_b(r_c(1):r_c(2),r_c(3):r_c(4), 1:end);
		cvlib_mex('addweighted',img_a,(1 - transparency), tmp, transparency)	% In-place
	end

	if (~isempty(alpha_b))
		for (k = 1:size(img_a,3))
			tmp = img_a(:,:,k);				tmp2 = img_b(r_c(1):r_c(2),r_c(3):r_c(4),k);
			tmp(alpha_b) = tmp2(alpha_b);	img_a(:,:,k) = tmp;
		end
	end
	img_b(r_c(1):r_c(2),r_c(3):r_c(4), 1:end) = img_a;


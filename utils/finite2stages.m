function stages = finite2stages(lon, lat, omega, t_start, half, side)
% Convert finite rotations to backwards stage rotations for backtracking
%
% LON, LAT, OMEGA & T_START are the finite rotation Euler pole parameters and age of pole
% Alternatively LON may be a Mx4 matrix with columns LON, LAT, OMEGA & T_START
% STAGES is a Mx5 matrix of stage pole (Euler) with the following format:
% lon(deg)  lat(deg)  tstart(Ma)  tstop(Ma)  ccw-angle(deg)
% stage records go from oldest to youngest rotation
%
% HALF = 1|2 If == 1 full angles are returned (good for plate reconstructions).
%            Else (== 2) compute half angles (good for flow lines in a single plate)
%
% NOTE: the notation is that the finite pole is b_ROT_a - Where B is the fixed plate
% The signal of HALF is used to compute b_STAGE_a (default) or a_STAGE_b (if HALF < 0)
%
% SIDE = 1  -> poles in the northern hemisphere
% SIDE = -1 -> poles in the southern hemisphere
% SIDE = 0  -> report positive rotation angles
%
% The libspotter always report the angles going from Old to New but we use in telha_m and other places
% the assumption is that we go from Young to Old. That means the angles are symetric of what we need.
% I'll set it here for the time being as -1 coefficient but in future this should be an input arg.
%
% Translated from C code of libspotter (Paul Wessel - GMT) (Now LGPL)
% Joaquim Luis 21-4-2005

%	Coffeeright (c) 2004-2013 by J. Luis
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

% $Id$

	n_args = nargin;
	if (~(n_args == 1 || n_args == 3 || n_args == 6))
		error('Wrong number of arguments')
	elseif (n_args == 1 || n_args == 3)
		if (n_args == 3),       half = lat;     side = omega;
		else                    half = 2;       side = 1;		% Default to half angles & North hemisphere poles
		end
		t_start = lon(:,4);     omega = lon(:,3);
		lat = lon(:,2);         lon = lon(:,1);
	end

	t_old = 0;
	R_young = eye(3);
	elon = zeros(1,length(lon));    elat = elon;    ew = elon;  t_stop = elon;
	for i = 1:length(lon)
		R_old = make_rot_matrix (lon(i), lat(i), omega(i));     % Get rotation matrix from pole and angle
		if (half > 0)                                           % the stages come in the reference b_STAGE_a
			R_stage = R_young * R_old;                          % This is R_stage = R_young^t * R_old
		else                                                    % the stages come in the reference a_STAGE_b
			R_stage = R_old * R_young;                          % This is R_stage = R_old * R_young^t
			R_stage = R_stage';
		end
		[elon(i), elat(i), ew(i)] = matrix_to_pole(R_stage,side); % Get rotation parameters from matrix
		if (elon(i) > 180), elon(i) = elon(i) - 360;     end    % Adjust lon
		R_young = R_old';                                       % Sets R_young = transpose (R_old) for next round
		t_stop(i) = t_old;
		t_old = t_start(i);
	end

	% Flip order since stages go from oldest to youngest
	% And for the reason explained above in the man section. Change the sing of the output angles
	ang_sign = -1;
	stages = flipud([elon(:) elat(:) t_start(:) t_stop(:) ew(:) / abs(half) * ang_sign]);

% --------------------------------------------------------
function R = make_rot_matrix (lonp, latp, w)
% lonp, latp	Euler pole in degrees
% w		angular rotation in degrees
% R		the rotation matrix

	D2R = pi / 180;
	[E0,E1,E2] = sph2cart(lonp*D2R,latp*D2R,1);

	sin_w = sin(w * D2R);
	cos_w = cos(w * D2R);
	c = 1 - cos_w;

	E_x = E0 * sin_w;
	E_y = E1 * sin_w;
	E_z = E2 * sin_w;
	E_12c = E0 * E1 * c;
	E_13c = E0 * E2 * c;
	E_23c = E1 * E2 * c;

	R(1,1) = E0 * E0 * c + cos_w;
	R(1,2) = E_12c - E_z;
	R(1,3) = E_13c + E_y;

	R(2,1) = E_12c + E_z;
	R(2,2) = E1 * E1 * c + cos_w;
	R(2,3) = E_23c - E_x;

	R(3,1) = E_13c - E_y;
	R(3,2) = E_23c + E_x;
	R(3,3) = E2 * E2 * c + cos_w;

% --------------------------------------------------------
function [plon,plat,w] = matrix_to_pole (T,side)
	R2D = 180 / pi;
	T13_m_T31 = T(1,3) - T(3,1);
	T32_m_T23 = T(3,2) - T(2,3);
	T21_m_T12 = T(2,1) - T(1,2);
	H = T32_m_T23 * T32_m_T23 + T13_m_T31 * T13_m_T31;
	L = sqrt (H + T21_m_T12 * T21_m_T12);
	H = sqrt (H);
	tr = T(1,1) + T(2,2) + T(3,3);

	plon = atan2(T13_m_T31, T32_m_T23) * R2D;
	%if (plon < 0)     plon = plon + 360;  end
	plat = atan2(T21_m_T12, H) * R2D;
	w = atan2(L, (tr - 1)) * R2D;

	if ((side == 1 && plat < 0) || (side == -1 && plat > 0))
		plat = -plat;
		plon = plon + 180;
		if (plon > 360),    plon = plon - 360;  end
		w = -w;
	end

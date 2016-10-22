function line_of_sight(hand, observerAltitude)
% Compute the local horizon from a DEM.
%
% The Terrain elevation (in degrees) functions are from FEX: 59421 - 'sunMask' (BSD License)
% from Benjamin Pillot. The rest, namely the part that computes the horizon polygon,
% make use of several Mirone functions. Namely, the Ellipsoide model used is the one
% stored in handles.DefineEllipsoide
%

%	Coffeeright (c) 2004-2016 by J. Luis
%
% License: BSD
%
%	Contact info: w3.ualg.pt/~jluis/mirone
% --------------------------------------------------------------------

% $Id: $

% This is a non-finished code. The computations works fine but the API needs to be decided.

	handles = guidata(hand);
	D2R = pi / 180;
	npt_circ = 360;
	rad = 0.449371337890625;		% ~ 50 km
	observerAltitude = 2;
	[X,Y,Z] = load_grd(handles);
	if isempty(Z),		return,		end		% An error message was already issued
	x_lim = get(handles.axes1,'xlim');    y_lim = get(handles.axes1,'ylim');

	% Get the observer point and move it to bin center (on top of grid node)
	x_click = get(hand,'XData');	y_click = get(hand,'YData');
	col = aux_funs('getPixel_coords', size(Z,2), [X(1) X(end)], x_click);
	row = aux_funs('getPixel_coords', size(Z,1), [Y(1) Y(end)], y_click);
	x_click = X(col);		y_click = Y(row);

	%rad = (3.57 * sqrt(Z(row, col) + observerAltitude) + 1) / 6371;		% Approximate horizon at observer altitude

	pad_x = handles.head(8);		pad_y = handles.head(9);
	[latc, lonc] = circ_geo(y_click, x_click, rad, [], npt_circ, handles.geog == 2);
	x_min_c = max(min(lonc)-pad_x, x_lim(1));		x_max_c = min(max(lonc)+pad_x, x_lim(2));
	y_min_c = max(min(latc)-pad_y, y_lim(1));		y_max_c = min(max(latc)+pad_y, y_lim(2));

	% Cut the smallest area need for computation
	[X_rec, Y_rec, Z_rec, head_rec] = mirone('ImageCrop_CB', handles, [x_min_c x_max_c y_min_c y_max_c], 'CropaGrid_pure');
 	Z_rec = double(Z_rec);

	% Geo-referencing meshgrid
	[longitude, latitude] = meshgrid(X_rec, Y_rec);

	[zObs, row, col] = bi_linear(X_rec, Y_rec, Z_rec, x_click, y_click);
	zObs = Z_rec(row, col) + observerAltitude;
	Z_rec(row, col) = Z_rec(row, col) + observerAltitude;

	% Corresponding elevation angle
	alpha = getElevationAngle(zObs, Z_rec, y_click * D2R, latitude * D2R, x_click * D2R, longitude * D2R, handles.DefineEllipsoide);
	alpha(row, col) = NaN;

	tic
	alpha = single(alpha);		% grid_profiler on singles is mutch faster
	hand_fake = struct('X',X_rec, 'Y',Y_rec, 'Z',alpha, 'head',head_rec, 'figure1',handles.figure1);
	x = zeros(npt_circ,1);		y = zeros(npt_circ,1);
	for (k = 1:npt_circ)
		[xx, yy, zz] = grid_profiler(hand_fake, [x_click; lonc(k)], [y_click; latc(k)], false);
		[ma, ind] = max(zz);
		x(k) = xx(ind);		y(k) = yy(ind);
	end
	toc

	h = line('Parent',handles.axes1, 'XData', x, 'YData', y);
	draw_funs(h,'line_uicontext')
	
% 	h = mirone(X_rec, Y_rec, alpha);
% 	handles2 = guidata(h);
% 	h = line('Parent',handles2.axes1, 'XData', x, 'YData', y);
% 	draw_funs(h,'line_uicontext')
% 	h = line('Parent',handles.axes1, 'XData', xx, 'YData', yy);
% 	draw_funs(h,'line_uicontext')

	c = sin(y_click*D2R).*sin(y*D2R) + cos(y_click*D2R).*cos(y*D2R).*cos((x-x_click)*D2R);
	ang = max(acos(c) / D2R);
	if ((rad - ang) < min(handles.head(8:9)))
		warndlg('Warning: Limiting radius might have been too small. Horizon was probably clipped by it.', 'Warning')
	end

% -------------------------------------------------------------------------------------------
function alpha = getElevationAngle(h_A, h_B, latitude_A, latitude_B, longitude_A, longitude_B, ellipsoid)
% Compute cartesian coordinates of point A and B located at altitude
% h_A and h_B from the ellipsoid surface (ellipsoidal heights)
	[x_A, y_A, z_A] = geographic2cartesian(latitude_A, longitude_A, h_A, ellipsoid);
	[x_B, y_B, z_B] = geographic2cartesian(latitude_B, longitude_B, h_B, ellipsoid);

	x_B = x_B - x_A;	y_B = y_B - y_A;	z_B = z_B - z_A;
	% Scalar product between AB and normal to the point A
	innerProduct = x_B .* cos(longitude_A) .* cos(latitude_A) + y_B .* ...
		sin(longitude_A) .* cos(latitude_A) + z_B .* sin(latitude_A);

	% Angular elevation computation
	norm = sqrt(x_B .^2 + y_B .^2 + z_B .^2);
	alpha = asin(innerProduct./norm) * 180/pi;

% -------------------------------------------------------------------------------------------
function [x, y, z] = geographic2cartesian_(latitude, longitude, altitude, ellipsoid)
% Compute semi-major axis and eccentricity of the specified ellipsoid
% Devec version that is still slower than the (ultageous memory consumig) vectorized one
	a = ellipsoid(1);			% semiMajorAxis
	e2 = 2*ellipsoid(3) - ellipsoid(3)*ellipsoid(3);		% e^2 = 2f - f^2

	% Compute ellipsoid normal
	N = zeros(size(altitude,1),1);
	for (m = 1:size(altitude,1))
		N(m) = a./sqrt(1 - e2 * (sin(latitude(m,1))).^2);
	end

	% Compute cartesian coordinates from geographic coordinates
	x = zeros(size(altitude));	y = zeros(size(altitude));	z = zeros(size(altitude));	
	cos_lon = cos(longitude(1,:));		sin_lon = sin(longitude(1,:));
	cos_lat = cos(latitude(:,1));		sin_lat = sin(latitude(:,1));
	for (n = 1:size(altitude,2))
		for (m = 1:size(altitude,1))
			x(m,n) = (N(m) + altitude(m,n)) * cos_lon(n) * cos_lat(m);
			y(m,n) = (N(m) + altitude(m,n)) * sin_lon(n) * cos_lat(m);
			z(m,n) = (N(m) * (1 - e2) + altitude(m,n))  * sin_lat(m);
		end
	end

% -------------------------------------------------------------------------------------------
function [x, y, z] = geographic2cartesian(latitude, longitude, altitude, ellipsoid)
% Compute semi-major axis and eccentricity of the specified ellipsoid
	a = ellipsoid(1);			% semiMajorAxis
	e2 = 2*ellipsoid(3) - ellipsoid(3)*ellipsoid(3);		% e^2 = 2f - f^2

	% Compute ellipsoid normal
	N = a./sqrt(1 - e2 * (sin(latitude)).^2);

	% Compute cartesian coordinates from geographic coordinates
	x = (N + altitude) .* cos(longitude) .* cos(latitude);
	y = (N + altitude) .* sin(longitude) .* cos(latitude);
	z = (N * (1 - e2) + altitude) .* sin(latitude);

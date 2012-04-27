function area = area_geo(lat,lon,ellipsoid)
% AREA_GEO  Calculates the spherical surface area of a polygon
%
%  area = AREA_GEO(lat,lon, ellips) calculates the spherical surface area of the polygon
%  specified by the input vectors lat, lon. The calculation uses a line integral
%  approach. The output is the surface area fraction covered by the polygon on a
%  unit sphere, or on a ellipsoide is a third input arg is transmitted.
%  ELLIPSOIDE (optional) is either a [a eccentricity] or [a b invflattening] vector
%  The second case (a 1x3 vector) the semi-minor axis 'b' value is ignored but this
%  form is practical as it is what is used in handles.DefineEllipsoide.
%
%  Accuracy of the integration method is inversely proportional to
%  the distance between lat/lon points.
%
%  This function uses the same algorithm as Matlab's AREAINT

%	Coffeeright (c) 2004-2012 by J. Luis
%
%	Contact info: w3.ualg.pt/~jluis/mirone
% --------------------------------------------------------------------

	D2R = pi/180;
	bodyRadius = 0;
	lat  = lat(:) * D2R;   lon = lon(:) * D2R;			% Make sure they are column vectors
	if (lat(1) ~= lat(end) || lon(1) ~= lon(end))		% Make sure that the polygon is closed.
		lat = [lat; lat(1)];   lon = [lon; lon(1)];
	end

	if (nargin == 3)				% Third argin must be a 2 or 3 elements vector, but not fully tested.
		a = ellipsoid(1);
		if (numel(ellipsoid) == 2)
			ecc = ellipsoid(2);
		elseif (numel(ellipsoid) == 3)
			f = ellipsoid(3);		% reciprocal flattening
			ecc = sqrt(2*f - f^2);	% first eccentricity
		end
		if ecc > 0
			f1 = a^2 / (2);
			f2 = (1 - ecc^2) / (2*ecc);
			f3 = log((1+ecc) / (1-ecc));
			bodyRadius = sqrt(f1 * (1 + f2 * f3));
		else
			bodyRadius = a;
		end

		f1 = ecc^2 /3 + 31*ecc^4 /180 + 59*ecc^6 /560;
		f2 = 17*ecc^4 /360 + 61*ecc^6 /1260;
		f3 = 383*ecc^6 /45360;

		%  Truncated series expansion.
		lat = lat - f1*sin(2*lat) + f2*sin(4*lat) - f3*sin(6*lat);
	end

	origin = zeros(size(lat));                      % Set origin for integration. (0,0) is as good as any other
	% Get colatitude and azimuth 
	colat = acos(sin(origin).*sin(lat) + cos(origin).*cos(lat).*cos(lon-origin));
	azim  = azimuth_geo(origin,origin,lat/D2R,lon/D2R) * D2R;   % azimuth_geo deals with angles in degrees

	diff_azim = diff(azim);
	% truncate angles into the [-pi pi] range
	diff_azim = pi*((abs(diff_azim)/pi) - 2*ceil(((abs(diff_azim)/pi)-1)/2)) .* sign(diff_azim);

	% Determine average surface distance for each step
	diff_colat = diff(colat)/2;
	colat2 = colat(1:length(colat)-1) + diff_colat;

	tmp = (1-cos(colat2)).*diff_azim;   % Integral over azimuth is 1-cos(colatitude)
	area = abs(sum(tmp))/(4*pi);        % Integrate. Area is as a fraction of that of the unit sphere.

	if (bodyRadius)						% If we know thw body geom, compute the surface area.
		area = area * 4 * pi * (bodyRadius ^2);
	end

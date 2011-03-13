function area = area_geo(lat,lon)
% AREA_GEO  Calculates the spherical surface area of a polygon
%
%  area = AREA_GEO(lat,lon) calculates the spherical surface area of the polygon
%  specified by the input vectors lat, lon. The calculation uses a line integral
%  approach. The output is the surface area fraction covered by the polygon on a
%  unit sphere. Accuracy of the integration method is inversely proportional to
%  the distance between lat/lon points.
%
%  This function uses the same algorithm as Matlab's AREAINT

%	Copyright (c) 2004-2011 by J. Luis
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

D2R = pi/180;
lat = lat(:);       lon = lon(:);               % Make sure that lat and lon are column vectors
lat  = lat * D2R;   lon = lon * D2R;
lat = [lat; lat(1)];   lon = [lon; lon(1)];     % Make sure that the polygon is closed.
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

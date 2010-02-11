function [lon_s,lat_s,ang_s] = add_poles(lon1,lat1,ang1,lon2,lat2,ang2)
% Add two finite rotation Euler poles.
% All angles are assumed to be in degrees

%	Copyright (c) 2004-2006 by J. Luis
%
%	This program is free software; you can redistribute it and/or modify
%	it under the terms of the GNU General Public License as published by
%	the Free Software Foundation; version 2 of the License.
%
%	This program is distributed in the hope that it will be useful,
%	but WITHOUT ANY WARRANTY; without even the implied warranty of
%	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%	GNU General Public License for more details.
%
%	Contact info: w3.ualg.pt/~jluis/mirone
% --------------------------------------------------------------------

D2R = pi / 180;
lon1 = lon1 * D2R;      lat1 = lat1 * D2R;      ang1 = ang1 * D2R;
lon2 = lon2 * D2R;      lat2 = lat2 * D2R;      ang2 = ang2 * D2R;

% Build the first rotation matrix
[px,py,pz] = sph2cart(lon1,lat1,1);               % Convert first pole to cartesian
num = 1.0 - cos(ang1);
r1(1,1) = px * px * num + cos(ang1);
r1(1,2) = px * py * num - pz * sin(ang1);
r1(1,3) = px * pz * num + py * sin(ang1);

r1(2,1) = py * px * num + pz * sin(ang1);
r1(2,2) = py * py * num + cos(ang1);
r1(2,3) = py * pz * num - px * sin(ang1);

r1(3,1) = pz * px * num - py * sin(ang1);
r1(3,2) = pz * py * num + px * sin(ang1);
r1(3,3) = pz * pz * num + cos(ang1);


% Build the second rotation matrix
[px,py,pz] = sph2cart(lon2,lat2,1);               % Convert the second pole to cartesian
num = 1.0 - cos(ang2);
r2(1,1) = px * px * num + cos(ang2);
r2(1,2) = px * py * num - pz * sin(ang2);
r2(1,3) = px * pz * num + py * sin(ang2);

r2(2,1) = py * px * num + pz * sin(ang2);
r2(2,2) = py * py * num + cos(ang2);
r2(2,3) = py * pz * num - px * sin(ang2);

r2(3,1) = pz * px * num - py * sin(ang2);
r2(3,2) = pz * py * num + px * sin(ang2);
r2(3,3) = pz * pz * num + cos(ang2);

% Compute the total rotation matrix T = R2 * R1
T = r2 * r1;

% Compute the Euler pole whose rotation matrix is T
lon_s = atan2( (T(1,3) - T(3,1)), (T(3,2) - T(2,3)) ) / D2R;        % Lon
tmp = sqrt( (T(3,2)-T(2,3))^2 + (T(1,3)-T(3,1))^2 + (T(2,1)-T(1,2))^2);
lat_s = asin( (T(2,1) - T(1,2)) / tmp) / D2R;                      % Lat
ang_s = atan(tmp / (T(1,1)+T(2,2)+T(3,3) -1 )) / D2R;

% Make sure that 0 <= ang_s <= 180
if (ang_s < 0),     ang_s = ang_s + 180;    end

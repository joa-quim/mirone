function [plon,plat,omega] = calcBoninEulerPole(xx,yy,x2,y2)
% Computes an Euler pole using the simplified (AKA Bonin) method.
%   calcBoninEulerPole(xx,yy) is used with an "Euler trapezium" (length(xx) MUST be = 5)
%   calcBoninEulerPole(xx,yy,x2,y2) computes the pole that brings together the
%   lines defined by xx,yy and x2,y2. NOTE first and second line must be given in the same order.
% If no output is requested the pole is written on the command line

plon = NaN;     plat = NaN;     omega = NaN;
if (nargin == 4)
    % No testing that xx, yy, & x2, y2 are inded the coords of two lines
    xx = [xx(:)' x2(2) x2(1)];      % Do this to reuse the trapezium algo
    yy = [yy(:)' y2(2) y2(1)];
else
    % Called by "Euler trapezium" so length(xx) MUST be = 5
    if (length(xx) ~= 5),    return;     end
end

D2R = pi / 180;
sa1 = sin(yy(1)*D2R);    sa2 = sin(yy(4)*D2R);
ca1 = cos(yy(1)*D2R);    ca2 = cos(yy(4)*D2R);
so1 = sin(xx(1)*D2R);    so2 = sin(xx(4)*D2R);
co1 = cos(xx(1)*D2R);    co2 = cos(xx(4)*D2R);
c_teta1 = sa1*sa2 + ca1*ca2*cos((xx(4)-xx(1))*D2R);
teta1 = acos(c_teta1);
sf = 1/sin(teta1);
x = (ca1*so1*sa2 - ca2*so2*sa1) * sf;       y = (ca2*co2*sa1 - ca1*co1*sa2) * sf;
z = (ca1*co1*ca2*so2 - ca2*co2*ca1*so1) * sf;
lon = atan2(y,x);
lat = atan(z/sqrt(x.^2 + y.^2));
[px,py,pz] = sph2cart(lon,lat,1);               % Convert pole to cartesian
[ax,ay,az] = sph2cart(xx(1)*D2R,yy(1)*D2R,1);   % Convert first point to cartesian
a = [ax;ay;az];
[ax,ay,az] = sph2cart(xx(2)*D2R,yy(2)*D2R,1);   % Convert second point to cartesian
b = [ax;ay;az];

% Build the rotation matrix
num = 1.0 - cos(teta1);
r(1,1) = px * px * num + cos(teta1);
r(1,2) = px * py * num - pz * sin(teta1);
r(1,3) = px * pz * num + py * sin(teta1);

r(2,1) = py * px * num + pz * sin(teta1);
r(2,2) = py * py * num + cos(teta1);
r(2,3) = py * pz * num - px * sin(teta1);

r(3,1) = pz * px * num - py * sin(teta1);
r(3,2) = pz * py * num + px * sin(teta1);
r(3,3) = pz * pz * num + cos(teta1);

% Rotate the side deffined by the first 2 points about the just computed 1st pole
a_rotated = r * a;
b_rotated = r * b;

% Get the rotated side 1 back into spherical coords
[lon_ar,lat_ar] = cart2sph(a_rotated(1),a_rotated(2),a_rotated(3));
[lon_br,lat_br] = cart2sph(b_rotated(1),b_rotated(2),b_rotated(3));

% Compute the new angle between the rotated first side and the original second side

% Make the 4th vertice be the new North Pole
orig = new_Npole(yy(4),xx(4));
[lat_brPN,lon_brPN] = rotate(lat_br/D2R,lon_br/D2R,orig);
[yy(3),xx(3)] = rotate(yy(3),xx(3),orig);
teta2 = (xx(3) - lon_brPN) * D2R;

[px,py,pz] = sph2cart(xx(4)*D2R,yy(4)*D2R,1);               % Convert the second pole to cartesian
% Build the rotation matrix
num = 1.0 - cos(teta2);
r2(1,1) = px * px * num + cos(teta2);
r2(1,2) = px * py * num - pz * sin(teta2);
r2(1,3) = px * pz * num + py * sin(teta2);

r2(2,1) = py * px * num + pz * sin(teta2);
r2(2,2) = py * py * num + cos(teta2);
r2(2,3) = py * pz * num - px * sin(teta2);

r2(3,1) = pz * px * num - py * sin(teta2);
r2(3,2) = pz * py * num + px * sin(teta2);
r2(3,3) = pz * pz * num + cos(teta2);

% Compute the total rotation matrix T = R2 * R
T = r2 * r;

% Compute the Euler pole whose rotation matrix is T
%plon = atan( (T(1,3) - T(3,1)) / (T(3,2) - T(2,3)) ) / D2R;        % Lon
plon = atan2( (T(1,3) - T(3,1)), (T(3,2) - T(2,3)) ) / D2R;        % Lon
tmp = sqrt( (T(3,2)-T(2,3))^2 + (T(1,3)-T(3,1))^2 + (T(2,1)-T(1,2))^2);
plat = asin( (T(2,1) - T(1,2)) / tmp) / D2R;                      % Lat
omega = atan(tmp / (T(1,1)+T(2,2)+T(3,3) -1 )) / D2R;

% Make sure that 0 <= omega <= 180
if (omega < 0),     omega = omega + 180;    end

% Now add the two rotations
% disp(['FIRST ROT -> Lon = ' num2str(lon/D2R) '  Lat = ' num2str(lat/D2R) '  Omega = ' num2str(teta1/D2R)])
% disp(['SECOND ROT -> Lon = ' num2str(xx(4)) '  Lat = ' num2str(yy(4)) '  Omega = ' num2str(teta2/D2R)])
if (nargout == 0)
    disp(['Lon = ' num2str(plon) '  Lat = ' num2str(plat) '  Omega = ' num2str(omega)])
end

%--------------------------------------------------------------------------------------------------
function origin = new_Npole(polelat,polelon)
%NEW_NPOLE  Computes the origin vector to place a point at the North pole (angles are in degrees)
% Copyright 1996-2003 The MathWorks, Inc.

%  Transform input data to radians
D2R = pi / 180;
polelat = polelat(:) * D2R;
polelon = polelon(:) * D2R;

%  Get the indices for the northern and southern hemisphere new poles
ind1 = find(polelat >= 0);    ind2 = find(polelat <  0);

origlat = zeros(size(polelat));
origlon = zeros(size(polelon));
orient  = zeros(size(polelat));

%  Compute the origin for northern hemisphere poles
if ~isempty(ind1)
    origlat(ind1) = pi/2 - polelat(ind1);
    origlon(ind1) = pi*((abs(polelon(ind1)+pi)/pi) - ...        % Make sure that angles are in [-180 180] range
        2*ceil(((abs(polelon(ind1)+pi)/pi)-1)/2)) .* sign(polelon(ind1)+pi);
    ind3 = find(polelat == pi/2);    %  Correct for any poles staying at the north pole
    if ~isempty(ind3)
        origlon(ind3) = pi*((abs(polelon(ind3)+pi)/pi) - ...    % Make sure that angles are in [-180 180] range
            2*ceil(((abs(polelon(ind3)+pi)/pi)-1)/2)) .* sign(polelon(ind3)+pi);
    end
end

%  Compute the origin for southern hemisphere poles
if ~isempty(ind2)
    origlat(ind2) = pi/2 + polelat(ind2);
    origlon(ind2) = pi*((abs(polelon(ind2)+pi)/pi) - ...        % Make sure that angles are in [-180 180] range
        2*ceil(((abs(polelon(ind2)+pi)/pi)-1)/2)) .* sign(polelon(ind2)+pi);
	orient(ind2)  = -pi;
end

origin = [origlat origlon orient] / D2R;    %  Transform back to degrees

%--------------------------------------------------------------------------------------------------
function [lat1,lon1] = rotate(lat,lon,orig)
%ROTATE  Rotate data for specified orig and orientation (angles are in degrees)
%  Copyright 1996-2003 The MathWorks, Inc.
	
D2R = pi / 180;
% Convert to radians
lat = lat * D2R;    lon = lon * D2R;    orig = orig * D2R;

rot1 = [cos(orig(2)) sin(orig(2))  0        % Rotation matrix about x axis
       -sin(orig(2)) cos(orig(2))  0
	    0            0             1];
rot2 = [cos(orig(1)) 0 sin(orig(1))         % Rotation matrix about y axis
        0            1 0
	   -sin(orig(1)) 0 cos(orig(1))];
rot3 = [1  0            0
        0  cos(orig(3)) sin(orig(3))        % Rotation matrix about z axis
        0 -sin(orig(3)) cos(orig(3))];

rot = rot3 * rot2 * rot1;                   % Euler rotation matrix

%  Move pi/2 points epsilon inward to prevent round-off problems with pi/2 points.
epsilon = 1e-5 * D2R;
indx = find(abs(pi/2 - abs(lat)) <= epsilon);
if ~isempty(indx)
	lat(indx) = (pi/2 - epsilon) * sign(lat(indx));
end

%  Prevent possible confusion with points at +180 or -180 degrees
lon = atan2(sin(lon*(1 - 1e-6)),cos(lon*(1 - 1e-6)));

%  Compute the new x,y,z point in cartesian space
xyz = rot * ([cos(lat).*cos(lon) cos(lat).*sin(lon) sin(lat)]');

% epsilon = 1.0E-8;
% indx = find(abs(xyz(1)) <= epsilon & abs(xyz(2)) <= epsilon);   % Be careful with x & y ~= 0 in atan2
% if ~isempty(indx);   x(indx) = 0;  y(indx) = 0;   end

[lon1, lat1] = cart2sph(xyz(1),xyz(2),xyz(3));  % Transform to spherical coordinates
lat1 = lat1 / D2R;      lon1 = lon1 / D2R;      % Transform back to degrees

function s = vdist(lat1,lon1,lat2,lon2,ellipsoide)
% VDIST - compute distance between points on the WGS-84 ellipsoidal Earth
%         to within a few millimeters of accuracy using Vincenty's algorithm
%
% s = vdist(lat1,lon1,lat2,lon2)
%
% s = distance in meters
% lat1 = GEODETIC latititude of first point (degrees)
% lon1 = longitude of first point (degrees)
% lat2, lon2 = second point (degrees)
%
% s = vdist(lat1,lon1,lat2,lon2,ellipsoide)
%  Transmits ellipsoide deffinition (either as [a,b] or [a,f]) as fifth argument ELLIPSOIDE
%
%  Original algorithm source:
%  T. Vincenty, "Direct and Inverse Solutions of Geodesics on the Ellipsoid
%  with Application of Nested Equations", Survey Review, vol. 23, no. 176,
%  April 1975, pp 88-93
%
% Notes: (1) Error correcting code, convergence failure traps, antipodal corrections,
%            polar error corrections, WGS84 ellipsoid perameters, testing, and comments
%            written by Michael Kleder, 2004.
%        (2) Vincenty describes his original algorithm as precise to within
%            0.01 millimeters, subject to the ellipsoidal model.
%        (3) Essentially anitpodal points are treated as exactly antipodal,
%            potentially reducing accuracy by a small amount.
%        (4) Failures for points exactly at the poles are eliminated by
%            moving the points by 0.6 millimeters.
%        (5) Vincenty's azimuth formulas are not implemented in this
%            version, but are available as comments in the code.
%        (6) The Vincenty procedure was transcribed verbatim by Peter Cederholm,
%            August 12, 2003. It was modified and translated to English by Michael Kleder.
%            Mr. Cederholm's website is http://www.plan.aau.dk/~pce/

% Input check:
if ( (abs(lat1-lat2) < 1e-7) && (abs(lon1-lon2) < 1e-7) )
    s = 0;  return
end

if (abs(lat1) > 90 || abs(lat2) > 90)
    errordlg('Input latitudes must be between -90 and 90 degrees, inclusive.','Error');
    s = 0;  return
end
if (nargin == 5)    % An ellipsoide vector (with a & b OR a & f) was supplyied
    a = ellipsoide(1);    %b = ellipsoide(2);
    if (ellipsoide(2) < 1)      % Second ellipsoide argument contains flattening instead of minor axis
        f = ellipsoide(2);  b = a * (1 - f);
    else                        % Second ellipsoide argument contains minor axis
        f = (a-ellipsoide(2))/a;
    end
else
    % Supply WGS84 earth ellipsoid axis lengths in meters:
    a = 6378137; % definitionally
    b = 6356752.31424518; % computed from WGS84 earth flattening coefficient definition
    f = (a-b)/a;
end
% convert inputs in degrees to radians:
lat1 = lat1 * 0.0174532925199433;
lon1 = lon1 * 0.0174532925199433;
lat2 = lat2 * 0.0174532925199433;
lon2 = lon2 * 0.0174532925199433;
% correct for errors at exact poles by adjusting 0.6 millimeters:
if (abs(pi/2-abs(lat1)) < 1e-10)
    lat1 = sign(lat1)*(pi/2-(1e-10));
end
if (abs(pi/2-abs(lat2)) < 1e-10)
    lat2 = sign(lat2)*(pi/2-(1e-10));
end
%f = (a-b)/a;
U1 = atan((1-f)*tan(lat1));
U2 = atan((1-f)*tan(lat2));
lon1 = mod(lon1,2*pi);
lon2 = mod(lon2,2*pi);
L = abs(lon2-lon1);
if (L > pi),    L = 2*pi - L;    end
lambda = L;
lambdaold = 0;
itercount = 0;
while ~itercount || abs(lambda-lambdaold) > 1e-12  % force at least one execution
    itercount = itercount+1;
    if itercount > 50
        warning('Points are essentially antipodal. Precision may be reduced slightly.');
        lambda = pi;
        break
    end
    lambdaold = lambda;
    sinsigma = sqrt((cos(U2)*sin(lambda))^2+(cos(U1)*sin(U2)-sin(U1)*cos(U2)*cos(lambda))^2);
    cossigma = sin(U1)*sin(U2)+cos(U1)*cos(U2)*cos(lambda);
    sigma = atan2(sinsigma,cossigma);
    alpha = asin(cos(U1)*cos(U2)*sin(lambda)/sin(sigma));
    cos2sigmam = cos(sigma)-2*sin(U1)*sin(U2)/cos(alpha)^2;
    C = f/16*cos(alpha)^2*(4+f*(4-3*cos(alpha)^2));
    lambda = L+(1-C)*f*sin(alpha)*(sigma+C*sin(sigma)*...
        (cos2sigmam+C*cos(sigma)*(-1+2*cos2sigmam^2)));
    % correct for convergence failure in the case of essentially antipodal points
    if lambda > pi
        warndlg('Points are essentially antipodal. Precision may be reduced slightly.','Warning');
        lambda = pi;        break
    end
end
u2 = cos(alpha)^2*(a^2-b^2)/b^2;
A = 1+u2/16384*(4096+u2*(-768+u2*(320-175*u2)));
B = u2/1024*(256+u2*(-128+u2*(74-47*u2)));
deltasigma = B*sin(sigma)*(cos2sigmam+B/4*(cos(sigma)*(-1+2*cos2sigmam^2)...
    -B/6*cos2sigmam*(-3+4*sin(sigma)^2)*(-3+4*cos2sigmam^2)));
s = b*A*(sigma-deltasigma);

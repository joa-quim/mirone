function [az, az2] = azimuth_geo(lat1,lon1,lat2,lon2,unites)
% Calculates great circle azimuth on a sphere. 
% lon1, lat1, ... are in degrees (unless, see below) and may be scalars or column vectors. 
% If UNITES exists and == 'radians', lats & lons are assumed to be in radians and so will be az

% No error testing.

do_conv = 1;
if (nargin == 5 && strcmp(unites,'radians')),   do_conv = 0;	end
D2R = pi/180;
if (do_conv)
    lat1 = lat1*D2R;    lon1 = lon1*D2R;
    lat2 = lat2*D2R;    lon2 = lon2*D2R;
end

f1 = cos(lat2) .* sin(lon2-lon1);
f2 = cos(lat1) .* sin(lat2);
f3 = sin(lat1) .* cos(lat2) .* cos(lon2-lon1);
az = atan2(f1,f2-f3);

% Make sure that angles are in the [0 2pi] interval
az = pi*((abs(az)/pi) - 2*ceil(((abs(az)/pi)-1)/2)) .* sign(az);
epsilon = -1e-6*D2R;            %  Allow points near zero
indx = find(az < epsilon);

%  Shift the points in [-pi 0] to [pi 2pi] range
if ~isempty(indx);  az(indx) = az(indx) + 2*pi;  end;

indx = find(az < 0);            %  Reset near zero points
if (~isempty(indx));  az(indx) = zeros(size(indx));  end

if (do_conv)        % Convert back to degrees
    az = az / D2R;
end

if (nargout == 2),		az2 = az;	end		% For convenience of buffer_j()
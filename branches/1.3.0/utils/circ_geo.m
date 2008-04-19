function [latc,lonc] = circ_geo(lat,lon,rng,azim,np)
% Small circle defined by its center, range and azimuth
% circ_geo(lat,lon,rng,azim,np)
% azim and np are optional arguments. "azim" is a one or two-column vector. 
% For single column, returns the arc between 0 and azim. For two columns, returns
% the arc between azim(1) and azim(2).
% "np" specifies the number of output points [default = 180].
% All angles are in degrees.

D2R = pi/180;   npts  = 180;       az    = [];
if (nargin < 3)
    errordlg('Error calling circ_geo. Must give at least 3 arguments','Error')
elseif (nargin == 4)
    az = azim;
elseif (nargin == 5)
    az = azim;    npts = np;
end
        
%  Allow for multiple circles starting from the same point
if length(lat) == 1 && length(lon) == 1 && length(rng) ~= 0
    lat = lat(ones(size(rng)));   lon = lon(ones(size(rng)));
end

if isempty(az),    az = [0 360];    az = az(ones([size(lat,1) 1]), :);  end

% convert to radians
lat = lat * D2R;    lon = lon * D2R;    az = az * D2R;    rng = rng * D2R;

%  Expand the azimuth inputs
if size(az,2) == 1              % Single column azimuth inputs
    az_neg = zeros(size(az));    az_pos = az;
else                            % Two column azimuth inputs
    az_neg = az(:,1);            az_pos = az(:,2);
end

az = zeros([size(az_neg,1) npts]);
for i = 1:size(az,1)
    % Handle case where limits give more than half of the circle. Go clockwise from start to end.
	if az_neg(i) > az_pos(i);     az_pos(i) = az_pos(i)+2*pi;	end
	az(i,:) = linspace(az_neg(i),az_pos(i),real(npts));	
end

% Each circle occupies a row of the output matrices.
lat = lat(:,ones([1,size(az,2)]) );
lon = lon(:,ones([1,size(az,2)]) );
rng = rng(:,ones([1,size(az,2)]) );

%  Ensure correct azimuths at either pole.
epsilon = 1e-07;     % Set tolerance
ind = find(lat >= pi/2-epsilon);    az(ind) = pi;    % starting at north pole
ind = find(lat <= epsilon-pi/2);    az(ind) = 0;     % starting at south pole

temp1  = sin(lat).*cos(rng);        temp2  = cos(lat).*sin(rng).*cos(az);
latc   = asin(temp1+temp2);
temp1  = sin(rng).*sin(az);         temp2  = cos(lat).*cos(rng);
temp3  = sin(lat).*sin(rng).*cos(az);
lonc   = lon + atan2(temp1,temp2-temp3);

% Truncate angles into the [-pi pi] range
lonc = pi*((abs(lonc)/pi) - 2*ceil(((abs(lonc)/pi)-1)/2)) .* sign(lonc);

%  Convert back to degrees
latc = latc / D2R;      lonc = lonc / D2R;

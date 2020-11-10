function  varargout = map_funs(fun,varargin)
% ...

	[varargout{1:nargout}] = feval(fun, varargin{:});

% ---------------------------------------------------------------------
function [lat2,lon2] = polyjoin(latcells,loncells)
% POLYJOIN convert polygon segments from cell array to vector format
% 
% [lat2,lon2] = POLYJOIN(latcells,loncells) converts polygons from cell
% array format to vector format. In cell array format, each element of 
% the cell array is a separate polygon. Each polygon may consist of an 
% outer contour followed by holes separated with NaNs. In vector format, 
% each vector may contain multiple faces separated by NaNs. There is no
% distinction between outer contours and holes.

%  Written by: W. Stumpf
%  Copyright 1996-2003 The MathWorks, Inc.
%  $Revision: 1.5 $  $Date: 2003/05/20 21:00:38 $

	if (nargin < 2),	error('Incorrect number of input arguments.'),	end

	if (~isa(latcells,'cell') || ~isa(loncells,'cell'))
		error('Inputs must be cell arrays.')
	end

	if isempty(latcells);
		[lat2,lon2] = deal([],[]);
		return
	else
		[lat2,lon2] = deal(latcells{1},loncells{1});
	end

	for i=2:length(latcells)
		[lat2,lon2] = deal([lat2;NaN;latcells{i}], [lon2;NaN;loncells{i}]);
	end

% ---------------------------------------------------------------------
function [latcells,loncells] = polysplit(lat,lon)

%POLYSPLIT extracts segments of NaN-delimited polygon vectors to cell arrays
%
% [latcells,loncells] = POLYSPLIT(lat,lon) returns the NaN-delimited segments
% of the vectors lat and lon as cell arrays. Each element of the cell array
% contains one segment.

%  Written by:  W. Stumpf
%  Copyright 1996-2003 The MathWorks, Inc.
%  $Revision: 1.5 $  $Date: 2003/05/20 21:00:39 $

if (nargin < 2),	error('Incorrect number of input arguments.'), end

if ~isnumeric(lat) || ~isnumeric(lon)
	error('Inputs must be numeric arrays.')
end

%  Input dimension tests
if any([min(size(lat)) min(size(lon))] ~= 1) || any([ndims(lat) ndims(lon)] > 2)
	error('Latitude and longitude inputs must a vector')
elseif ~isequal(size(lat),size(lon))
	error('Inconsistent dimensions on lat and lon input')
end

%  Ensure at a terminating NaN in the vectors
if ~isnan( lat(numel(lat)) );    lat = [lat; NaN];   end
if ~isnan( lon(numel(lon)) );    lon = [lon; NaN];   end

%  Ensure vectors don't begin with NaNs
if isnan(lat(1)) || isnan(lon(1))
	lat = lat(2:numel(lat));
	lon = lon(2:numel(lon));
end

%  Find segment demarcations
indx  = find(isnan(lat));
indx2 = find(isnan(lon));

if ~isequal(indx,indx2)
	error('Lat and lon must have NaNs in same locations')
end

%  Extract each segment
latcells = cell(1,numel(indx));
loncells = cell(1,numel(indx));
ind = false(1, numel(indx));
for (i = 1:numel(indx))	% Pull segment out of main vectors
	if (i > 1)
		latcells{i} = lat(indx(i-1)+1:indx(i)-1);
		loncells{i} = lon(indx(i-1)+1:indx(i)-1);
		if (isempty(latcells{i})),	ind(i) = true;	end			% They may be empty if consecutive NaNs
	else
		latcells{i} = lat(1:indx(i)-1);
		loncells{i} = lon(1:indx(i)-1);
	end
end
if (any(ind))			% if we have empties, remove them
	latcells(ind) = [];
	loncells(ind) = [];
end

% ---------------------------------------------------------------------
function [xv,yv] = polycut(xp1,yp1)
%POLYCUT  Polygon branch cuts for holes.
%
%   [lat2,long2] = POLYCUT(lat,long) connects the contour and holes of polygons
%   using optimal branch cuts.  Polygons are input as NaN-delimited vectors, or 
%   as cell arrays containing individual polygons in each element with the outer 
%   face separated from the subsequent inner faces by NaNs. Multiple polygons 
%   outputs are separated by NaNs. 

%  Written by:  A. Kim
%  Copyright 1996-2003 The MathWorks, Inc.
%  $Revision: 1.6 $  $Date: 2003/05/20 21:00:37 $

if (nargin < 2),	error('Incorrect number of input arguments.'),	end

if isempty(xp1) && isempty(yp1)
	xv = [];	yv = [];
	return
end

% convert inputs to cell
if ~isa(xp1,'cell')
	[xp1,yp1]=deal({xp1},{yp1});
end

for i=1:length(xp1)
	[xp1{i},yp1{i}] = closefaces(xp1{i},yp1{i});		% close open polygons
end


% loop through each polygon input
xv = [];  yv = [];
for n=1:length(xp1)
   % extract polygon types and convert to structure
	[xc,yc,xh,yh,nanindx] = extractpoly(xp1(n),yp1(n));
	polystruc(1).x = xc;  polystruc(1).y = yc;
	for m=1:length(nanindx)-1
		indx = nanindx(m)+1:nanindx(m+1)-1;
		polystruc(m+1).x = xh(indx);
		polystruc(m+1).y = yh(indx);
	end
	icont = 1;  ihole = 2:length(nanindx);

	if isempty(ihole)

		x = [];  y = [];
		for k=1:length(polystruc)
			x = [x; nan; polystruc(k).x];  y = [y; nan; polystruc(k).y];
		end
		x(1) = [];  y(1) = [];
		if all(isnan(x)) && all(isnan(y))
			x = [];  y = [];
		end

	else

%***** polygon branch cut algorithm *****
		Np = length(polystruc);  pvec = (1:Np)';

% starting with the contour polygon, find the hole polygon closet to it,
% repeat for each succesive polygon
		polynum = icont;
		for k=1:Np-1
			pvec(pvec==polynum(k)) = [];
			x1 = polystruc(polynum(k)).x;  y1 = polystruc(polynum(k)).y;
			for m=1:length(pvec)
				x2 = polystruc(pvec(m)).x;  y2 = polystruc(pvec(m)).y;
				[i1(m),i2(m),dmin(m)] = mindist(x1,y1,x2,y2);
			end
			imin= find(dmin==min(dmin));  polynum(k+1,1) = pvec(imin(1));
			pt1(k+1,1) = i2(imin(1));  pt2(k,1) = i1(imin(1));
			clear i1 i2 dmin
		end
		pt1(1) = 1;  pt2(Np) = pt1(Np);

% record path of polygons along branch cuts
		x1 = [];  y1 = [];  x2 = [];  y2 = [];
		for k=1:Np
			poly = polynum(k);  i1 = pt1(k);  i2 = pt2(k);
			xp = polystruc(poly).x;  yp = polystruc(poly).y;
			if pt1(k)<pt2(k)
				x1 = [x1; xp(i1:i2)];			  y1 = [y1; yp(i1:i2)];
				x2 = [xp(i2:end); xp(2:i1); x2];  y2 = [yp(i2:end); yp(2:i1); y2];
			else
				x1 = [x1; xp(i1:end); xp(2:i2)];  y1 = [y1; yp(i1:end); yp(2:i2)];
				x2 = [xp(i2:i1); x2];			  y2 = [yp(i2:i1); y2];
			end
		end
		x = [x1; x2];  y = [y1; y2];

	end

	xv = [xv; x; nan];
	yv = [yv; y; nan];
	clear x y polystruc

end

xv(end) = [];  yv(end) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [i1,i2,dmin] = mindist(x1,y1,x2,y2)
%MINDIST  Minimum polygon distance points.
%   [I1,I2] = MINDIST(X1,Y1,X2,Y2) returns the indices of the closest
%   points between two polygons.
%
%   [I1,I2,DMIN] = MINDIST(X1,Y1,X2,Y2) also returns the minimum distances.
%  Written by:  A. Kim

% columnize inputs
x1 = x1(:);  y1 = y1(:);  x2 = x2(:);  y2 = y2(:);

% if self-enclosed, remove last point
if (x1(1) == x1(end) && y1(1) == y1(end))
	x1 = x1(1:end-1);	y1 = y1(1:end-1);
end
if (x2(1) == x2(end) && y2(1) == y2(end))
	x2 = x2(1:end-1);  y2 = y2(1:end-1);
end

% tile matrices for vectorized distance calculations
n1 = numel(x1);  n2 = numel(x2);
X1 = reshape(repmat(x1,1,n2)',1,n1*n2)';
Y1 = reshape(repmat(y1,1,n2)',1,n1*n2)';
I1 = reshape(repmat((1:n1)',1,n2)',1,n1*n2)';
X2 = repmat(x2,n1,1);
Y2 = repmat(y2,n1,1);
I2 = repmat((1:n2)',n1,1);

% calculate distances between points
dist = sqrt((X1-X2).^2 + (Y1-Y2).^2);
dmin = unique(min(dist));
imin = find(dist==dmin);
i1 = I1(imin(1));  i2 = I2(imin(1));

% ---------------------------------------------------------------------
function [lat,lon] = closefaces(latin,lonin)
%CLOSEFACES closes all faces of a polygon
%  $Revision: 1.5 $  $Date: 2003/07/17 15:25:52 $

[lat,lon] = deal(latin,lonin);

if ~iscell(lat),	[lat,lon] = polysplit(lat,lon);		end
   
for i=1:length(lat)
	[latface,lonface] = deal(lat{i},lon{i});
	[latfacecell,lonfacecell] = polysplit(latface,lonface);
	
	for (j=1:length(latfacecell))	% close open polygons
		if latfacecell{j}(1) ~= latfacecell{j}(end) | lonfacecell{j}(1) ~= lonfacecell{j}(end)
			latfacecell{j}(end+1) = latfacecell{j}(1);
			lonfacecell{j}(end+1) = lonfacecell{j}(1);
		end	
	end
	[lat{i},lon{i}]=polyjoin(latfacecell,lonfacecell);
end

if ~iscell(latin)
	[lat,lon] = polyjoin(lat,lon);
end

% ---------------------------------------------------------------------
function [xc,yc,xh,yh,nanindx] = extractpoly(xp,yp)
% EXTRACTPOLY  Extract contour and holes from polygon element.
% [XC,YC,XH,YH,NANINDX] = EXTRACTPOLY(XP,YP) extracts the contour and holes data from a polygon element.
%  Copyright 1996-2003 The MathWorks, Inc.
%  $Revision: 1.5 $  $Date: 2003/07/17 15:25:54 $

	% convert data to matrix
	if iscell(xp),		xp = xp{1};  yp = yp{1};	end
	
	% append a nan at the end
	x = [xp; nan];  y = [yp; nan];
	indx = find(isnan(x));
	
	% contour
	indxc = 1:indx(1)-1;
	xc = x(indxc);  yc = y(indxc);
	
	% holes
	xh = [];  yh = [];  nanindx = [];
	if length(indx)>1
		indxh = indx(1):indx(end);
		xh = x(indxh);  yh = y(indxh);
		nanindx = find(isnan(xh));
	end

% ---------------------------------------------------------------------
function [lat,lon,Z] = trimwrap(lat0,lon0,latlim,lonlim, Z0, opt)
%  [lat,lon] = TRIMWRAP(lat0,lon0,latlim,lonlim) trims a line map
%  to a region specified by latlim and lonlim.  Latlim and lonlim
%  are two element vectors, defining the latitude and longitude limits
%  respectively.  The inputs lat0 and lon0 must be vectors.
%  [lat,lon] = TRIMWRAP(lat0,lon0,latlim,lonlim, opt) wraps lines
%  at a 2PI discontinuity. NaN are inserted at each discontinuity point.
%  [lat,lon,Z] = TRIMWRAP(lat0,lon0,latlim,lonlim, Z0, opt) apply trim or wrap to Z0 too (can be empty)

if nargin < 4;   error('Incorrect number of arguments');   end

%  Test the inputs
latlim = sort(latlim);		lonlim = sort(lonlim);
if  ~isequal(size(latlim), size(lonlim), [1 2])
	error('Lat and lon limit inputs must be 2 element vectors')
end
if ~isequal(size(lat0),size(lon0))
     error('Inconsistent dimensions on input data')
end

trim = true;
if ( (nargin == 5 && ischar(Z0)) || (nargin == 6) )
	trim = false;
	if (nargin == 5)	Z0 = [];	end			% Input Z0 did its job, now we make it indicate that no Z in input
end
if (~isempty(Z0) && nargout == 2)				% Idiot choice. Z0 was transmitted, but Z not required on output.
	Z0 = [];
end
if (isempty(Z0) && nargout == 3)				% Here it's licit because Z0 doesn't need to be known in advance.
	Z = [];
end
%if (nargin == 5),	trim = false;		end		% Wrap arround the world

%  Ensure column vectors on input
if (size(lat0,2) > 1)
	lat0 = lat0(:);		lon0 = lon0(:);
end

%  Get the corners of the submap region
up    = latlim(2);		low  = latlim(1);
right = lonlim(2);		left = lonlim(1);

if (trim)
	% Determine the points which lie outside the region of interest
	indx = find(lat0 < low | lat0 > up | lon0 < left | lon0 > right);
else		% Wrap
	% Determine the points which fall at the end of the world
	dife = diff(lon0);
	indx = (dife > 180 | dife < -180);
	if (~any(indx)),
		indx = [];
	else
		indx = [false; indx];
	end
end

%  Extract the found points by first replacing the points outside the map
%  with NaNs.  Then eliminate multiple NaNs in the vector.  This is
%  necessary incase a line segement exits and enters the trim box, so
%  as to NaN clip the exit/enter point.

if (~isempty(indx) && trim)				% Trim
	lat = lat0;             lon = lon0;
	if (~isempty(Z0))		Z = Z0;		end
	lat(indx) = NaN;        lon(indx) = NaN;	
	if (~isempty(Z0))		Z(indx) = NaN;		end
	nanloc = isnan(lat);	[r,c] = size(nanloc);
	nanloc = find(nanloc(1:r-1,:) & nanloc(2:r,:));
	lat(nanloc) = [];  lon(nanloc) = [];
	if (~isempty(Z0))		Z(nanloc) = [];		end
elseif (~isempty(indx) && ~trim)		% Wrap
	nanloc = find(indx);
	% Create vectors with the total number of points
	lon = zeros(numel(lon0)+numel(nanloc),1);		lat = zeros(numel(lon0)+numel(nanloc),1);
	if (~isempty(Z0))		Z = lat;	end
	nanloc = [1; nanloc];

	for (k = 1:numel(nanloc)-1)
		off = (k-1);
		lon( nanloc(k)+off:nanloc(k+1)-1+off ) = lon0( nanloc(k):nanloc(k+1)-1 );		lon(nanloc(k+1)+off) = NaN;
		lat( nanloc(k)+off:nanloc(k+1)-1+off ) = lat0( nanloc(k):nanloc(k+1)-1 );		lat(nanloc(k+1)+off) = NaN;
		if (~isempty(Z0))
			Z( nanloc(k)+off:nanloc(k+1)-1+off ) = Z0( nanloc(k):nanloc(k+1)-1 );		Z(nanloc(k+1)+off) = NaN;
		end
	end
	% Copy the oints since last NaN till end of arrays
	lon( nanloc(end)+1+off:numel(lon) ) = lon0( nanloc(end):numel(lon0) );
	lat( nanloc(end)+1+off:numel(lon) ) = lat0( nanloc(end):numel(lat0) );
	if (~isempty(Z0))
		Z( nanloc(end)+1+off:numel(lon) ) = Z0( nanloc(end):numel(Z0) );
	end
else
	lat = lat0;             lon = lon0;
	if (~isempty(Z0))		Z = Z0;		end
end

% ---------------------------------------------------------------------
function [lat,lon] = maptrimp(lat0,lon0,latlim,lonlim)
%  [lat,lon] = MAPTRIMP(lat0,lon0,latlim,lonlim) trims a patch map
%  to a region specified by latlim and lonlim.  Latlim and lonlim
%  are two element vectors, defining the latitude and longitude limits
%  respectively.  The inputs lat0 and lon0 must be vectors representing
%  patch map vector data.

%  Copyright 1996-2003 The MathWorks, Inc.
%  Written by:  E. Byrns, E. Brown


if nargin < 4;   error('Incorrect number of arguments');   end

%  Test the inputs
if  ~isequal(sort(size(latlim)),sort(size(lonlim)),[1 2])
	error('Lat and lon limit inputs must be 2 element vectors')
end

%  Get the corners of the submap region
up    = max(latlim);   low  = min(latlim);
right = max(lonlim);   left = min(lonlim);

%  Copy the input data and ensure column vectors.
lat = lat0(:);   lon = lon0(:);

%  Get the vector of patch items and remove any NaN padding
%  at the beginning or end of the column.  This eliminates potential
%  multiple NaNs at the beginning and end of the patches.
while (isnan(lat(1)) || isnan(lon(1)))
	lat(1) = [];   lon(1) = [];
end
while ( isnan(lat(length(lat))) || isnan(lon(length(lon))) )   
	lat(length(lat)) = [];   lon(length(lon)) = [];
end

%  Add a NaN to the end of the data vector.  Necessary for processing of multiple patches.
lat(length(lat)+1) = NaN;   lon(length(lon)+1) = NaN;

%  Find the individual patches and then trim the data
indx = find(isnan(lon) | isnan(lat));
if isempty(indx);   indx = length(lon)+1;   end

for i = 1:length(indx)

	if (i == 1),	startloc = 1;
	else			startloc = indx(i-1)+1;
	end
	endloc   = indx(i)-1;
	
	indices = (startloc:endloc)';   %  Indices will be empty if NaNs are
                                    %  neighboring in the vector data.
	if ~isempty(indices)            %  Should not happen, but test just in case

		%  Patches which lie completely outside the trim window.  Replace
		%  with NaNs and then eliminate it entirely later.  Replacing with
		%  NaNs is useful so that the indexing with indices is not messed
		%  up if an entire patch is eliminated at this point.
		
		%  If at least one point of the patch edge does not lie with the
		%  specified window limits, then the entire patch is trimmed.
	
         if ~any(lon(indices) >= left & lon(indices) <= right & ...
			 	lat(indices) >= low  & lat(indices) <= up)
				lon(indices) = NaN;    lat(indices) = NaN;
         end
	
		%  Need to only test along edge since patch must lie somehow within the window.
	
		%  Points off the bottom
		loctn = find( lon(indices) < left );
		if ~isempty(loctn);  lon(indices(loctn)) = left;     end
	
		%  Points off the top
		loctn = find( lon(indices) > right );
		if ~isempty(loctn);   lon(indices(loctn)) = right;   end
	
		%  Points off the left
		loctn = find( lat(indices) < low );
		if ~isempty(loctn);   lat(indices(loctn)) = low;     end
	
		%  Points off the right
		loctn = find( lat(indices) > up );
		if ~isempty(loctn);   lat(indices(loctn)) = up;      end
	end
end


%  Eliminate multiple NaNs in the vector.  Will occur if a patch
%  lies entirely outside the window of interest.
if ~isempty(lat)
	nanloc = isnan(lat);	[r,c] = size(nanloc);
	nanloc = find(nanloc(1:r-1,:) & nanloc(2:r,:));
	lat(nanloc) = [];  lon(nanloc) = [];
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
function [yvec,xvec,trimpts] = trimpatch(yvec,ylim,xvec,xlim)
%TRIMPATCH  Trim patch vector data exceeding projection limits
%
%  Purpose
%  TRIMPATCH will identify points in patch vector data
%  which exceed projection limits.  The projection limits
%  are defined by the lower and upper inputs.  Points outside
%  the projection range are replaced with their respective limits,
%  thereby trimming them from the display.  If a patch lies completely
%  outside the trim limits, it is completely replaced with NaNs.
%
%  Synopsis
%           [ymat,xmat,trimpts] = trimpatch(ymat,ylim,xmat,xlim)
%
%  Copyright 1996-2003 The MathWorks, Inc.

if nargin ~= 4;  error('Incorrect number of arguments');  end

%  Dimension tests
if ndims(xvec) > 2 || ndims(yvec) > 2
    error('Pages not allowed for xmat or ymat')
elseif ~isequal(size(xvec),size(yvec))
    error('Inconsistent xmat and ymat dimensions')
end

xvec = xvec(:);
yvec = yvec(:);

%  Find the individual patches
indx = find(isnan(xvec) | isnan(yvec));
if isempty(indx);   indx = length(xvec)+1;   end

trimpts = [];

for i = 1:length(indx)
    if (i == 1),    startloc = 1;
    else            startloc = indx(i-1)+1;
    end
    endloc   = indx(i)-1;

    indices = (startloc:endloc)';   %  Indices will be empty if NaNs are
	                                %  neighboring in the vector data.
    if ~isempty(indices)            %  Should not happen, but test just in case

        %  Patches which lie completely outside the trim window.
        %  If at least one point of the patch edge does not lie with the
        %  specified window limits, then the entire patch is trimmed.
	
        if ~any(xvec(indices) >= min(xlim) & xvec(indices) <= max(xlim) & ...
             yvec(indices) >= min(ylim)  & yvec(indices) <= max(ylim))
            trimpts = [trimpts; indices ones(size(indices)) yvec(indices) xvec(indices)];
		    xvec(indices) = NaN;    yvec(indices) = NaN;
        end
	
        %  Need to only test along edge since patch must lie somehow within
        %  the window.  Make sure that the original data is saved before
        %  the points are replaced with the limit data.
	
        %  Points off the bottom
        loctn = find( xvec(indices) < min(xlim) );
        if ~isempty(loctn)
            trimpts = [trimpts; indices(loctn) ones(size(loctn)) yvec(indices(loctn)) xvec(indices(loctn))];
			xvec(indices(loctn)) = min(xlim);
        end
	
        %  Points off the top
        loctn = find( xvec(indices) > max(xlim) );
        if ~isempty(loctn)
            trimpts = [trimpts; indices(loctn) ones(size(loctn)) yvec(indices(loctn)) xvec(indices(loctn))];
			xvec(indices(loctn)) = max(xlim);
        end
	
        %  Points off the left
        loctn = find( yvec(indices) < min(ylim) );
        if ~isempty(loctn)
            trimpts = [trimpts; indices(loctn) ones(size(loctn)) yvec(indices(loctn)) xvec(indices(loctn))];
            yvec(indices(loctn)) = min(ylim);
        end
	
        %  Points off the right
        loctn = find( yvec(indices) > max(ylim) );
        if ~isempty(loctn)
            trimpts = [trimpts; indices(loctn) ones(size(loctn)) yvec(indices(loctn)) xvec(indices(loctn))];
            yvec(indices(loctn)) = max(ylim);
        end
    end
end

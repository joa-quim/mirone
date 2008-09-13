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
if any([min(size(lat)) min(size(lon))] ~= 1) | any([ndims(lat) ndims(lon)] > 2)
	error('Latitude and longitude inputs must a vector')
elseif ~isequal(size(lat),size(lon))
	error('Inconsistent dimensions on lat and lon input')
end

%  Ensure at a terminating NaN in the vectors
if ~isnan( lat(length(lat)) );    lat = [lat; NaN];   end
if ~isnan( lon(length(lon)) );    lon = [lon; NaN];   end

%  Ensure vectors don't begin with NaNs
if isnan(lat(1)) || isnan(lon(1))
	lat = lat(2:length(lat));
	lon = lon(2:length(lon));
end

%  Find segment demarcations
indx  = find(isnan(lat));
indx2 = find(isnan(lon));

if ~isequal(indx,indx2)
	error('Lat and lon must have NaNs in same locations')
end

%  Extract each segment
for (i = 1:length(indx))	% Pull segment out of main vectors
	if (i > 1)
		latcells{i} = lat(indx(i-1)+1:indx(i)-1);
		loncells{i} = lon(indx(i-1)+1:indx(i)-1);
	else
		latcells{i} = lat(1:indx(i)-1);
		loncells{i} = lon(1:indx(i)-1);
	end
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

if isempty(xp1) & isempty(yp1)
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

	if length(ihole)==0

		x = [];  y = [];
		for k=1:length(polystruc)
			x = [x; nan; polystruc(k).x];  y = [y; nan; polystruc(k).y];
		end
		x(1) = [];  y(1) = [];
		if all(isnan(x)) & all(isnan(y))
			x = [];  y = [];
		end

	else

%***** polygon branch cut algorithm *****
		Np = length(polystruc);  pvec = (1:Np)';

% starting with the contour polygon, find the hole polygon closet to it,
% repeat for each succesive polygon
		polynum = icont;
		for k=1:Np-1
			pvec(find(pvec==polynum(k))) = [];
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

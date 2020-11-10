function [xi,yi,zi] = griddata_j(x,y,z,xi,yi,method,options)
%GRIDDATA Data gridding and surface fitting.
%   ZI = GRIDDATA(X,Y,Z,XI,YI) fits a surface of the form Z = F(X,Y) to the
%   data in the (usually) nonuniformly-spaced vectors (X,Y,Z). GRIDDATA
%   interpolates this surface at the points specified by (XI,YI) to produce
%   ZI.  The surface always goes through the data points. XI and YI are
%   usually a uniform grid (as produced by MESHGRID) and is where GRIDDATA
%   gets its name.
%
%   XI can be a row vector, in which case it specifies a matrix with
%   constant columns. Similarly, YI can be a column vector and it specifies
%   a matrix with constant rows.
%
%   [XI,YI,ZI] = GRIDDATA(X,Y,Z,XI,YI) also returns the XI and YI formed
%   this way (the results of [XI,YI] = MESHGRID(XI,YI)).
%
%   [...] = GRIDDATA(X,Y,Z,XI,YI,METHOD) where METHOD is one of
%       'linear'    - Triangle-based linear interpolation (default)
%       'cubic'     - Triangle-based cubic interpolation
%       'nearest'   - Nearest neighbor interpolation
%   defines the type of surface fit to the data. The 'cubic' and 'v4'
%   methods produce smooth surfaces while 'linear' and 'nearest' have
%   discontinuities in the first and zero-th derivative respectively.  All
%   the methods except 'v4' are based on a Delaunay triangulation of the
%   data.
%   If METHOD is [], then the default 'linear' method will be used.
%
%   [...] = GRIDDATA(X,Y,Z,XI,YI,METHOD,OPTIONS) specifies a cell array 
%   of strings OPTIONS to be used as options in Qhull via DELAUNAYN. 
%   If OPTIONS is [], the default DELAUNAYN options will be used.
%   If OPTIONS is {''}, no options will be used, not even the default.
%
%   See also GRIDDATA3, GRIDDATAN, DELAUNAY, INTERP2, MESHGRID, DELAUNAYN.

%   Copyright 1984-2003 The MathWorks, Inc. 
%   $Revision: 5.33.4.3 $  $Date: 2004/01/16 20:04:55 $

error(nargchk(5,7,nargin))

[msg,x,y,z,xi,yi] = xyzchk(x,y,z,xi,yi);
if ~isempty(msg), error(msg); end

if ( nargin < 6 || isempty(method) ),  method = 'linear'; end
if ~ischar(method), 
  error('MATLAB:griddata:InvalidMethod',...
        'METHOD must be one of ''linear'',''cubic'',''nearest''.');
end

if (nargin == 7)
    if ~iscellstr(options)
        error('MATLAB:OptsNotStringCell','OPTIONS should be cell array of strings.');           
    end
    opt = options;
else
    opt = [];
end

% Sort x and y so duplicate points can be averaged before passing to delaunay

%Need x,y and z to be column vectors
sz = numel(x);
x = reshape(x,sz,1);
y = reshape(y,sz,1);
z = reshape(z,sz,1);
sxyz = sortrows([x y z],[2 1]);
x = sxyz(:,1);
y = sxyz(:,2);
z = sxyz(:,3);
myeps = max(max(abs(x)),max(abs(y)))*eps^(1/3);
ind = [0; ((abs(diff(y)) < myeps) & (abs(diff(x)) < myeps)); 0];

if (sum(ind) > 0)
  warning('MATLAB:griddata:DuplicateDataPoints',['Duplicate x-y data points ' ...
            'detected: using average of the z values.']);
  fs = find(ind(1:end-1) == 0 & ind(2:end) == 1);
  fe = find(ind(1:end-1) == 1 & ind(2:end) == 0);
  for i = 1 : length(fs)
    % averaging z values
    z(fe(i)) = mean(z(fs(i):fe(i)));
  end
  x = x(~ind(2:end));
  y = y(~ind(2:end));
  z = z(~ind(2:end));
end

switch lower(method),
  case 'linear'
    zi = linear(x,y,z,xi,yi,opt);
  case 'cubic'
    zi = cubic(x,y,z,xi,yi,opt);
  case 'nearest'
    zi = nearest(x,y,z,xi,yi,opt);
  otherwise
    error('MATLAB:griddata:UnknownMethod', 'Unknown method.');
end
  
if nargout<=1, xi = zi; end


%------------------------------------------------------------
function zi = linear(x,y,z,xi,yi,opt)
%LINEAR Triangle-based linear interpolation

%   Reference: David F. Watson, "Contouring: A guide
%   to the analysis and display of spacial data", Pergamon, 1994.

siz = size(xi);
xi = xi(:); yi = yi(:); % Treat these as columns
x = x(:); y = y(:); % Treat these as columns

% Triangularize the data
if isempty(opt)
    tri = delaunayn([x y]);
else
    tri = delaunayn([x y],opt);
end
    
if isempty(tri),
  warning('MATLAB:griddata:CannotTriangulate','Data cannot be triangulated.');
  zi = repmat(NaN,size(xi));
  return
end

% Find the nearest triangle (t)
t = tsearch(x,y,tri,xi,yi);

% Only keep the relevant triangles.
out = find(isnan(t));
if ~isempty(out), t(out) = ones(size(out)); end
tri = tri(t,:);

% Compute Barycentric coordinates (w).  P. 78 in Watson.
del = (x(tri(:,2))-x(tri(:,1))) .* (y(tri(:,3))-y(tri(:,1))) - ...
      (x(tri(:,3))-x(tri(:,1))) .* (y(tri(:,2))-y(tri(:,1)));
w(:,3) = ((x(tri(:,1))-xi).*(y(tri(:,2))-yi) - ...
          (x(tri(:,2))-xi).*(y(tri(:,1))-yi)) ./ del;
w(:,2) = ((x(tri(:,3))-xi).*(y(tri(:,1))-yi) - ...
          (x(tri(:,1))-xi).*(y(tri(:,3))-yi)) ./ del;
w(:,1) = ((x(tri(:,2))-xi).*(y(tri(:,3))-yi) - ...
          (x(tri(:,3))-xi).*(y(tri(:,2))-yi)) ./ del;
w(out,:) = zeros(length(out),3);

z = z(:).'; % Treat z as a row so that code below involving
            % z(tri) works even when tri is 1-by-3.
zi = sum(z(tri) .* w,2);

zi = reshape(zi,siz);

if ~isempty(out), zi(out) = NaN; end
%------------------------------------------------------------

%------------------------------------------------------------
function zi = cubic(x,y,z,xi,yi,opt)
%TRIANGLE Triangle-based cubic interpolation

%   Reference: T. Y. Yang, "Finite Element Structural Analysis",
%   Prentice Hall, 1986.  pp. 446-449.
%
%   Reference: David F. Watson, "Contouring: A guide
%   to the analysis and display of spacial data", Pergamon, 1994.

% Triangularize the data
if isempty(opt)
    tri = delaunayn([x(:) y(:)]);
else
    tri = delaunayn([x(:) y(:)],opt);
end
if isempty(tri), 
  warning('MATLAB:griddata:CannotTriangulate','Data cannot be triangulated.');
  zi = repmat(NaN,size(xi));
  return
end

% Find the nearest triangle (t)
t = tsearch(x,y,tri,xi,yi);

zi = cubicmx(x,y,z,xi,yi,tri,t);
%------------------------------------------------------------

%------------------------------------------------------------
function zi = nearest(x,y,z,xi,yi,opt)
%NEAREST Triangle-based nearest neightbor interpolation

%   Reference: David F. Watson, "Contouring: A guide
%   to the analysis and display of spacial data", Pergamon, 1994.

siz = size(xi);
xi = xi(:); yi = yi(:); % Treat these a columns
x = x(:); y = y(:); z = z(:); % Treat these as columns

% Triangularize the data
if isempty(opt)
    tri = delaunayn([x y]);
else
    tri = delaunayn([x y],opt);
end
if isempty(tri), 
  warning('MATLAB:griddata:CannotTriangulate','Data cannot be triangulated.');
  zi = repmat(NaN,size(xi));
  return
end

% Find the nearest vertex
k = dsearch(x,y,tri,xi,yi);
zi = k;
d = find(isfinite(k));
zi(d) = z(k(d));
zi = reshape(zi,siz);

%----------------------------------------------------------
function t = delaunayn(x,options)
%DELAUNAYN  N-D Delaunay tessellation.
%   T = DELAUNAYN(X) returns a set of simplices such that no data points of
%   X are contained in any circumspheres of the simplices. The set of
%   simplices forms the Delaunay tessellation. X is an m-by-n array
%   representing m points in n-D space. T is a numt-by-(n+1) array where
%   each row is the indices into X of the vertices of the corresponding
%   simplex. When the simplices cannot be computed (such as when X is
%   degenerate, or X is empty), an empty matrix is returned.
%
%   DELAUNAYN uses Qhull. 
%
%   T = DELAUNAYN(X,OPTIONS) specifies a cell array of strings OPTIONS to
%   be used as options in Qhull. The default options are:
%                                 {'Qt','Qbb','Qc'} for 2D and 3D input,
%                                 {'Qt','Qbb','Qc','Qx'} for 4D and higher.  
%   If OPTIONS is [], the default options will be used.
%   If OPTIONS is {''}, no options will be used, not even the default.
%   For more information on Qhull options, see http://www.qhull.org.
%
%   Example:
%      X = [-0.5 -0.5  -0.5;
%           -0.5 -0.5   0.5;
%           -0.5  0.5  -0.5;
%           -0.5  0.5   0.5;
%            0.5 -0.5  -0.5;
%            0.5 -0.5   0.5;
%            0.5  0.5  -0.5;
%            0.5  0.5   0.5];
%      T = delaunayn(X);
%   errors, but hints that adding 'Qz' to the default options might help.
%      T = delaunayn(X,{'Qt','Qbb','Qc','Qz'});
%   To visualize this answer you can use the TETRAMESH function:
%      tetramesh(T,X)
%
%   See also QHULL, VORONOIN, CONVHULLN, DELAUNAY, DELAUNAY3, TETRAMESH.

%   Copyright 1984-2003 The MathWorks, Inc.
%   $Revision: 1.20.4.5 $ $Date: 2004/01/16 20:04:53 $

if nargin < 1
    error('MATLAB:delaunayn:NotEnoughInputs', 'Needs at least 1 input.');
end
if isempty(x), t = []; return; end
[x,idx,jdx] = unique(x,'rows');
idx = idx';

[m,n] = size(x);

if m < n+1,
  error('MATLAB:delaunayn:NotEnoughPtsForTessel',...
        'Not enough unique points to do tessellation.');
end
if any(isinf(x(:)) | isnan(x(:)))
  error('MATLAB:delaunayn:CannotTessellateInfOrNaN',...
        'Data containing Inf or NaN cannot be tessellated.');
end
if (m == n+1)
    t = idx(1:n+1);
  return;
end

%default options
if (n >= 4),   opt = 'Qt Qbb Qc Qx';
else            opt = 'Qt Qbb Qc';  end

if ( nargin > 1 && ~isempty(options) )
    if ~iscellstr(options)
        error('MATLAB:delaunayn:OptsNotStringCell',...
              'OPTIONS should be cell array of strings.');
    end
    sp = {' '};
    c = strcat(options,sp);
    opt = cat(2,c{:});   
end

t = qhullmx(x', 'd ', opt);

% try to get rid of zero volume simplices. They are generated
% because of the fuzzy jiggling.

[mt, nt] = size(t);
v = true(mt,1);

seps = eps^(4/5)*max(abs(x(:)));
try
    for i=1:mt
        val = abs(det(x(t(i,1:nt-1),:)-x(t(i,nt)*ones([nt-1, 1]),:))); 
        if val < seps
           v(i) = false;
        end
    end
catch
    error('MATLAB:delaunayn:InvalidOpts',...
          'Bad options choice. Try using different options.');
end

t = t(v,:);
t = idx(t);
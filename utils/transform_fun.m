%#mex
function  varargout = transform_fun(opt,varargin)

switch opt
    case 'cp2tform'
        varargout{1} = cp2tform(varargin{:});
    case 'imrotate'
        varargout{1} = imrotate(varargin{:});
    case 'imtransform'
        [B,XDATA,YDATA] = imtransform(varargin{:});
        varargout{1} = B;   varargout{2} = XDATA;   varargout{3} = YDATA;
    case 'maketform'
        varargout{1} = maketform(varargin{:});
    case 'makeresampler'
        varargout{1} = makeresampler(varargin{:});
    case 'tformarray'
        varargout{1} = tformarray(varargin{:});
    case 'tforminv'
        varargout = tform('inv', nargout, varargin{:});
    case 'tformfwd'
        varargout = tform('fwd', nargout, varargin{:});
end

function [trans,uv,xy,uv_dev,xy_dev] = cp2tform(varargin)
%CP2TFORM Infer spatial transformation from control point pairs.
%   CP2TFORM takes pairs of control points and uses them to infer a
%   spatial transformation. 
%
%   TFORM = CP2TFORM(INPUT_POINTS,BASE_POINTS,TRANSFORMTYPE) returns a TFORM
%   structure containing a spatial transformation. INPUT_POINTS is an M-by-2
%   double matrix containing the X and Y coordinates of control points in
%   the image you want to transform. BASE_POINTS is an M-by-2 double matrix
%   containing the X and Y coordinates of control points in the base
%   image. TRANSFORMTYPE can be 'linear conformal', 'affine', 'projective',
%   'polynomial', 'piecewise linear' or 'lwm'. See the reference page for
%   CP2TFORM for information about choosing TRANSFORMTYPE.
%
%   TFORM = CP2TFORM(CPSTRUCT,TRANSFORMTYPE) works on a CPSTRUCT structure
%   that contains the control point matrices for the input and base
%   images. The Control Point Selection Tool, CPSELECT, creates the
%   CPSTRUCT.
%
%   [TFORM,INPUT_POINTS,BASE_POINTS] = CP2TFORM(CPSTRUCT,...) returns the
%   control points that were actually used in INPUT_POINTS, and
%   BASE_POINTS. Unmatched and predicted points are not used. See
%   CPSTRUCT2PAIRS.
%
%   TFORM = CP2TFORM(INPUT_POINTS,BASE_POINTS,'polynomial',ORDER) 
%   ORDER specifies the order of polynomials to use. ORDER can be 2, 3, 4.
%   If you omit ORDER, it defaults to 3.
%
%   TFORM = CP2TFORM(CPSTRUCT,'polynomial',ORDER) works on a CPSTRUCT
%   structure.
%
%   TFORM = CP2TFORM(INPUT_POINTS,BASE_POINTS,'piecewise linear') Creates a
%   delaunay triangulation of the base control points, and maps
%   corresponding input control points to the base control points. The
%   mapping is linear (affine) for each triangle, and continuous across the
%   control points, but not continuously differentiable as each triangle has
%   its own mapping.
%
%   TFORM = CP2TFORM(CPSTRUCT,'piecewise linear') works on a CPSTRUCT
%   structure.
%
%   TFORM = CP2TFORM(INPUT_POINTS,BASE_POINTS,'lwm',N) The local weighted
%   mean (lwm) method creates a mapping, by inferring a polynomial at each
%   control point using neighboring control points. The mapping at any
%   location depends on a weighted average of these polynomials.  You can
%   optionally specify the number of points, N, used to infer each
%   polynomial. The N closest points are used to infer a polynomial of order
%   2 for each control point pair. If you omit N, it defaults to 12. N can
%   be as small as 6, BUT making N small risks generating ill-conditioned
%   polynomials.
%
%   TFORM = CP2TFORM(CPSTRUCT,'lwm',N) works on a CPSTRUCT structure.
%
%   [TFORM,INPUT_POINTS,BASE_POINTS,INPUT_POINTS_BAD,BASE_POINTS_BAD] = ...
%        CP2TFORM(INPUT_POINTS,BASE_POINTS,'piecewise linear') 
%   returns the control points that were actually used in INPUT_POINTS and
%   BASE_POINTS, and the control points that were eliminated because they
%   were middle vertices of degenerate fold-over triangles in
%   INPUT_POINTS_BAD and BASE_POINTS_BAD.
%
%   [TFORM,INPUT_POINTS,BASE_POINTS,INPUT_POINTS_BAD,BASE_POINTS_BAD] = ...
%        CP2TFORM(CPSTRUCT,'piecewise linear') works on a CPSTRUCT structure.
%
%   TRANSFORMTYPE
%   -------------
%   CP2TFORM requires a minimum number of control point pairs to infer a
%   TFORM structure of each TRANSFORMTYPE:
%
%       TRANSFORMTYPE         MINIMUM NUMBER OF PAIRS
%       -------------         -----------------------
%       'linear conformal'               2 
%       'affine'                         3 
%       'projective'                     4 
%       'polynomial' (ORDER=2)           6
%       'polynomial' (ORDER=3)          10
%       'polynomial' (ORDER=4)          15
%       'piecewise linear'               4
%       'lwm'                            6
%      
%   When TRANSFORMTYPE is 'linear conformal', 'affine', 'projective', or
%   'polynomial', and INPUT_POINTS and BASE_POINTS (or CPSTRUCT) have the
%   minimum number of control points needed for a particular transformation,
%   the coefficients are found exactly. If INPUT_POINTS and BASE_POINTS have
%   more than the minimum, a least squares solution is found. See MLDIVIDE.
%
%   Note
%   ----
%   When either INPUT_POINTS or BASE_POINTS has a large offset with
%   respect to ther origin (relative to range of values that it spans), the
%   points are shifted to center their bounding box on the origin before
%   fitting a TFORM structure.  This enhances numerical stability and
%   is handled transparently by wrapping the origin-centered TFORM within a
%   custom TFORM that automatically applies and undoes the coordinate shift
%   as needed. This means that fields(T) may give different results for
%   different coordinate inputs, even for the same transformtype.
%
%   Example
%   -------
%   I = checkerboard;
%   J = imrotate(I,30);
%   base_points = [11 11; 41 71];
%   input_points = [14 44; 70 81];
%   cpselect(J,I,input_points,base_points);
%
%   t = cp2tform(input_points,base_points,'linear conformal');
%
%   % Recover angle and scale by checking how a unit vector 
%   % parallel to the x-axis is rotated and stretched. 
%   u = [0 1]; 
%   v = [0 0]; 
%   [x, y] = tformfwd(t, u, v); 
%   dx = x(2) - x(1); 
%   dy = y(2) - y(1); 
%   angle = (180/pi) * atan2(dy, dx) 
%   scale = 1 / sqrt(dx^2 + dy^2) 

%   Copyright 1993-2005 The MathWorks, Inc. 
%   $Revision: 1.10.4.7 $  $Date: 2005/11/15 00:58:13 $

[uv, xy, method, options] = cp2tform_ParseInputs(varargin{:});

% initialize deviation matrices
xy_dev = [];    uv_dev = [];

% Assign function according to method and 
% set K = number of control point pairs needed. 
switch method
  case 'affine',            findT_fcn = @findAffineTransform;
  case 'projective',        findT_fcn = @findProjectiveTransform;
  case 'linear conformal',  findT_fcn = @findLinearConformalTransform;
  case 'polynomial',        findT_fcn = @findPolynomialTransform;
  case 'piecewise linear',  findT_fcn = @findPiecewiseLinear;
  case 'lwm',               findT_fcn = @findLWM;
end    

% error if user enters too few control point pairs
M = size(uv,1);
if (M < options.K)
    [msg, eid] = CountError(options.K,method);
    error(eid,msg);
end

% get offsets to apply to before/after spatial transformation
uvShift = getShift(uv);
xyShift = getShift(xy);
needToShift = any([uvShift xyShift] ~= 0);

if ~needToShift     % infer transform
    [trans, output] = feval(findT_fcn,uv,xy,options);
else                % infer transform for shifted data
    [tshifted, output] = feval(findT_fcn,applyShift(uv,uvShift),applyShift(xy,xyShift),options);

    % construct custom tform with tshifted between forward and inverse shifts
    tdata = struct('uvShift',uvShift,'xyShift',xyShift,'tshifted',tshifted);
    trans = maketform('custom',2,2,@fwd,@inverse,tdata);
end

if strcmp(method,'piecewise linear')
    uv = undoShift(output.uv,uvShift);  % No-ops if needToShift
    xy = undoShift(output.xy,xyShift);  % is false.
    uv_dev = output.uv_dev;
    xy_dev = output.xy_dev;
end

%-------------------------------
function shift = getShift(points)
tol = 1e+3;
minPoints = min(points);
maxPoints = max(points);
center = (minPoints + maxPoints) / 2;
span = maxPoints - minPoints;
if (span(1) > 0 && abs(center(1))/span(1) > tol) ||...
   (span(2) > 0 && abs(center(2))/span(2) > tol)
    shift = center;
else
    shift = [0 0];
end

%-------------------------------
function shiftedPoints = applyShift(points,shift)
shiftedPoints = points - repmat(shift,[size(points,1) 1]);

%-------------------------------
function points = undoShift(shiftedPoints,shift)
points = shiftedPoints + repmat(shift,[size(shiftedPoints,1) 1]);

%-------------------------------
function x = fwd(u,t)
x = undoShift(tformfwd(applyShift(u,t.tdata.uvShift),t.tdata.tshifted),t.tdata.xyShift);

%-------------------------------
function u = inverse(x,t)
u = undoShift(tforminv(applyShift(x,t.tdata.xyShift),t.tdata.tshifted), t.tdata.uvShift);

%--------------------- Function  findAffineTransform -----------------------------------
function [trans, output] = findAffineTransform(uv,xy,options)
%
% For an affine transformation:
%
%
%                     [ A D 0 ]
% [u v 1] = [x y 1] * [ B E 0 ]
%                     [ C F 1 ]
%
% There are 6 unknowns: A,B,C,D,E,F
%
% Another way to write this is:
%
%                   [ A D ]
% [u v] = [x y 1] * [ B E ]
%                   [ C F ]
%
% Rewriting the above matrix equation:
% U = X * T, where T = reshape([A B C D E F],3,2)
%
% With 3 or more correspondence points we can solve for T,
% T = X\U which gives us the first 2 columns of T, and
% we know the third column must be [0 0 1]'.

K = options.K;
M = size(xy,1);
X = [xy ones(M,1)];

U = uv;         % just solve for the first two columns of T

% We know that X * T = U
if rank(X) >= K
    Tinv = X \ U;
else
    [msg, eid] = RankError(K,'affine');    error(eid,msg);
end

Tinv(:,3) = [0 0 1]';       % add third column

T = inv(Tinv);
T(:,3) = [0 0 1]';

trans = maketform('affine', T);
output = [];

%--------------- Function  findProjectiveTransform ----------------------------
function [trans, output] = findProjectiveTransform(uv,xy,options)
%
% For a projective transformation:
%
% u = (Ax + By + C)/(Gx + Hy + I)
% v = (Dx + Ey + F)/(Gx + Hy + I)
%
% Assume I = 1, multiply both equations, by denominator:
%
% u = [x y 1 0 0 0 -ux -uy] * [A B C D E F G H]'
% v = [0 0 0 x y 1 -vx -vy] * [A B C D E F G H]'
%
% With 4 or more correspondence points we can combine the u equations and
% the v equations for one linear system to solve for [A B C D E F G H]:
%
% [ u1  ] = [ x1  y1  1  0   0   0  -u1*x1  -u1*y1 ] * [A]
% [ u2  ] = [ x2  y2  1  0   0   0  -u2*x2  -u2*y2 ]   [B]
% [ u3  ] = [ x3  y3  1  0   0   0  -u3*x3  -u3*y3 ]   [C]
% [ u1  ] = [ x4  y4  1  0   0   0  -u4*x4  -u4*y4 ]   [D]
% [ ... ]   [ ...                                  ]   [E]
% [ un  ] = [ xn  yn  1  0   0   0  -un*xn  -un*yn ]   [F]
% [ v1  ] = [ 0   0   0  x1  y1  1  -v1*x1  -v1*y1 ]   [G]
% [ v2  ] = [ 0   0   0  x2  y2  1  -v2*x2  -v2*y2 ]   [H]
% [ v3  ] = [ 0   0   0  x3  y3  1  -v3*x3  -v3*y3 ]
% [ v4  ] = [ 0   0   0  x4  y4  1  -v4*x4  -v4*y4 ]
% [ ... ]   [ ...                                  ]  
% [ vn  ] = [ 0   0   0  xn  yn  1  -vn*xn  -vn*yn ]
%
% Or rewriting the above matrix equation:
% U = X * Tvec, where Tvec = [A B C D E F G H]'
% so Tvec = X\U.

K = options.K;          M = size(xy,1);
x = xy(:,1);            y = xy(:,2);
vec_1 = ones(M,1);      vec_0 = zeros(M,1);
u = uv(:,1);            v = uv(:,2);

U = [u; v];
X = [x      y      vec_1  vec_0  vec_0  vec_0  -u.*x  -u.*y;
     vec_0  vec_0  vec_0  x      y      vec_1  -v.*x  -v.*y  ];

% We know that X * Tvec = U
if rank(X) >= 2*K 
    Tvec = X \ U;    
else
    [msg, eid] = RankError(K,'projective');
    error(eid,msg);
end

Tvec(9) = 1;        % We assumed I = 1;

Tinv = reshape(Tvec,3,3);
T = inv(Tinv);

trans = maketform('projective', T);
output = [];

%-------------------- Function  findLinearConformalTransform ------------
function [trans, output] = findLinearConformalTransform(uv,xy,options)
% For a linear conformal transformation:
%
% let sc = s*cos(theta)
% let ss = s*sin(theta)
%
%                   [ sc -ss
% [u v] = [x y 1] *   ss  sc
%                     tx  ty]
%
% There are 4 unknowns: sc,ss,tx,ty.
%
% Another way to write this is:
%
% u = [x y 1 0] * [sc
%                  ss
%                  tx
%                  ty]
%
% v = [y -x 0 1] * [sc
%                   ss
%                   tx
%                   ty]
%
% With 2 or more correspondence points we can combine the u equations and
% the v equations for one linear system to solve for sc,ss,tx,ty.
%
% [ u1  ] = [ x1  y1  1  0 ] * [sc]
% [ u2  ]   [ x2  y2  1  0 ]   [ss]
% [ ... ]   [ ...          ]   [tx]
% [ un  ]   [ xn  yn  1  0 ]   [ty]
% [ v1  ]   [ y1 -x1  0  1 ]   
% [ v2  ]   [ y2 -x2  0  1 ]    
% [ ... ]   [ ...          ]
% [ vn  ]   [ yn -xn  0  1 ]
%
% Or rewriting the above matrix equation:
% U = X * r, where r = [sc ss tx ty]'
% so r = X\U.

K = options.K;
M = size(xy,1);
x = xy(:,1);
y = xy(:,2);
X = [x   y  ones(M,1)   zeros(M,1);
     y  -x  zeros(M,1)  ones(M,1)  ];

u = uv(:,1);
v = uv(:,2);
U = [u; v];

% We know that X * r = U
if rank(X) >= 2*K 
    r = X \ U;    
else
    eid = sprintf('Images:%s:twoUniquePointsReq',mfilename);
    error(eid,'At least 2 unique points needed to infer linear conformal transform.');
end

sc = r(1);      ss = r(2);
tx = r(3);      ty = r(4);

Tinv = [sc -ss 0;
        ss  sc 0;
        tx  ty 1];

T = inv(Tinv);
T(:,3) = [0 0 1]';

trans = maketform('affine', T);
output = [];

%-------------------- Function  findPolynomialTransform -------------------
function [trans, output] = findPolynomialTransform(uv,xy,options)
% For a polynomial transformation: 
%
% u = X*A, v = X*B, solve for A and B:
%     A = X\u;
%     B = X\v;   
% 
% The matrix X depends on the order of the polynomial.
% X will be M-by-K, where K = (order+1)*(order+2)/2;
% so A and B will be vectors of length K.
%
%   order = 2   
%     X = [ones(M,1),  x,  y,  x.*y,  x.^2,  y.^2];
%     so X is an M-by-6 matrix 
%
%   order = 3
%     X = [ones(M,1),  x,  y,  x.*y,  x.^2,  y.^2, ...
%          (x.^2).*y,  (y.^2).*x,  x.^3,  y.^3];
%     so X is an M-by-10 matrix 
%
%   order = 4
%     X = [ones(M,1),  x,  y,  x.*y,  x.^2,  y.^2, ...
%          (x.^2).*y,  (y.^2).*x,  x.^3,  y.^3, ...
%          (x.^3).*y,  (x.^2).*(y.^2),  x.*(y.^3),  x.^4,  y.^4];
%     so X is an M-by-15 matrix 
%
%

K = options.K;
order = options.order;

u = uv(:,1);    v = uv(:,2);
X = getTerms(order,xy);

if (rank(X) >= K)    % u = X*A, v = X*B, solve for A and B:
    A = X\u;
    B = X\v;   
else
    [msg, eid] = RankError(K,'polynomial');
    error(eid,msg);
end

trans = maketform('custom',2,2,[],@inv_polynomial,[A B]);
output = [];

%-------------------------------
function uv = inv_polynomial(xy,t)
% xy must be M-by-2
% t.tdata must be nterms-by-2, where nterms is 6, 10, or 15.

if size(xy,2) ~= 2
    eid = sprintf('Images:%s:xyMustHave2Cols',mfilename);
    error(eid,'XY must have two columns.');
end

if size(xy,1) < 1
    eid = sprintf('Images:%s:xyMustHaveOneRow',mfilename);
    error(eid,'XY must have at least one row.');
end

if size(t.tdata,2) ~= 2
    eid = sprintf('Images:%s:tdataMustHave2Cols',mfilename);
    error(eid,'TDATA must have two columns.');
end

nterms = size(t.tdata,1);

switch nterms
  case 6,    order = 2;
  case 10,   order = 3;
  case 15,   order = 4;
  otherwise
    eid = sprintf('Images:%s:tdataInvalidNumOfRows',mfilename);
    error(eid,'TDATA must have 6, 10, or 15 rows.');
end

X = getTerms(order,xy);

uv = X * t.tdata;

%-------------------------------
function X = getTerms(order,xy)

M = size(xy,1);
x = xy(:,1);        y = xy(:,2);

switch order
  case 2   
    X = [ones(M,1),  x,  y,  x.*y,  x.^2,  y.^2];
  case 3
    X = [ones(M,1),  x,  y,  x.*y,  x.^2,  y.^2, ...
         (x.^2).*y,  (y.^2).*x,  x.^3,  y.^3];
  case 4
    X = [ones(M,1),  x,  y,  x.*y,  x.^2,  y.^2, ...
         (x.^2).*y,  (y.^2).*x,  x.^3,  y.^3, ...
         (x.^3).*y,  (x.^2).*(y.^2),  x.*(y.^3),  x.^4,  y.^4];    
  otherwise
    eid = sprintf('Images:%s:invalidOrder',mfilename);
    error(eid,'ORDER must be 2, 3, or 4.');
end

%-------------------- Function  findPiecewiseLinear -----------------------
function [trans,output] = findPiecewiseLinear(uv,xy,options)
% For each triangle, finding a plane to warp xy into uv is an affine transformation.
%
% For an affine transformation:
%
%                   [ A D ]
% [u v] = [x y 1] * [ B E ]
%                   [ C F ]
%
% [ u1 v1 ]   [ x1 y1 1 ]   [ A D ]
% [ u2 v2 ] = [ x2 y2 1 ] * [ B E ]
% [ u3 v3 ]   [ x3 y3 1 ]   [ C F ]
%
% Rewriting the above matrix equation:
% U = X * T, where T = [A B C; D E F]'
%
% With the 3 correspondence points of each triangle, we can solve for T,
% T = X\U 
%
% see "Piecewise linear mapping functions for image registration" Ardeshir
% Goshtasby, Pattern Recognition, Vol 19, pp. 459-466, 1986.

% initialize
output.uv = uv;         output.xy = xy;
output.uv_dev = [];     output.xy_dev = [];
x = xy(:,1);            y = xy(:,2);

% Need to pass DELAUNAY options 'QJ' and 'Pp' to avoid qhull errors and 
% warnings for 2-D triangulations.
%tri = delaunay(x,y,{'QJ','Pp'}); 
tri = delaunayn([x(:) y(:)]);       % <- It will default to {'QJ','Pp'} options
ntri = size(tri,1);

if (ntri<2)
    [msg, eid] = RankError(options.K,'piecewise linear');
    error(eid,msg);
end

% Find all inside-out triangles
bad_triangles =  FindInsideOut(xy,uv,tri);
if ~isempty(bad_triangles)
    
    % find bad_vertices, eliminate bad_vertices
    num_bad_triangles = length(bad_triangles);
    bad_vertices = zeros(num_bad_triangles,1);
    for i = 1:num_bad_triangles
        bad_vertices(i) = FindBadVertex(x,y,tri(bad_triangles(i),:));
    end
    bad_vertices = unique(bad_vertices);
    num_bad_vertices = length(bad_vertices);    

    output.xy_dev = xy(bad_vertices,:); % update to return
    output.uv_dev = uv(bad_vertices,:);
    
    xy(bad_vertices,:) = []; % eliminate bad ones
    uv(bad_vertices,:) = [];
    nvert = size(xy,1);
    
    output.xy = xy; % update to return
    output.uv = uv;

    if (nvert < options.K)
        msg1 = sprintf('Eliminated %d control point pair(s).\n', num_bad_vertices);
        msg1 = [msg1 sprintf('Only %d control points remain.\n',nvert)];
        [eid,msg2] = CountError(options.K,'piecewise linear');
        error(eid,[msg1 msg2]);
    end
    x = xy(:,1);
    y = xy(:,2);
    tri = delaunay(x,y);
    ntri = size(tri,1);
    
    % see if there are any more bad triangles
    more_bad_triangles = FindInsideOut(xy,uv,tri);
    if ~isempty(more_bad_triangles)
        msg = sprintf('Eliminated %d control point pair(s).\n',...
                      num_bad_vertices);
        msg = [msg 'Fold-over triangles remain. See CP2TFORM reference page.'];
        eid = sprintf('Images:%s:foldoverTriangles',mfilename);
        error(eid,msg)
    end
    
    % print report about triangles and how many points were eliminated
    tri_list = sprintf('%d  ', bad_triangles);
    line1 =  ['Fold-over triangle(s):  ' tri_list sprintf('\n')];
    vert_list = sprintf('%d  ', bad_vertices);
    line2 = ['Compensated by eliminating control point pair(s):  ', vert_list];
    warnmsg = [line1 line2];
    warning('Images:cp2tform:foldOverTriangles', '%s', warnmsg)    
end

% calculate reverse mapping for each triangle
T = zeros(3,2,ntri);
for itri = 1:ntri
    X = [ xy( tri(itri,:), : ) ones(3,1)];
    U =   uv( tri(itri,:), : );
    if (rank(X) >= 3)
        T(:,:,itri) = X\U;
    else
        [msg, eid] = RankError(options.K,'piecewise linear');
        error(eid,msg);
    end
end

% Create TriangleGraph which is a sparse connectivity matrix
nxy = size(xy,1);
S = sparse( repmat((1:ntri)',1,3), tri, 1, ntri, nxy);

% Create OnHull to be 1 for ControlPoints on convex hull, 0 for
% interior points.
hull_indices = convhull(x,y);
OnHull = zeros(size(x));
OnHull(hull_indices) = 1;

tdata.Triangles = tri;
tdata.ControlPoints = xy;
tdata.OnHull = OnHull;
tdata.ConvexHullVertices = hull_indices;
tdata.TriangleGraph = S;
tdata.PiecewiseLinearTData = T;

trans = maketform('custom',2,2,[],@inv_piecewiselinear,tdata);

%-------------------- Function FindInsideOut ----------------------------
function index = FindInsideOut(xy,uv,tri)
% look for inside-out triangles using line integrals
x = xy(:,1);        y = xy(:,2);
u = uv(:,1);        v = uv(:,2);

p = size(tri,1);

xx = reshape(x(tri),p,3)';
yy = reshape(y(tri),p,3)';
xysum = sum( (xx([2 3 1],:) - xx).* (yy([2 3 1],:) + yy), 1 );

uu = reshape(u(tri),p,3)';
vv = reshape(v(tri),p,3)';
uvsum = sum( (uu([2 3 1],:) - uu).* (vv([2 3 1],:) + vv), 1 );
    
index = find(xysum.*uvsum<0);

%-------------------- Function FindBadVertex -----------------------------
function vertex = FindBadVertex(x,y,vertices)
% Get middle vertex of triangle where "middle" means the largest angle,
% which will have the smallest cosine.

vx = x(vertices)';      vy = y(vertices)';
abc = [ vx - vx([3 1 2]); vy - vy([3 1 2]) ];
a = abc(:,1);           b = abc(:,2);
c = abc(:,3);

% find cosine of angle between 2 vectors
vcos(1) = get_cos(-a, b);
vcos(2) = get_cos(-b, c);
vcos(3) = get_cos( a,-c);

[dum index] = min(vcos);
vertex = vertices(index);

%-------------------- Function get_cos -----------------------------------
function vcos = get_cos(a,b)
mag_a = sqrt( a(1)*a(1) + a(2)*a(2) );
mag_b = sqrt( b(1)*b(1) + b(2)*b(2) );
vcos = dot(a,b) / (mag_a*mag_b);


%-------------------- Function  findLWM ----------------------------------
function [trans,output] = findLWM(uv,xy,options)
% For a polynomial transformation: 
%
% u = X*A, v = X*B, solve for A and B:
%     A = X\u;
%     B = X\v;   
% 
% The matrix X depends on the order of the polynomial.
% X will be M-by-K, where K = (order+1)*(order+2)/2;
% so A and B will be vectors of length K.
%
%   order = 2   
%     X = [ones(M,1),  x,  y,  x.*y,  x.^2,  y.^2];
%     so X is an M-by-6 matrix 
%
% see "Image registration by local approximation methods" Ardeshir
% Goshtasby, Image and Vision Computing, Vol 6, p. 255-261, 1988.

if (options.order ~= 2)
    eid = sprintf('Images:%s:internalProblemPolyOrd',mfilename);
    error(eid,'Internal problem: only polynomials of order=2 are supported.');
end

output = [];
N = options.N;  M = size(xy,1);

x = xy(:,1);    y = xy(:,2);
u = uv(:,1);    v = uv(:,2);

T = zeros(options.K,2,M);
radii = zeros(M,1);
for icp = 1:M
    
    % find N closest points
    distcp = sqrt( (x-x(icp)).^2 + (y-y(icp)).^2 );
    [dist_sorted,indx] = sort(distcp);
    radii(icp) = dist_sorted(N);
    neighbors = indx(1:N);        
    neighbors = sort(neighbors);
    xcp = x(neighbors);
    ycp = y(neighbors);    
    ucp = u(neighbors);
    vcp = v(neighbors);    

    % set up matrix eqn for polynomial of order=2
    X = [ones(N,1),  xcp,  ycp,  xcp.*ycp,  xcp.^2,  ycp.^2];

    if rank(X)>=options.K
        % u = X*A, v = X*B, solve for A and B:
        A = X\ucp;
        B = X\vcp;
        T(:,:,icp) = [A B];
    else
        [msg, eid] = RankError(options.K,'polynomial');
        error(eid,msg);
    end

end

tdata.LWMTData = T;
tdata.ControlPoints = xy;
tdata.RadiiOfInfluence = radii;

trans = maketform('custom',2,2,[],@inv__lwm,tdata);

%-------------------- Function  ParseInputs -----------------------------
function [uv, xy, method, options] = cp2tform_ParseInputs(varargin)

% defaults
options.order = 3;      options.K = [];     N = [];
checknargin(2,4,nargin,mfilename);

% figure out if syntax is
% CP2TFORM(CPSTRUCT,TRANSFORMTYPE,...) or
% CP2TFORM(INPUT_POINTS,BASE_POINTS,TRANSFORMTYPE,...)

if isa(varargin{1},'struct')
    % TRANS = CP2TFORM(CPSTRUCT,TRANSFORMTYPE)
    % TRANS = CP2TFORM(CPSTRUCT,'polynomial',ORDER)        
    % TRANS = CP2TFORM(CPSTRUCT,'lwm',N)    
    checknargin(2,3,nargin,mfilename);
    
    [uv,xy] = cpstruct2pairs(varargin{1});
    method = GetMethod(varargin{2});

    nargs_to_go = nargin - 2;
    if (nargs_to_go > 0),   args = varargin(3:end);    end
else
    % TRANS = CP2TFORM(INPUT_POINTS,BASE_POINTS,TRANSFORMTYPE)
    % TRANS = CP2TFORM(INPUT_POINTS,BASE_POINTS,'polynomial',ORDER)        
    % TRANS = CP2TFORM(INPUT_POINTS,BASE_POINTS,'lwm',N)    

    checknargin(3,4,nargin,mfilename);
    uv = varargin{1};
    xy = varargin{2};
    method = GetMethod(varargin{3});

    nargs_to_go = nargin - 3;
    if (nargs_to_go > 0),   args = varargin(4:end);    end
end
    
if size(uv,2) ~= 2 || size(xy,2) ~= 2
    eid = sprintf('Images:%s:invalidControlPointMatrix', mfilename);
    error(eid,'In function %s, control point matrices must be M-by-2.', mfilename);
end

if size(uv,1) ~= size(xy,1)
    eid = sprintf('Images:%s:needSameNumControlPoints', mfilename);
    error(eid,'In function %s, INPUT and BASE images need same number of control points.', mfilename);
end

switch nargs_to_go
  case 0
    % TRANS = CP2TFORM(INPUT_POINTS,BASE_POINTS,TRANSFORMTYPE)
    % TRANS = CP2TFORM(CPSTRUCT,TRANSFORMTYPE)
  case 1
    if strcmp(method,'polynomial')
        % TRANS = CP2TFORM(INPUT_POINTS,BASE_POINTS,'polynomial',ORDER)        
        % TRANS = CP2TFORM(CPSTRUCT,'polynomial',ORDER)
        options.order = args{1};    
    elseif strcmp(method,'lwm')
        % TRANS = CP2TFORM(INPUT_POINTS,BASE_POINTS,'lwm',N)        
        % TRANS = CP2TFORM(CPSTRUCT,'lwm',N)                
        N = args{1};
    else 
        eid = sprintf('Images:%s:tooManyInputs', mfilename);
        error(eid,'Too many input arguments were passed to %s.', mfilename);
    end
  otherwise
    eid = sprintf('Images:%s:tooManyInputs', mfilename);
    error(eid,'Too many input arguments were passed to %s.', mfilename);
end

switch method
  case 'affine',            options.K = 3; 
  case 'projective',        options.K = 4;
  case 'linear conformal',  options.K = 2;
  case 'polynomial'
    order = options.order;
    % validate order
    if ~isnumeric(order) || numel(order)~=1 || sum(order==2:4)~=1
        eid = sprintf('Images:%s:invalidPolynomialOrder', mfilename);
        error(eid,'In function %s, polynomial order must be 2, 3, 4.', mfilename);
    end
    options.K = (order+1)*(order+2)/2;
  case 'piecewise linear'
    options.K = 4;      % for 'piecewise linear' need at least 4 points, 2 triangles
  case 'lwm'
    order = 2;
    options.K = (order+1)*(order+2)/2;
    options.order = order;

    if isempty(N)       % conservative default N protects user from ill-conditioned polynomials
        N = 2*options.K; 
    else                % validate N
        if ~isnumeric(N) || numel(N)~=1 || rem(N,1)~=0 || N<options.K
            eid = sprintf('Images:%s:invalidInputN', mfilename);
            error(eid,'In function %s, N must be an integer greater than or equal to %d.',...
                          mfilename, options.K);
        end
    end
    options.N = N;
  otherwise
    eid = sprintf('Images:%s:internalProblem', mfilename);
    error(eid,'Function %s has an internal problem: unrecognized method.', mfilename);
end    

%-------------------- Function  GetMethod ------------------------------------------
function method = GetMethod(method_string)

method_string = lower(method_string);

% Figure out which method to use
methods = {'affine', 'linear conformal', 'projective', 'polynomial', ...
           'piecewise linear', 'lwm'};
if ischar(method_string)
    indx = strmatch(method_string, methods);
    switch length(indx)
      case 0
        eid = sprintf('Images:%s:unrecognizedTransformType',mfilename);
        error(eid,'Unrecognized TRANSFORMTYPE ''%s'' in function %s.',method_string,mfilename);
      case 1
        method = methods{indx};
      otherwise
        eid = sprintf('Images:%s:ambiguousTransformType',mfilename);
        error(eid,'Ambiguous TRANSFORMTYPE ''%s'' in function %s.',method_string,mfilename);
    end
else
    eid = sprintf('Images:%s:transformTypeIsNotString',mfilename);
    error(eid,'TRANSFORMTYPE must be a string in function %s.',mfilename);
end	

%-------------------- Function  CountError ----------------------------------
function [msg, eid] = CountError(K,transform_string)
msg = sprintf('At least %d points needed to infer %s transform.', K,transform_string);
eid = sprintf('Images:%s:atLeast%dPointsReq',mfilename,K);

%-------------------- Function  RankError -----------------------------------
function [msg, eid] = RankError(K,transform_string)
msg = 'At least %d non-collinear points needed to infer %s transform.';
msg = sprintf(msg,K,transform_string);
eid = sprintf('Images:%s:atLeast%dNonColinearPointsReq',mfilename,K);

%----------------------------------------------------------------------------------
function checknargin(low, high, numInputs, function_name)
%CHECKNARGIN Check number of input arguments.
%   Copyright 1993-2003 The MathWorks, Inc.
%   $Revision: 1.3 $  $Date: 2003/01/17 16:28:15 $

% Input arguments are not checked for validity.
if numInputs < low
  msgId = sprintf('Images:%s:tooFewInputs', function_name);
  if low == 1
    msg1 = sprintf('Function %s expected at least 1 input argument', upper(function_name));
  else
    msg1 = sprintf('Function %s expected at least %d input arguments', upper(function_name), low);
  end
  
  if numInputs == 1
    msg2 = 'but was called instead with 1 input argument.';
  else
    msg2 = sprintf('but was called instead with %d input arguments.', numInputs);
  end
  
  error(msgId, '%s\n%s', msg1, msg2);
elseif numInputs > high
  msgId = sprintf('Images:%s:tooManyInputs', function_name);

  if high == 1
    msg1 = sprintf('Function %s expected at most 1 input argument', upper(function_name));
  else
    msg1 = sprintf('Function %s expected at most %d input arguments', upper(function_name), high);
  end
  
  if numInputs == 1
    msg2 = 'but was called instead with 1 input argument.';
  else
    msg2 = sprintf('but was called instead with %d input arguments.', numInputs);
  end
  error(msgId, '%s\n%s', msg1, msg2);
end

% -----------------------------------------------------------------------------------
function t = maketform( varargin )
% MAKETFORM Create spatial transformation structure (TFORM).  
%   T = MAKETFORM(TRANSFORMTYPE,...) creates a multidimensional spatial
%   transformation structure (a 'TFORM struct') that can be used with
%   TFORMFWD, TFORMINV, FLIPTFORM, IMTRANSFORM, or TFORMARRAY.
%   TRANSFORMTYPE can be 'affine', 'projective', 'custom', 'box', or
%   'composite'. Spatial transformations are also called geometric
%   transformations.
%
%   T = MAKETFORM('affine',A) builds a TFORM struct for an N-dimensional
%   affine transformation.  A is a nonsingular real (N+1)-by-(N+1) or
%   (N+1)-by-N matrix.  If A is (N+1)-by-(N+1), then the last column
%   of A must be [zeros(N,1); 1].  Otherwise, A is augmented automatically
%   such that its last column is [zeros(N,1); 1].  A defines a forward
%   transformation such that TFORMFWD(U,T), where U is a 1-by-N vector,
%   returns a 1-by-N vector X such that X = U * A(1:N,1:N) + A(N+1,1:N).
%   T has both forward and inverse transformations.
%
%   T = MAKETFORM('projective',A) builds a TFORM struct for an N-dimensional
%   projective transformation.  A is a nonsingular real (N+1)-by-(N+1)
%   matrix.  A(N+1,N+1) cannot be 0.  A defines a forward transformation
%   such that TFORMFWD(U,T), where U is a 1-by-N vector, returns a 1-by-N
%   vector X such that X = W(1:N)/W(N+1), where W = [U 1] * A.  T has
%   both forward and inverse transformations.
%   
%   T = MAKETFORM('affine',U,X) builds a TFORM struct for a
%   two-dimensional affine transformation that maps each row of U
%   to the corresponding row of X.  U and X are each 3-by-2 and
%   define the corners of input and output triangles.  The corners
%   may not be collinear.
%
%   T = MAKETFORM('projective',U,X) builds a TFORM struct for a
%   two-dimensional projective transformation that maps each row of U
%   to the corresponding row of X.  U and X are each 4-by-2 and
%   define the corners of input and output quadrilaterals.  No three
%   corners may be collinear.
%
%   T = MAKETFORM('custom',NDIMS_IN,NDIMS_OUT,FORWARD_FCN,INVERSE_FCN,
%   TDATA) builds a custom TFORM struct based on user-provided function
%   handles and parameters.  NDIMS_IN and NDIMS_OUT are the numbers of
%   input and output dimensions.  FORWARD_FCN and INVERSE_FCN are
%   function handles to forward and inverse functions.  Those functions
%   must support the syntaxes X = FORWARD_FCN(U,T) and U =
%   INVERSE_FCN(X,T), where U is a P-by-NDIMS_IN matrix whose rows are
%   points in the transformation's input space, and X is a
%   P-by-NDIMS_OUT matrix whose rows are points in the transformation's
%   output space.  TDATA can be any MATLAB array and is typically used to
%   store parameters of the custom transformation.  It is accessible to
%   FORWARD_FCN and INVERSE_FNC via the "tdata" field of T.  Either
%   FORWARD_FCN or INVERSE_FCN can be empty, although at least
%   INVERSE_FCN must be defined to use T with TFORMARRAY or IMTRANSFORM.
%
%   T = MAKETFORM('composite',T1,T2,...,TL) or T = MAKETFORM('composite',
%   [T1 T2 ... TL]) builds a TFORM whose forward and inverse functions
%   are the functional compositions of the forward and inverse functions
%   of the T1, T2, ..., TL.  For example, if L = 3 then TFORMFWD(U,T) is
%   the same as TFORMFWD(TFORMFWD(TFORMFWD(U,T3),T2),T1).  The components
%   T1 through TL must be compatible in terms of the numbers of input and
%   output dimensions.  T has a defined forward transform function only
%   if all of the component transforms have defined forward transform
%   functions.  T has a defined inverse transform function only if all of
%   the component functions have defined inverse transform functions.
%
%   T = MAKETFORM('box',TSIZE,LOW,HIGH) or T = MAKETFORM('box',INBOUNDS,
%   OUTBOUNDS) builds an N-dimensional affine TFORM struct, T.  TSIZE is
%   an N-element vector of positive integers, and LOW and HIGH are also
%   N-element vectors.  The transformation maps an input "box" defined
%   by the opposite corners ONES(1,N) and TSIZE or, alternatively, by
%   corners INBOUNDS(1,:) and INBOUND(2,:) to an output box defined by
%   the opposite corners LOW and HIGH or OUTBOUNDS(1,:) and OUTBOUNDS(2,:).
%   LOW(K) and HIGH(K) must be different unless TSIZE(K) is 1, in which
%   case the affine scale factor along the K-th dimension is assumed to be
%   1.0.  Similarly, INBOUNDS(1,K) and INBOUNDS(2,K) must be different
%   unless OUTBOUNDS(1,K) and OUTBOUNDS(1,K) are the same, and vice versa.
%   The 'box' TFORM is typically used to register the row and column
%   subscripts of an image or array to some "world" coordinate system.
%
%   Example
%   -------
%   Make and apply an affine transformation.
%
%       T = maketform('affine',[.5 0 0; .5 2 0; 0 0 1]);
%       tformfwd([10 20],T);
%       I = imread('cameraman.tif');
%       transformedI = imtransform(I,T);
%       figure, imshow(I), figure, imshow(transformedI)
%
%   See also FLIPTFORM, IMTRANSFORM, TFORMARRAY, TFORMFWD, TFORMINV.

%   Copyright 1993-2004 The MathWorks, Inc.
%   $Revision: 1.10.4.6 $ $Date: 2005/12/12 23:20:40 $

% Testing notes
% Syntaxes
%---------
% T = MAKETFORM( 'affine', A )
%
% A:        Numeric, non-singular, real square matrix (no Infs or Nans).
%           Last column must be be zero except for a one in the lower right corner.

msgstruct = nargchk(1,Inf,nargin);
if ~isempty(msgstruct)
    eid = sprintf('Images:%s:invalidNumInputs',mfilename);
    error(eid,'%s',msgstruct.message);
end

transform_type = getTransformType(varargin{1});

switch transform_type
    case 'affine',        fcn = @affine;
    case 'projective',    fcn = @projective;
    case 'composite',     fcn = @composite;
    case 'custom',        fcn = @custom;
    case 'box',           fcn = @box;
    otherwise
        eid = sprintf('Images:%s:unknownTransformType',mfilename);
        error(eid,'%s',sprintf('Unknown TRANSFORMTYPE: %s.', varargin{1}));
end

t = feval(fcn,varargin{2:end});

%--------------------------------------------------------------------------
function t = assigntform(ndims_in, ndims_out, forward_fcn, inverse_fcn, tdata)

% Use this function to ensure consistency in the way we assign
% the fields of each TFORM struct.

t.ndims_in    = ndims_in;
t.ndims_out   = ndims_out;
t.forward_fcn = forward_fcn;
t.inverse_fcn = inverse_fcn;
t.tdata       = tdata;

%--------------------------------------------------------------------------
function t = affine( varargin )
% Build an affine TFORM struct.

msgstruct = nargchk(1,2,nargin);
if ~isempty(msgstruct)
    eid = sprintf('Images:%s:invalidNumInputsForAffine',mfilename);    
    error(eid,'%s',msgstruct.message);
end

if nargin == 2
    % Construct a 3-by-3 2-D affine transformation matrix A
    % that maps the three points in X to the three points in U.
    U = varargin{1};
    X = varargin{2};
    A = construct_matrix( U, X, 'affine' ); 
    A(:,3) = [0 0 1]';  % Clean up noise before validating A.
else
    A = varargin{1};
end
A = validate_matrix( A, 'affine' );

N = size(A,2) - 1;
tdata.T    = A;
tdata.Tinv = inv(A);

% In case of numerical noise, coerce the inverse into the proper form.
tdata.Tinv(1:end-1,end) = 0;
tdata.Tinv(end,end) = 1;

t = assigntform(N, N, @fwd_affine, @inv_affine, tdata);

%--------------------------------------------------------------------------
function U = inv_affine( X, t )
% INVERSE affine transformation 
%
% T is an affine transformation structure. X is the row vector to
% be transformed, or a matrix with a vector in each row.

U = trans_affine(X, t, 'inverse');

%--------------------------------------------------------------------------
function X = fwd_affine( U, t)
% FORWARD affine transformation 
%
% T is an affine transformation structure. U is the row vector to
% be transformed, or a matrix with a vector in each row.

X = trans_affine(U, t, 'forward');

%--------------------------------------------------------------------------
function U = trans_affine( X, t, direction )
% Forward/inverse affine transformation method
%
% T is an affine transformation structure. X is the row vector to
% be transformed, or a matrix with a vector in each row.
% DIRECTION is either 'forward' or 'inverse'.

if strcmp(direction,'forward')
    M = t.tdata.T;
elseif strcmp(direction,'inverse')
    M = t.tdata.Tinv;
else
    eid = sprintf('Images:%s:invalidDirection',mfilename);
    error(eid,'%s','DIRECTION must be either ''forward'' or ''inverse''.');
end

X1 = [X ones(size(X,1),1)];   % Convert X to homogeneous coordinates
U1 = X1 * M;                  % Transform in homogeneous coordinates
U  = U1(:,1:end-1);           % Convert homogeneous coordinates to U

%--------------------------------------------------------------------------
function t = projective( varargin )
% Build a projective TFORM struct.

msgstruct = nargchk(1,2,nargin);
if ~isempty(msgstruct)
    eid = sprintf('Images:%s:invalidNumInputsForProjective',mfilename);    
    error(eid,'%s',msgstruct.message);
end

if nargin == 2
    % Construct a 3-by-3 2-D projective transformation matrix A
    % that maps the four points in U to the four points in X.
    U = varargin{1};
    X = varargin{2};
    A = construct_matrix( U, X, 'projective' ); 
else
    A = varargin{1};
end
A = validate_matrix( A, 'projective' );

N = size(A,2) - 1;
tdata.T    = A;
tdata.Tinv = inv(A);

t = assigntform(N, N, @fwd_projective, @inv_projective, tdata);

%--------------------------------------------------------------------------
function U = inv_projective( X, t )
% INVERSE projective transformation 
%
% T is an projective transformation structure. X is the row vector to
% be transformed, or a matrix with a vector in each row.
U = trans_projective(X, t, 'inverse');

%--------------------------------------------------------------------------
function X = fwd_projective( U, t)
% FORWARD projective transformation 
% T is an projective transformation structure. U is the row vector to
% be transformed, or a matrix with a vector in each row.
X = trans_projective(U, t, 'forward');

%--------------------------------------------------------------------------
function U = trans_projective( X, t, direction )
% Forward/inverse projective transformation method
%
% T is an projective transformation structure. X is the row vector to
% be transformed, or a matrix with a vector in each row.
% DIRECTION is either 'forward' or 'inverse'.

if strcmp(direction,'forward')
    M = t.tdata.T;
elseif strcmp(direction,'inverse')
    M = t.tdata.Tinv;
else
    eid = sprintf('Images:%s:invalidDirection',mfilename);
    error(eid,'%s','DIRECTION must be either ''forward'' or ''inverse''.');
end

N  = t.ndims_in;
X1 = [X ones(size(X,1),1)];   % Convert X to homogeneous coordinates
U1 = X1 * M;                  % Transform in homogeneous coordinates
UN = repmat(U1(:,end),[1 N]); % Replicate the last column of U
U  = U1(:,1:end-1) ./ UN;     % Convert homogeneous coordinates to U

%---------------------------------------------------------------
function A = validate_matrix( A, transform_type )

% Make sure A is double and real-valued.
if ~isa(A,'double') || ~isreal(A)
    eid = sprintf('Images:%s:invalidA',mfilename);
    error(eid,'%s','A must be a real-valued matrix of class double.');
end

% Make sure A is finite.
if ~all(isfinite(A))
    eid = sprintf('Images:%s:aContainsInfs',mfilename);
    error(eid,'%s','All elements of A must be finite.');
end

% Make sure A is (N + 1)-by-(N + 1).  Append a column if needed for 'affine'.
N = size(A,1) - 1;
if strcmp(transform_type,'affine') && size(A,2) == N
    A(:,N+1) = [zeros(N,1); 1];
end
if N < 1 || size(A,2) ~= N + 1
    eid = sprintf('Images:%s:invalidA',mfilename);
    error(eid,'%s','A must be square and at least 2-by-2.');
end

switch transform_type
    case 'affine'
      % Validate the final column of A.
      if any(A(:,N+1) ~= [zeros(N,1); 1])
           eid = sprintf('Images:%s:invalidAForAffine',mfilename);
           msg1 = 'The final column of A must consist of zeroes, except ';
           error(eid,'%s%s',msg1,'for a one in the last row.');
      end

    case 'projective'
      % Validate lower right corner of A
      if abs(A(N+1,N+1)) <= 100 * eps * norm(A)
            wid = sprintf('Images:%s:lastElementInANearZero',mfilename);
            warning(wid,'%s','The last element of A is very close to zero.');
      end
end

if cond(A) > 1e9
    wid = sprintf('Images:%s:conditionNumberofAIsHigh',mfilename);
    warning(wid,'%s',['The condition number of A is ' num2str(cond(A)) '.']);
end

%---------------------------------------------------------------
function A = construct_matrix( U, X, transform_type )
% Construct a 3-by-3 2-D transformation matrix A
% that maps the points in U to the points in X.

switch transform_type
  case 'affine'
    nPoints = 3;
    unitFcn = @UnitToTriangle;
  case 'projective'
    nPoints = 4;
    unitFcn = @UnitToQuadrilateral;
end

if any(size(U) ~= [nPoints 2])
    eid = sprintf('Images:%s:invalidU',mfilename);
    error(eid,'%s',sprintf('U must be %d-by-2.',nPoints));
end

if any(size(X) ~= [nPoints 2])
    eid = sprintf('Images:%s:invalidX',mfilename);
    error(eid,'%s',sprintf('X must be %d-by-2.',nPoints));
end

if ~isa(U,'double') || ~isreal(U) || ~all(isfinite(U(:)))
    eid = sprintf('Images:%s:invalidU',mfilename);
    error(eid,'%s','U must be real-valued (class double) and finite.');
end

if ~isa(X,'double') || ~isreal(X) || ~all(isfinite(X(:)))
    eid = sprintf('Images:%s:invalidX',mfilename);
    error(eid,'%s','X must be real-valued (class double) and finite.');
end

Au = feval(unitFcn,U);
if cond(Au) > 1e9
    wid = sprintf('Images:%s:conditionNumberOfUIsHigh',mfilename);
    warning(wid,'%s','Points in U are nearly collinear and/or coincidental.');
end

Ax = feval(unitFcn,X);
if cond(Ax) > 1e9
    wid = sprintf('Images:%s:conditionNumberOfXIsHigh',mfilename);
    warning(wid,'%s','Points in X are nearly collinear and/or coincidental.');
end

% (unit shape) * Au = U
% (unit shape) * Ax = X
%
% U * inv(Au) * Ax = (unit shape) * Ax = X and U * A = X,
% so inv(Au) * Ax = A, or Au * A = Ax, or A = Au \ Ax.

A = Au \ Ax;

if any(~isfinite(A(:)))
    eid = sprintf('Images:%s:collinearPointsinUOrX',mfilename);
    error(eid,'%s','Collinear points in U or X; cannot continue.');
end

A = A / A(end,end);

%---------------------------------------------------------------
function A = UnitToTriangle( X )

% Computes the 3-by-3 two-dimensional affine transformation
% matrix A that maps the unit triangle ([0 0], [1 0], [0 1])
% to a triangle with corners (X(1,:), X(2,:), X(3,:)).
% X must be 3-by-2, real-valued, and contain three distinct
% and non-collinear points. A is a 3-by-3, real-valued matrix.

A = [ X(2,1) - X(1,1)   X(2,2) - X(1,2)    0; ...
      X(3,1) - X(1,1)   X(3,2) - X(1,2)    0; ...
           X(1,1)            X(1,2)        1  ];

%---------------------------------------------------------------
function A = UnitToQuadrilateral( X )

% Computes the 3-by-3 two-dimensional projective transformation
% matrix A that maps the unit square ([0 0], [1 0], [1 1], [0 1])
% to a quadrilateral corners (X(1,:), X(2,:), X(3,:), X(4,:)).
% X must be 4-by-2, real-valued, and contain four distinct
% and non-collinear points.  A is a 3-by-3, real-valued matrix.
% If the four points happen to form a parallelogram, then
% A(1,3) = A(2,3) = 0 and the mapping is affine.
%
% The formulas below are derived in
%   Wolberg, George. "Digital Image Warping," IEEE Computer
%   Society Press, Los Alamitos, CA, 1990, pp. 54-56,
% and are based on the derivation in
%   Heckbert, Paul S., "Fundamentals of Texture Mapping and
%   Image Warping," Master's Thesis, Department of Electrical
%   Engineering and Computer Science, University of California,
%   Berkeley, June 17, 1989, pp. 19-21.

x = X(:,1);
y = X(:,2);

dx1 = x(2) - x(3);
dx2 = x(4) - x(3);
dx3 = x(1) - x(2) + x(3) - x(4);

dy1 = y(2) - y(3);
dy2 = y(4) - y(3);
dy3 = y(1) - y(2) + y(3) - y(4);

if dx3 == 0 && dy3 == 0
    % Parallelogram: Affine map
    A = [ x(2) - x(1)    y(2) - y(1)   0 ; ...
          x(3) - x(2)    y(3) - y(2)   0 ; ...
          x(1)           y(1)          1 ];
else
    % General quadrilateral: Projective map
    a13 = (dx3 * dy2 - dx2 * dy3) / (dx1 * dy2 - dx2 * dy1);
    a23 = (dx1 * dy3 - dx3 * dy1) / (dx1 * dy2 - dx2 * dy1);
    
    A = [x(2) - x(1) + a13 * x(2)   y(2) - y(1) + a13 * y(2)   a13 ;...
         x(4) - x(1) + a23 * x(4)   y(4) - y(1) + a23 * y(4)   a23 ;...
         x(1)                       y(1)                       1   ];
end

%---------------------------------------------------------------
function t = composite( varargin )
% Construct COMPOSITE transformation structure.

% Create TDATA as a TFORM structure array.
if nargin == 0
    eid = sprintf('Images:%s:tooFewTformStructs',mfilename);
    error(eid,'%s','At least one TFORM struct is required.');
elseif nargin == 1
    if length(varargin{1}) == 1
        % One TFORM input, just copy to t
        t = varargin{1};
        if ~istform(t)
            eid = sprintf('Images:%s:invalidTformStruct',mfilename);
            error(eid,'%s','T1 must be a TFORM struct.');
        end
        return;
    else
        % An array of TFORMs
        tdata = varargin{1};
        if ~istform(tdata)
            eid = sprintf('Images:%s:invalidTformStructArray',mfilename);
            error(eid,'%s','[T1, T2, ... TN] must be an array of TFORM structs.');
        end
    end
else
    % A list of TFORMs
    for k = 1:nargin
        if ~istform(varargin{k}) || length(varargin{k}) ~= 1
            eid = sprintf('Images:%s:invalidTformsInArray',mfilename);
            error(eid,'%s','T1, T2, ... TN must each be a TFORM struct.');
        end
    end
    tdata = [varargin{:}];
end

% Check for consistency of dimensions
N = length(tdata);
ndims_in =  [tdata.ndims_in];
ndims_out = [tdata.ndims_out];
if any(ndims_in(1:N-1) ~= ndims_out(2:N))
    eid = sprintf('Images:%s:tFormsDoNotHaveSameDimension',mfilename);
    error(eid,'%s','Input TFORM objects have inconsistent dimensions.');
end

% Check existence of forward and inverse function handles 
if any(cellfun('isempty',{tdata.forward_fcn}))
    forward_fcn = [];
else
    forward_fcn = @fwd_composite;
end

if any(cellfun('isempty',{tdata.inverse_fcn}))
    inverse_fcn = [];
else
    inverse_fcn = @inv_composite;
end

if (isempty(forward_fcn) && isempty(inverse_fcn))
    eid = sprintf('Images:%s:invalidForwardOrInverseFunction',mfilename);
    error(eid,'%s','Unable to compose either a forward or inverse transformation.');
end

t = assigntform(tdata(N).ndims_in, tdata(1).ndims_out, forward_fcn, inverse_fcn, tdata);

%---------------------------------------------------------------
function X = fwd_composite( U, t )
% FORWARD composite transformation 
% U is the row vector to be transformed, or a matrix with a vector in each row.
X = U;
for i = length(t.tdata):-1:1
    X = feval(t.tdata(i).forward_fcn, X, t.tdata(i));
end

%---------------------------------------------------------------
function U = inv_composite( X, t )
% INVERSE composite transformation 
% X is the row vector to be transformed, or a matrix with a vector in each row.
U = X;
for i = 1:length(t.tdata)
    U = feval(t.tdata(i).inverse_fcn, U, t.tdata(i));
end

%--------------------------------------------------------------------------
function t = custom( varargin )
msgstruct = nargchk(5,5,nargin);
if ~isempty(msgstruct)
    eid = sprintf('Images:%s:invalidNumInputsForCustom',mfilename);
    error(eid,'%s',msgstruct.message);
end

ndims_in    = varargin{1};      ndims_out   = varargin{2};
forward_fcn = varargin{3};      inverse_fcn = varargin{4};
tdata       = varargin{5};

% Validate sizes and types
if length(ndims_in) ~= 1 || ~isdoubleinteger(ndims_in)
    eid = sprintf('Images:%s:invalidNDims_In',mfilename);
    msg = 'NDIMS_IN must be a finite, integer-valued double.';
    error(eid,'%s',msg);
end

if length(ndims_out) ~= 1 || ~isdoubleinteger(ndims_out)
    eid = sprintf('Images:%s:invalidNDims_Out',mfilename);
    error(eid,'%s','NDIMS_OUT must be a finite, integer-valued double.');
end

if ndims_in < 1
    eid = sprintf('Images:%s:nDimsInIsNotPositive',mfilename);
    error(eid,'%s','NDIMS_IN must be positive.');
end

if ndims_out < 1
    eid = sprintf('Images:%s:nDimsOutIsNotPositive',mfilename);
    error(eid,'%s','NDIMS_OUT must be positive.');
end

if ~isempty(forward_fcn)
    if length(forward_fcn) ~= 1 || ~isa(forward_fcn,'function_handle')
        eid = sprintf('Images:%s:invalidForwardFcn',mfilename);
        error(eid,'%s','FORWARD_FCN must be a function handle.');
    end
end

if ~isempty(inverse_fcn)
    if length(inverse_fcn) ~= 1 || ~isa(inverse_fcn,'function_handle')
        eid = sprintf('Images:%s:invalidInverseFcn',mfilename);
        error(eid,'%s','INVERSE_FCN must be a function handle.');  
    end
end

if isempty(forward_fcn) && isempty(inverse_fcn)
    eid = sprintf('Images:%s:emptyFowardAndInverseFcn',mfilename);
    error(eid,'%s','FORWARD_FCN and INVERSE_FCN cannot both be empty.');
end

t = assigntform(ndims_in, ndims_out, forward_fcn, inverse_fcn, tdata);

%--------------------------------------------------------------------------
function t = box( varargin )
msgstruct = nargchk(2,3,nargin);
if ~isempty(msgstruct)
    eid = sprintf('Images:%s:invalidNumInputsForBox',mfilename);
    error(eid,'%s',msgstruct.message);
end

if nargin == 3
    % Construct an affine TFORM struct that maps a box bounded by 1 and TSIZE(k)
    % in dimension k to a box bounded by LO(k) and HI(k) in dimension k.
    % Construct INBOUNDS and OUTBOUNDS arrays, then call BOX2.
    
    tsize = varargin{1};    lo = varargin{2};   hi = varargin{3};
    tsize = tsize(:);       lo = lo(:);         hi = hi(:);
    
    if ~isdoubleinteger(tsize)
        eid = sprintf('Images:%s:invalidTSize',mfilename);
        error(eid,'%s','TSIZE must be integer-valued (class double) and finite.');
    end
    
    if any(tsize < 1 )
        eid = sprintf('Images:%s:tSizeIsNotPositive',mfilename);
        error(eid,'%s','All values in TSIZE must be greater than or equal to 1.');
    end
    
    if ~isa(lo,'double') || ~isreal(lo) || ~all(isfinite(lo))
        eid = sprintf('Images:%s:invalidLo',mfilename);
        error(eid,'%s','LO must be real-valued (class double) and finite.');
    end
    
    if ~isa(hi,'double') || ~isreal(hi) || ~all(isfinite(hi))
        eid = sprintf('Images:%s:invalidHi',mfilename);
        error(eid,'%s','HI must be real-valued (class double) and finite.');  
    end
    
    N = length(tsize);
    if length(lo) ~= N || length(hi) ~= N
        eid = sprintf('Images:%s:unequalLengthsForLoHiAndTSize',mfilename);
        error(eid,'%s','TSIZE, LO, and HI must be the same length.');  
    end

    if any(lo == hi & ~(tsize == 1))  % ok
        eid = sprintf('Images:%s:invalidLoAndHi',mfilename);
        msg1 = 'The corresponding entries in LO and HI may not be equal ';
        error(eid,'%s%s',msg1,'unless TSIZE is 1.');  
    end
    
    inbounds  = [ones(1,N); tsize'];
    outbounds = [lo'; hi'];
else
    inbounds  = varargin{1};
    outbounds = varargin{2};
end

t = box2(inbounds,outbounds);

%--------------------------------------------------------------------------
function t = box2( inBounds, outBounds )
% Construct an affine TFORM struct that maps a box bounded by INBOUNDS(1,k)
% and INBOUNDS(2,k) in dimensions k to a box bounded by OUTBOUNDS(1,k) and
% OUTBOUNDS(2,k).
%
% inBounds:   2-by-N
% outBounds:  2-by-N

if ~isfinitedouble(inBounds)
    eid = sprintf('Images:%s:invalidInbounds',mfilename);
    error(eid,'%s','INBOUNDS must be real and finite.');
end

if ~isfinitedouble(outBounds)
    eid = sprintf('Images:%s:invalidOutbounds',mfilename);
    error(eid,'%s','OUTBOUNDS must be real and finite.');
end

N = size(inBounds,2);
if (ndims(inBounds) ~= 2 || ndims(outBounds)  ~= 2 || size(inBounds,1) ~= 2 ...
   || size(outBounds,1) ~= 2 || size(outBounds,2) ~= N)
   eid = sprintf('Images:%s:invalidInboundsAndOutbounds',mfilename);
   error(eid,'%s','INBOUNDS and OUTBOUNDS must be 2-by-N (for the same N).');
end

qDegenerate  = (inBounds(1,:) == inBounds(2,:));
if any((outBounds(1,:) == outBounds(2,:)) ~= qDegenerate)
   eid = sprintf('Images:%s:invalidInboundsAndOutbounds',mfilename);
   msg1 = 'OUTBOUNDS(1,k) may equal OUTBOUNDS(2,k) if and only if ';
   msg2 = 'INBOUNDS(1,k) equals INBOUNDS(2,k).''INBOUNDS and OUTBOUNDS ';
   error(eid,'%s%s%s',msg1,msg2,'must be 2-by-N (for the same N).');
end

num = outBounds(2,:) - outBounds(1,:);
den =  inBounds(2,:)  - inBounds(1,:);

% Arbitrarily set the scale to unity for degenerate dimensions.
num(qDegenerate) = 1;
den(qDegenerate) = 1;

tdata.scale = num ./ den;
tdata.shift = outBounds(1,:) - tdata.scale .* inBounds(1,:);

t = assigntform(N, N, @fwd_box, @inv_box, tdata);

%--------------------------------------------------------------------------
function X = inv_box( X, t )
% INVERSE box transformation 
% T is an box transformation structure. X is the row vector to
% be transformed, or a matrix with a vector in each row.

M = size(X,1);
if (size(X,2) == 2 && numel(t.tdata.scale) == 2)    % Play safe. This is an order of mag faster
    X(:,1) = (X(:,1) - t.tdata.shift(1)) / t.tdata.scale(1);
    X(:,2) = (X(:,2) - t.tdata.shift(2)) / t.tdata.scale(2);
else        % This is the original code. HORRIBLY INEFFICIENT
	scale = repmat(t.tdata.scale, [M 1]);
	shift = repmat(t.tdata.shift, [M 1]);
	X = (X - shift) ./ scale;
end
% scale = repmat(t.tdata.scale, [M 1]);
% shift = repmat(t.tdata.shift, [M 1]);
% U = (X - shift) ./ scale;


%--------------------------------------------------------------------------
function U = fwd_box( U, t)
% FORWARD box transformation 
% T is an box transformation structure. U is the row vector to
% be transformed, or a matrix with a vector in each row.

M = size(U,1);
if (size(U,2) == 2 && numel(t.tdata.scale) == 2)    % Play safe. This is an order of mag faster
    U(:,1) = (U(:,1) * t.tdata.scale(1)) + t.tdata.shift(1);
    U(:,2) = (U(:,2) * t.tdata.scale(2)) + t.tdata.shift(2);
else
    scale = repmat(t.tdata.scale, [M 1]);
    shift = repmat(t.tdata.shift, [M 1]);
    U = (U .* scale) + shift;
end
% scale = repmat(t.tdata.scale, [M 1]);
% shift = repmat(t.tdata.shift, [M 1]);
% X = (U .* scale) + shift;

%--------------------------------------------------------------------------
function transform_type = getTransformType(type)

if ischar(type)
    low_type = lower(type);
else
    eid = sprintf('Images:%s:invalidTransformType',mfilename);
    error(eid,'%s','TRANSFORMTYPE must be a string.');
end
    
transform_names = {'affine','projective','composite','custom','box'};

% try to recognize the TransformType
imatch = strmatch(low_type,transform_names);

switch length(imatch)
  case 1    % one match
    transform_type = transform_names{imatch};
  case 0    % no matches
    eid = sprintf('Images:%s:unknownTransformType',mfilename);
    error(eid,'%s',sprintf('Unrecognized TRANSFORMTYPE ''%s''.',type));
  otherwise % more than one match
    eid = sprintf('Images:%s:ambiguousTransformType',mfilename);
    error(eid,'%s',sprintf('Ambiguous TRANSFORMTYPE ''%s''.', type));
end

%--------------------------------------------------------------------------
function q = isfinitedouble( x )
% Return true iff x is a finite, real-valued double array.
q = isa(x,'double');
if q
    q = isreal(x) && all(isfinite(x(:)));
end

%--------------------------------------------------------------------------
function q = isdoubleinteger( x )
% Return true iff x is a double array containing (real-valued and finite) integers.
if isa(x,'double')
    qint = (x == floor(x));
    q = ~(~isreal(x) || ~all(qint(:)) || ~all(isfinite(x(:))));
else
    q = ~1;
end

% -----------------------------------------------------------------------------------
function [B,xdata,ydata] = imtransform(varargin)
%IMTRANSFORM Apply 2-D spatial transformation to image.
%   B = IMTRANSFORM(A,TFORM) transforms the image A according to the 2-D
%   spatial transformation defined by TFORM, which is a tform structure
%   as returned by MAKETFORM or CP2TFORM.  If ndims(A) > 2, such as for
%   an RGB image, then the same 2-D transformation is automatically
%   applied to all 2-D planes along the higher dimensions.
%
%   When you use this syntax, IMTRANSFORM automatically shifts the origin of
%   your output image to make as much of the transformed image visible as
%   possible. If you are using IMTRANSFORM to do image registration, this syntax
%   is not likely to give you the results you expect; you may want to set
%   'XData' and 'YData' explicitly. See the description below of 'XData' and
%   'YData' as well as Example 3.
%
%   B = IMTRANSFORM(A,TFORM,INTERP) specifies the form of interpolation to
%   use.  INTERP can be one of the strings 'nearest', 'bilinear', or
%   'bicubic'.  Alternatively INTERP can be a RESAMPLER struct as returned
%   by MAKERESAMPLER.  This option allows more control over how resampling
%   is performed.  The default value for INTERP is 'bilinear'.
%
%   [B,XDATA,YDATA] = IMTRANSFORM(...) returns the location of the output
%   image B in the output X-Y space.  XDATA and YDATA are two-element
%   vectors.  The elements of XDATA specify the x-coordinates of the first
%   and last columns of B.  The elements of YDATA specify the y-coordinates
%   of the first and last rows of B.  Normally, IMTRANSFORM computes XDATA
%   and YDATA automatically so that B contains the entire transformed image
%   A.  However, you can override this automatic computation; see below.
%
%   [B,XDATA,YDATA] = IMTRANSFORM(...,PARAM1,VAL1,PARAM2,VAL2,...)
%   specifies parameters that control various aspects of the spatial
%   transformation. Parameter names can be abbreviated, and case does not
%   matter.
%
%   Parameters include:
%
%   'UData'      Two-element real vector.
%   'VData'      Two-element real vector.
%                'UData' and 'VData' specify the spatial location of the
%                image A in the 2-D input space U-V.  The two elements of
%                'UData' give the u-coordinates (horizontal) of the first
%                and last columns of A, respectively.  The two elements
%                of 'VData' give the v-coordinates (vertical) of the
%                first and last rows of A, respectively.  
%
%                The default values for 'UData' and 'VData' are [1
%                size(A,2)] and [1 size(A,1)], respectively.
%
%   'XData'      Two-element real vector.
%   'YData'      Two-element real vector.
%                'XData' and 'YData' specify the spatial location of the
%                output image B in the 2-D output space X-Y.  The two
%                elements of 'XData' give the x-coordinates (horizontal)
%                of the first and last columns of B, respectively.  The
%                two elements of 'YData' give the y-coordinates
%                (vertical) of the first and last rows of B,
%                respectively.  
%
%                If 'XData' and 'YData' are not specified, then
%                IMTRANSFORM estimates values for them that will
%                completely contain the entire transformed output image. 
%
%   'XYScale'    A one- or two-element real vector.
%                The first element of 'XYScale' specifies the width of
%                each output pixel in X-Y space.  The second element (if
%                present) specifies the height of each output pixel.  If
%                'XYScale' has only one element, then the same value is
%                used for both width and height.  
%
%                If 'XYScale' is not specified but 'Size' is, then
%                'XYScale' is computed from 'Size', 'XData', and 'YData'.
%                If neither 'XYScale' nor 'Size' is provided, then
%                the scale of the input pixels is used for 'XYScale'.
%
%   'Size'       A two-element vector of nonnegative integers.
%                'Size' specifies the number of rows and columns of the
%                output image B.  For higher dimensions, the size of B is
%                taken directly from the size of A.  In other words,
%                size(B,k) equals size(A,k) for k > 2.
%
%                If 'Size' is not specified, then it is computed from
%                'XData', 'YData', and 'XYScale'.
%
%   'FillValues' An array containing one or several fill values.
%                Fill values are used for output pixels when the
%                corresponding transformed location in the input image is
%                completely outside the input image boundaries.  If A is
%                2-D then 'FillValues' must be a scalar.  However, if A's
%                dimension is greater than two, then 'FillValues' can be
%                an array whose size satisfies the following constraint:
%                size(fill_values,k) must either equal size(A,k+2) or 1.
%                For example, if A is a uint8 RGB image that is
%                200-by-200-by-3, then possibilities for 'FillValues'
%                include:
%
%                    0                 - fill with black
%                    [0;0;0]           - also fill with black
%                    255               - fill with white
%                    [255;255;255]     - also fill with white
%                    [0;0;255]         - fill with blue
%                    [255;255;0]       - fill with yellow
%
%                If A is 4-D with size 200-by-200-by-3-by-10, then
%                'FillValues' can be a scalar, 1-by-10, 3-by-1, or
%                3-by-10.
%
%   Notes
%   -----
%   - When you do not specify the output-space location for B using
%     'XData' and 'YData', IMTRANSFORM estimates them automatically using
%     the function FINDBOUNDS.  For some commonly-used transformations,
%     such as affine or projective, for which a forward-mapping is easily
%     computable, FINDBOUNDS is fast.  For transformations that do not
%     have a forward mapping, such as the polynomial ones computed by
%     CP2TFORM, FINDBOUNDS can take significantly longer.  If you can
%     specify 'XData' and 'YData' directly for such transformations,
%     IMTRANSFORM may run noticeably faster.
%
%   - The automatic estimate of 'XData' and 'YData' using FINDBOUNDS is
%     not guaranteed in all cases to completely contain all the pixels of
%     the transformed input image.
%
%   - The output values XDATA and YDATA may not exactly equal the input
%     'XData and 'YData' parameters.  This can happen either because of
%     the need for an integer number or rows and columns, or if you
%     specify values for 'XData', 'YData', 'XYScale', and 'Size' that
%     are not entirely consistent.  In either case, the first element of
%     XDATA and YDATA always equals the first element of 'XData' and
%     'YData', respectively.  Only the second elements of XDATA and YDATA
%     might be different.
%
%   - IMTRANSFORM assumes spatial-coordinate conventions for the
%     transformation TFORM.  Specifically, the first dimension of the
%     transformation is the horizontal or x-coordinate, and the second
%     dimension is the vertical or y-coordinate.  Note that this is the
%     reverse of MATLAB's array subscripting convention.
%
%   - TFORM must be a 2-D transformation to be used with IMTRANSFORM.
%     For arbitrary-dimensional array transformations, see TFORMARRAY.
%
%   Class Support
%   -------------
%   A can be of any nonsparse numeric class, real or complex.  It can also be
%   logical.  The class of B is the same as the class of A.
%
%   Example 1
%   ---------
%   Apply a horizontal shear to an intensity image.
%
%       I = imread('cameraman.tif');
%       tform = maketform('affine',[1 0 0; .5 1 0; 0 0 1]);
%       J = imtransform(I,tform);
%       figure, imshow(I), figure, imshow(J)
%
%   Example 2
%   ---------
%   A projective transformation can map a square to a quadrilateral.  In
%   this example, set up an input coordinate system so that the input
%   image fills the unit square and then transform the image from the
%   quadrilateral with vertices (0 0), (1 0), (1 1), (0 1) to the
%   quadrilateral with vertices (-4 2), (-8 -3), (-3 -5), and (6 3).  Fill
%   with gray and use bicubic interpolation.  Make the output size the
%   same as the input size.
%
%       I = imread('cameraman.tif');
%       udata = [0 1];  vdata = [0 1];  % input coordinate system
%       tform = maketform('projective',[ 0 0;  1  0;  1  1; 0 1],...
%                                      [-4 2; -8 -3; -3 -5; 6 3]);
%       [B,xdata,ydata] = imtransform(I,tform,'bicubic','udata',udata,...
%                                                       'vdata',vdata,...
%                                                       'size',size(I),...
%                                                       'fill',128);
%       imshow(I,'XData',udata,'YData',vdata), axis on
%       figure, imshow(B,'XData',xdata,'YData',ydata), axis on
%
%   Example 3
%   ---------  
%   Register an aerial photo to an orthophoto.
%  
%       unregistered = imread('westconcordaerial.png');
%       figure, imshow(unregistered)
%       figure, imshow('westconcordorthophoto.png')
%       load westconcordpoints % load some points that were already picked     
%       t_concord = cp2tform(input_points,base_points,'projective');
%       info = imfinfo('westconcordorthophoto.png');
%       registered = imtransform(unregistered,t_concord,...
%                                'XData',[1 info.Width], 'YData',[1 info.Height]);
%       figure, imshow(registered)                       
%
%   See also CHECKERBOARD, CP2TFORM, IMRESIZE, IMROTATE, MAKETFORM, MAKERESAMPLER, TFORMARRAY.

%   Copyright 1993-2004 The MathWorks, Inc.
%   $Revision: 1.18.4.4 $  $Date: 2004/08/10 01:40:34 $

% Input argument details
% ----------------------
% A              numeric nonsparse array, any dimension, real or complex
%                may be logical
%                may be empty
%                NaN's and Inf's OK
%                required
%
% TFORM          valid TFORM struct as returned by MAKETFORM
%                checked using private/istform
%                required
%
% INTERP         one of these strings: 'nearest', 'linear', 'cubic'
%                case-insensitive match
%                nonambiguous abbreviations allowed
%
%                OR a resampler structure as returned by makeresampler
%                checked using private/isresample
%
%                optional; defaults to 'linear'
%
% 'FillValues'   double real matrix
%                Inf's and NaN's allowed
%                may be []
%                optional; defaults to 0
%                
% 'UData'        2-element real double vector
%                No Inf's or NaN's
%                The two elements must be different unless A has only
%                    one column
%                optional; defaults to [1 size(A,2)]
%
% 'VData'        2-element real double vector
%                No Inf's or NaN's
%                The two elements must be different unless A has only
%                    one row
%                optional; defaults to [1 size(A,1)]
%
% 'XData'        2-element real double vector
%                No Inf's or NaN's
%                optional; if not provided, computed using findbounds
%
% 'YData'        2-element real double vector
%                No Inf's or NaN's
%                optional; if not provided, computed using findbounds
%
% 'XYScale'      1-by-2 real double vector
%                elements must be positive
%                optional; default is the horizontal and vertical
%                  scale of the input pixels.
%
% 'Size'         real double row-vector
%                elements must be nonnegative integers
%                Can be 1-by-2 or 1-by-numdims(A).  If it is
%                  1-by-numdims(A), sizeB(3:end) must equal
%                  sizeA(3:end).
%                optional; default computation:
%                  num_rows = ceil(abs(ydata(2) - ydata(1)) ./ xyscale(2)) + 1;
%                  num_cols = ceil(abs(xdata(2) - xdata(1)) ./ xyscale(1)) + 1;
%
% If Size is provided and XYScale is not provided, then the output xdata
% and ydata will be the same as the input XData and YData.  Otherwise,
% the output xdata and ydata must be modified to account for the fact
% that the output size is contrained to contain only integer values.
% The first elements xdata and ydata is left alone, but the second
% values may be altered to make them consistent with Size and XYScale.

args = parse_inputs(varargin{:});
args.tform = make_composite_tform(args);

% imtransform uses x-y convention for ordering dimensions.
tdims_a = [2 1];        tdims_b = [2 1];

tsize_b = args.size([2 1]);
tmap_b = [];

B = tformarray(args.A, args.tform, args.resampler, tdims_a, tdims_b, ...
               tsize_b, tmap_b, args.fill_values);

xdata = args.xdata;     ydata = args.ydata;

%--------------------------------------------------
function [xdata,ydata] = recompute_output_bounds(args)

% For the purpose of this computation, change any 0's in out_size to 1's.
out_size = max(args.size, 1);
xdata = args.xdata;     ydata = args.ydata;

xdata(2) = args.xdata(1) + (out_size(2) - 1) .* args.xyscale(1) .* ...
    sign(args.xdata(2) - args.xdata(1));

ydata(2) = args.ydata(1) + (out_size(1) - 1) .* args.xyscale(2) .* ...
    sign(args.ydata(2) - args.ydata(1));

%--------------------------------------------------
function xyscale = compute_xyscale(args)

size_A = size(args.A);
xscale = compute_scale(args.udata, size_A(2));
yscale = compute_scale(args.vdata, size_A(1));

% If the output size would otherwise be twice the input size
% (in both dimensions), then multiply the output scale by a
% factor greater than 1.  This makes the output pixels larger
% (as measured in the output transform coordinates) and limits
% the size of the output image.

minscale = min(abs(args.xdata(2) - args.xdata(1)) / (2*size_A(2)),...
               abs(args.ydata(2) - args.ydata(1)) / (2*size_A(1)));

if xscale > minscale && yscale > minscale
  xyscale = [xscale yscale];
else
  xyscale = [xscale yscale] * minscale / max(xscale,yscale);
end

%--------------------------------------------------
function scale = compute_scale(udata, N)

scale_numerator = udata(2) - udata(1);
scale_denominator = max(N - 1, 0);
if scale_denominator == 0
    if scale_numerator == 0
        scale = 1;
    else
        eid = sprintf('Images:%s:unclearSpatialLocation',mfilename);
        error(eid,'%s','''UData'' or ''VData'' is inconsistent with size(A).');
    end
else
    scale = scale_numerator / scale_denominator;
end

%--------------------------------------------------
function new_tform = make_composite_tform(args)

reg_b = maketform('box', fliplr(args.size(1:2)), ...
                  [args.xdata(1) args.ydata(1)], ...
                  [args.xdata(2) args.ydata(2)]);

in_size = size(args.A);
in_size = in_size(1:2);

reg_a = maketform('box', fliplr(in_size), ...
					[args.udata(1) args.vdata(1)], ...
					[args.udata(2) args.vdata(2)]);

new_tform = maketform('composite', fliptform(reg_b), args.tform, reg_a);

%--------------------------------------------------
function args = parse_inputs(varargin)

%iptchecknargin(2,Inf,nargin,mfilename);

args.resampler = [];    args.fill_values = 0;
args.udata = [];        args.vdata = [];
args.xdata = [];        args.ydata = [];
args.size = [];         args.xyscale = [];

args.A = check_A(varargin{1});
args.tform = check_tform(varargin{2});

if rem(nargin,2) == 1    % IMTRANSFORM(A,TFORM,INTERP,<<prop/value pairs>>)
    args.resampler = check_resampler(varargin{3});
    first_prop_arg = 4;
else                    % IMTRANSFORM(A,TFORM,<<prop/value pairs>>)
    args.resampler = makeresampler('linear', 'fill');
    first_prop_arg = 3;
end

for k = first_prop_arg:2:nargin
    prop_string = check_property_string(varargin{k});
    switch prop_string
      case 'fillvalues',    args.fill_values = check_fill_values(varargin{k+1},args.A,args.tform);
      case 'udata',         args.udata = check_udata(varargin{k+1}, args.A);
      case 'vdata',         args.vdata = check_vdata(varargin{k+1}, args.A);
      case 'xdata',         args.xdata = check_xydata(varargin{k+1});
      case 'ydata',         args.ydata = check_xydata(varargin{k+1});
      case 'size',          args.size = check_size(varargin{k+1}, args.A);
      case 'xyscale',       args.xyscale = check_xyscale(varargin{k+1});
      otherwise
        msg = sprintf('Internal problem: unrecognized property string: %s', prop_string);
        error(sprintf('Images:%s:internalError',mfilename),'%s',msg);
    end
end

% Provide default values that require calculation.
if isempty(args.udata)
    args.udata = [1 size(args.A, 2)];
end

if isempty(args.vdata)
    args.vdata = [1 size(args.A, 1)];
end

if (isempty(args.xdata) + isempty(args.ydata) == 1)
    eid = sprintf('Images:%s:missingXDataAndYData',mfilename);
    msg1 = 'If either ''XData'' or ''YData'' is specified, ';
    error(eid,'%s%s',msg1,'then both must be specified.');
end

if isempty(args.xdata)
    % Output bounds not provided - estimate them.
    input_bounds = [args.udata(1) args.vdata(1);
                    args.udata(2) args.vdata(2)];
    try
        output_bounds = findbounds(args.tform, input_bounds);
        args.xdata = [output_bounds(1,1) output_bounds(2,1)];
        args.ydata = [output_bounds(1,2) output_bounds(2,2)];
    catch
        msg1 = 'XData and YData could not be automatically determined.';
        msg2 = 'Try specifying XData and YData explicitly in the call';
        msg3 = ' to IMTRANSFORM.';
        error(sprintf('Images:%s:unestimableOutputBounds',mfilename),'%s%s%s',msg1,msg2,msg3);
    end
end

if ~isempty(args.size)          % Output size was provided.
    if ~isempty(args.xyscale)   % xy_scale was provided; recompute bounds.
        [args.xdata,args.ydata] = recompute_output_bounds(args);
    else
        % Do nothing.  Scale was not provided but it is not needed.
    end

else                            % Output size was not provided.
    if isempty(args.xyscale)
        % Output scale was not provided.  Use the scale of the input pixels.
        args.xyscale = compute_xyscale(args);
    end
    
    % Compute output size.
    num_rows = ceil(abs(args.ydata(2) - args.ydata(1)) ./ args.xyscale(2)) + 1;
    num_cols = ceil(abs(args.xdata(2) - args.xdata(1)) ./ args.xyscale(1)) + 1;
    args.size = [num_rows, num_cols];    
    [args.xdata,args.ydata] = recompute_output_bounds(args);
end

%--------------------------------------------------
function A = check_A(A)
if (~isnumeric(A) && ~islogical(A)) || issparse(A)
    eid = sprintf('Images:%s:invalidImage',mfilename);
    error(eid,'%s','A must be a nonsparse numeric array.');
end

%--------------------------------------------------
function tform = check_tform(tform)
if ~istform(tform)
    eid = sprintf('Images:%s:invalidTFORM',mfilename);
    error(eid,'%s','Second input argument must be a valid TFORM structure.');
end

%--------------------------------------------------
function r = check_resampler(r)

if ischar(r)
    valid_strings = {'nearest','linear','cubic','bilinear','bicubic'};
    idx = strmatch(lower(r), valid_strings);
    switch length(idx)
      case 0
        eid = sprintf('Images:%s:unknownINTERP',mfilename);
        error(eid,'%s',sprintf('Unrecognized interpolation option: ''%s''', r));
      case 1
        r = valid_strings{idx};
        if (strcmp(r,'bicubic')),   r = 'cubic';    end
        if strcmp(r,'bilinear'),    r = 'linear';   end
      otherwise
        eid = sprintf('Images:%s:ambiguousINTERP',mfilename);
        error(eid,'%s',sprintf('Ambiguous interpolation option: ''%s''', r)); 
    end
    
    r = makeresampler(r, 'fill');
else
    if ~isresampler(r)
        msg1 = 'INTERP must be ''nearest'', ''linear'', ''cubic'',';
        msg2 = 'or a valid resampler as returned by MAKERESAMPLER.';
        error(sprintf('Images:%s:invalidResampler',mfilename),'%s%s', msg1, msg2);
    end
end


%--------------------------------------------------
function prop_string = check_property_string(prop_string)

if ~ischar(prop_string)
    eid = sprintf('Images:%s:invalidPropertyName',mfilename);
    error(eid,'%s','Property names must be strings.');
end

valid_strings = {'fillvalues'
                 'udata'
                 'vdata'
                 'xdata'
                 'ydata'
                 'xyscale'
                 'size'};

idx = strmatch(lower(prop_string), valid_strings);
switch length(idx)
 case 0
    eid = sprintf('Images:%s:unknownPropertyName',mfilename);
    error(eid,'%s',sprintf('Unrecognized property name: ''%s''', prop_string));    
  case 1
    prop_string = valid_strings{idx};
  otherwise
    eid = sprintf('Images:%s:ambiguousPropertyName',mfilename);
    error(eid,'%s',sprintf('Ambiguous property name: ''%s''', prop_string));
end

%--------------------------------------------------
function fill_values = check_fill_values(fill_values, A, tform)

N = tform.ndims_in;
size_F = size(fill_values);
size_A = size(A);
osize = size_A(N+1:end);
size_diff = length(osize) - length(size_F);
if size_diff < 0
    osize = [osize ones(1, -size_diff)];
else
    size_F = [size_F ones(1, size_diff)];
end

idx = find(size_F ~= osize);
if ~isempty(idx)
    if ~all(size_F(idx) == 1)
        error(sprintf('Images:%s:invalidFillValues',mfilename),'%s','Invalid size for ''FillValues''.');
    end
end

%--------------------------------------------------
function udata = check_udata(udata, A)

[m,n] = size(udata);

if ~(ndims(udata) == 2) || ~isa(udata,'double') || ...
        (m ~= 1) || (n ~= 2) || ~isreal(udata)
    eid = sprintf('Images:%s:invalidUData',mfilename);
    error(eid,'%s','''UData'' must be a real, double, 1-by-2 vector.');
end

if any(~isfinite(udata(:)))
    eid = sprintf('Images:%s:uDataContainsNansOrInfs',mfilename);
    error(eid,'%s','''UData'' must not contain NaN''s or Inf''s.');
end

if (udata(1) == udata(2)) && (size(A,2) ~= 1)
    eid = sprintf('Images:%s:uDataContainsInvalidElements',mfilename);
    error(eid,'%s','''UData'' elements cannot be equal unless A has only one column.');
end

%--------------------------------------------------
function vdata = check_vdata(vdata, A)

[m,n] = size(vdata);

if ~(ndims(vdata) == 2) || ~isa(vdata,'double') || ...
        (m ~= 1) || (n ~= 2) || ~isreal(vdata)
    eid = sprintf('Images:%s:invalidVData',mfilename);
    error(eid,'%s','''VData'' must be a real, double, 1-by-2 vector.');
end

if any(~isfinite(vdata(:)))
    eid = sprintf('Images:%s:vDataContainsNansOrInfs',mfilename);
    error(eid,'%s','''VData'' must not contain NaN''s or Inf''s.');
end

if (vdata(1) == vdata(2)) && (size(A,1) ~= 1)
    eid = sprintf('Images:%s:vDataContainsInvalidElements',mfilename);
    error(eid,'%s','''VData'' elements cannot be equal unless A has only one row.');
end

%--------------------------------------------------
function xdata = check_xydata(xdata)

[m,n] = size(xdata);
if ~(ndims(xdata) == 2) || ~isa(xdata,'double') || ...
        (m ~= 1) || (n ~= 2) || ~isreal(xdata)
    eid = sprintf('Images:%s:invalidXDataYData',mfilename);
    error(eid,'%s','''XData'' and ''YData'' must be real, double, 1-by-2 vectors.');
end

if any(~isfinite(xdata(:)))
    eid = sprintf('Images:%s:xDataYDataContainsNansOrInfs',mfilename);
    error(eid,'%s','''XData'' and ''YData'' must not contain NaN''s or Inf''s.');
end

%--------------------------------------------------
function output_size = check_size(output_size, A)
  
[m,n] = size(output_size);
size_A = size(A);
if ~isa(output_size,'double') || (m ~= 1) || ~isreal(output_size)
    eid = sprintf('Images:%s:invalidOutputSize',mfilename);
    error(eid,'%s','''Size'' must be real double row vector.');
end

if (n ~= 2) && n ~= ndims(A)
    eid = sprintf('Images:%s:invalidNumElementsInOutputSize',mfilename);
    error(eid,'%s','''Size'' must have either 2 or ndims(A) elements.');
end

if (n > 2) && ~isequal(output_size(3:end), size_A(3:end))
    eid = sprintf('Images:%s:invalidOutputSizeGreaterThan2D',mfilename);
    error(eid,'%s','''Size'' must equal size(A) in the 3rd and higher dimensions.');
end

if any(~isfinite(output_size(:)))
    eid = sprintf('Images:%s:outputSizeContainsNansOrInfs',mfilename);
    error(eid,'%s','''Size'' must not contain NaN''s or Inf''s.');
end

if any(floor(output_size(:)) ~= output_size(:)) || any(output_size(:) < 0)
    eid = sprintf('Images:%s:outputSizeContainsInvalidElements',mfilename);
    error(eid,'%s','''Size'' must contain nonnegative integers.');
end

%--------------------------------------------------
function xyscale = check_xyscale(xyscale)

[m,n] = size(xyscale);

if ~isa(xyscale,'double') || (m ~= 1) || ((n ~= 1) && (n ~= 2)) || ~isreal(xyscale)
    eid = sprintf('Images:%s:invalidXYScale',mfilename);
    error(eid,'%s','''XYScale'' must be a real, double, 1-by-1 or 1-by-2 vector.');
end

if any(~isfinite(xyscale(:)))
    eid = sprintf('Images:%s:xyScaleContainsNansOrInfs',mfilename);
    error(eid,'%s','''XYScale'' must not contain NaN''s or Inf''s.');
end

if any(xyscale(:) <= 0)
    eid = sprintf('Images:%s:xyScaleHasNegativeValues',mfilename);
    error(eid,'%s','''XYScale'' must contain positive values.');
end

if length(xyscale) == 1
    xyscale = [xyscale xyscale];
end

%--------------------------------------------------------------------------
function q = istform(t)
%ISTFORM True for valid geometric transformation structure.
%   ISTFORM(T) returns 1 if R is a valid geometric transformation structure,
%   such as one created by MAKETFORM, and 0 otherwise.

q = isa(t,'struct') & isfield(t,'ndims_in') ...
    & isfield(t,'ndims_out') & isfield(t,'forward_fcn') ...
    & isfield(t,'inverse_fcn') & isfield(t,'tdata');

% --------------------------------------------------------------------------
function B = tformarray( A, T, R, tdims_A, tdims_B, tsize_B, tmap_B, F )
%TFORMARRAY Apply spatial transformation to N-D array.
%   B = TFORMARRAY(A,T,R,TDIMS_A,TDIMS_B,TSIZE_B,TMAP_B,F) applies
%   a geometric transformation to array A to produce array B.
%
%   TFORMARRAY is like IMTRANSFORM, but is intended for problems
%   involving higher-dimensioned arrays or mixed input/output
%   dimensionality, or requiring greater user control or customization.
%   (Anything -- and more -- that can be accomplished with IMTRANSFORM
%   can be accomplished with a combination of MAKETFORM, MAKERESAMPLER,
%   FINDBOUNDS, and TFORMARRAY, but for many tasks involving 2-D images,
%   IMTRANSFORM is simpler.)
%
%   Brief description of inputs
%   ---------------------------
%   A         Input array or image
%   T         Geometric transformation, typically created with
%             MAKETFORM
%   R         Resampler, typically created with MAKERESAMPLER
%   TDIMS_A   Row vector listing the input transform dimensions
%   TDIMS_B   Row vector listing the output transform dimensions
%   TSIZE_B   Output array size in the transform dimensions
%   TMAP_B    Array of point locations in output space; can be used as an
%             alternative way to specify a geometric transformation
%   F         Array of fill values
%
%   Detailed description of inputs
%   ------------------------------
%   A         A can be any nonsparse numeric array, and can be real or
%             complex. It can also be logical.
%
%   T         T is a structure that defines a particular geometric
%             transformation.  For each location in the output transform
%             subscript space (as defined by TDIMS_B and TSIZE_B),
%             TFORMARRAY uses T and the function TFORMINV to compute the
%             corresponding location in the input transform subscript
%             space (as defined by TDIMS_A and SIZE(A)).
%
%             If T is empty, then TFORMARRAY operates as a direct resampling
%             function, applying the resampler defined in R to compute
%             values at each transform space location defined in TMAP_B 
%             (if TMAP_B is non-empty) or at each location in the output
%             transform subscript grid.
%
%   R         R is a structure that defines how to interpolate values of
%             the input array at specified locations.  R is usually
%             created with MAKERESAMPLER, which allows fine control
%             over how to interpolate along each dimension, as well as
%             what input array values to use when interpolating close the
%             edge of the array.
%
%   TDIMS_A   TDIMS_A and TDIMS_B indicate which dimensions of the input
%   TDIMS_B   and output arrays are involved in the geometric
%             transformation.  Each element must be unique, and must be a
%             positive integer.  The entries need not be listed in
%             increasing order, but the order matters.  It specifies the
%             precise correspondence between dimensions of arrays A and B
%             and the input and output spaces of the transformer, T.
%             LENGTH(TDIMS_A) must equal T.ndims_in, and LENGTH(TDIMS_B)
%             must equal T.ndims_out.
%
%             Suppose, for example, that T is a 2-D transformation,
%             TDIMS_A = [2 1], and TDIMS_B = [1 2].  Then the column
%             dimension and row dimension of A correspond to the first
%             and second transformation input-space dimensions,
%             respectively.  The row and column dimensions of B
%             correspond to the first and second output-space
%             dimensions, respectively.
%
%   TSIZE_B   TSIZE_B specifies the size of the array B along the
%             output-space transform dimensions.  Note that the size of B
%             along nontransform dimensions is taken directly from the
%             size of A along those dimensions.  If, for example, T is a
%             2-D transformation, size(A) = [480 640 3 10], TDIMS_B is
%             [2 1], and TSIZE_B is [300 200], then size(B) is [200 300 3].
%
%   TMAP_B    TMAP_B is an optional array that provides an alternative
%             way of specifying the correspondence between the position
%             of elements of B and the location in output transform
%             space.  TMAP_B can be used, for example, to compute the
%             result of an image warp at a set of arbitrary locations in
%             output space.  If TMAP_B is not empty, then the size of
%             TMAP_B takes the form: 
%
%                 [D1 D2 D3 ... DN L]
%
%             where N equals length(TDIMS_B).  The vector [D1 D2 ... DN]
%             is used in place of TSIZE_B.  If TMAP_B is not empty, then
%             TSIZE_B should be [].
%
%             The value of L depends on whether or not T is empty.  If T
%             is not empty, then L is T.ndims_out, and each L-dimension
%             point in TMAP_B is transformed to an input-space location
%             using T.  If T is empty, then L is LENGTH(TDIMS_A), and
%             each L-dimensional point in TMAP_B is used directly as a
%             location in input space.
%
%   F         F is a double-precision array containing fill values.
%             The fill values in F may be used in three situations:
%
%             (1) When a separable resampler is created with MAKERESAMPLER
%                 and its PADMETHOD is set to either 'fill' or 'bound',
% 
%             (2) When a custom resampler is used that supports the 'fill'
%                 or 'bound' pad methods (with behavior that is specific
%                 to the customization), or
% 
%             (3) When the map from the transform dimensions of B
%                 to the transform dimensions of A is deliberately
%                 undefined for some points. Such points are encoded in
%                 the input transform space by NaNs in either TMAP_B or
%                 in the output of TFORMINV.
%
%             In the first two cases, fill values are used to compute values
%             for output locations that map outside or near edges of the input
%             array.  Fill values are copied into B when output locations map
%             well outside the input array.  Type 'help makeresampler' for
%             further details on 'fill' and 'bound'.
%
%             F can be a scalar (including NaN), in which case its value
%             is replicated across all the nontransform dimensions.  Or F
%             can be a nonscalar whose size depends on size(A) in the
%             nontransform dimensions.  Specifically, if K is the J-th
%             nontransform dimension of A, then SIZE(F,J) must be either
%             SIZE(A,K) or 1.  As a convenience to the user, TFORMARRAY
%             replicates F across any dimensions with unit size such that
%             after the replication SIZE(F,J) equals size(A,K).
%
%             For example, suppose A represents 10 RGB images and has
%             size 200-by-200-by-3-by-10, T is a 2-D transformation, and
%             TDIMS_A and TDIMS_B are both [1 2].  In other words,
%             TFORMARRAY will apply the same 2-D transform to each color
%             plane of each of the 10 RGB images.  In this situation you
%             have several options for F:
%
%             - F can be a scalar, in which case the same fill value
%               is used for each color plane of all 10 images.
%
%             - F can be a 3-by-1 vector, [R G B]'.  Then R, G, and B
%               will be used as the fill values for the corresponding
%               color planes of each of the 10 images.  This can be
%               interpreted as specifying an RGB "fill color," with
%               the same color used for all 10 images.
%
%             - F can be a 1-by-10 vector.  This can be interpreted as
%               specifying a different fill value for each of 10
%               images, with that fill value being used for all three
%               color planes.  (Each image gets a distinct grayscale
%               fill-color.)
%
%             - F can be a 3-by-10 matrix, which can be interpreted as
%               supplying a different RGB fill-color for each of the
%               10 images.
%
%   Example
%   -------
%   Create a 2-by-2 checkerboard image where each square is 20 pixels
%   wide, then transform it with a projective transformation.  Use a
%   pad method of 'circular' when creating a resampler, so that the
%   output appears to be a perspective view of an infinite checkerboard.
%   Swap the output dimensions.  Specify a 100-by-100 output image.
%   Leave TMAP_B empty, since TSIZE_B is specified.  Leave the fill
%   value empty, since it won't be needed.
%
%       I = checkerboard(20,1,1);
%       figure, imshow(I)
%       T = maketform('projective',[1 1; 41 1; 41 41;   1 41],...
%                                  [5 5; 40 5; 35 30; -10 30]);
%       R = makeresampler('cubic','circular');
%       J = tformarray(I,T,R,[1 2],[2 1],[100 100],[],[]);
%       figure, imshow(J)
%
%   See also IMTRANSFORM, MAKERESAMPLER, MAKETFORM, FINDBOUNDS.

%   Copyright 1993-2004 The MathWorks, Inc.
%   $Revision: 1.9.4.6 $ $Date: 2005/12/12 23:21:20 $

% Start checking the inputs.
message = nargchk(8,8,nargin);
if ~isempty(message)
    error(sprintf('Images:%s:nargchkError',mfilename),'%s',message);
end

% Construct a new tsize_B if tmap_B is non-empty.
tsize_B = CheckDimsAndSizes( tdims_A, tdims_B, tsize_B, tmap_B, T );

% Get the 'full sizes' of A and B and their non-transform sizes (osize).
[fsize_A, fsize_B, osize] = fullsizes(size(A), tsize_B, tdims_A, tdims_B);

% Finish checking the inputs.
CheckInputArray( A );
CheckResampler( R, tdims_A );
F = CheckFillArray( F, osize );

% Determine blocking, if any.
if ~isempty(tmap_B)
    nBlocks = 1;  % Must process in a single block if tmap_B is supplied.
else
    blockingfactors = GetBlocking(tsize_B,max(length(tdims_A),length(tdims_B)));
    nBlocks = prod(blockingfactors);
end

% If there is no tmap_B, process large arrays in multiple blocks to conserve
% the memory required by the output grid G and its mapping to the input
% space M.  Otherwise, do the resampling in one large block.
if nBlocks == 1
    % If not already supplied in tmap_B, construct a grid G mapping B's
    % transform subscripts to T's output space or (if T = []) directly to
    % A's transform subscripts.
    if ~isempty(tmap_B)
        G = tmap_B;
    else
        hi = tsize_B;
        lo = ones(1,length(hi));
        G = ConstructGrid(lo,hi);
    end

    % If there is a tform, use it to extend the map in G to A's transform subscripts.
    if ~isempty(T)
        M = tforminv(G,T);
    else
        M = G;
    end
    G = []; % Free memory used by G (which is no longer needed).

    % Resample A using the map M and resampler R.
    B = resample(A, M, tdims_A, tdims_B, fsize_A, fsize_B, F, R);
else
    % Pre-allocate B with size fsize_B and class(B) = class(A).
    B(prod(fsize_B)) = A(1);
    B = reshape(B,fsize_B);
    
    % Loop over blocks in output transform space... 
    [lo, hi] = GetBlockLimits(tsize_B, blockingfactors);
    for i = 1:nBlocks
        % Construct the geometric map for the current block.
        G = ConstructGrid(lo(i,:),hi(i,:));
        if ~isempty(T)
            M = tforminv(G,T);
        else
            M = G;
        end
        G = []; % Free memory used by G (which is no longer needed).
    
        % Construct size and subscript arrays for the block, then resample.
        [bsize_B, S] = ConstructSubscripts(fsize_B, tdims_B, lo(i,:), hi(i,:));
        B(S{:}) = resample(A, M, tdims_A, tdims_B, fsize_A, bsize_B, F, R);
    end
end

%--------------------------------------------------------------------------
function B = resample( A, M, tdims_A, tdims_B, fsize_A, fsize_B, F, R )
% Evaluates the resampling function defined in the resampler R at
% each point in the subscript map M and replicated across all the
% non-transform dimensions of A, producing a warped version of A in B.

B = feval(R.resamp_fcn, A, M, tdims_A, tdims_B, fsize_A, fsize_B, F, R );

%--------------------------------------------------------------------------
function [bsize_B, S] = ConstructSubscripts( fsize_B, tdims_B, lo, hi )

% Determines the size of the block of array B bounded by the values in
% vectors LO and HI, and constructs an array of subscripts for assigning
% an array into that block.

bsize_B = fsize_B;
bsize_B(tdims_B) = hi - lo + 1;

fdims_B = length(fsize_B);
S = repmat({':'},[1,fdims_B]);
for i = 1:length(tdims_B)
	S{tdims_B(i)} = lo(i) : hi(i);
end

%--------------------------------------------------------------------------
function G = ConstructGrid( lo, hi )
% Constructs a regular grid from the range of subscripts defined by
% vectors LO and HI, and packs the result into a single array, G.
% G is D(1) x D(2) x ... x D(N) x N where D(k) is the length of
% lo(k):hi(k) and N is the length of LO and HI.
N = length(hi);     E = cell(1,N);
for i = 1:N;
	E{i} = lo(i):hi(i);  % The subscript range for each dimension
end
G = CombinedGrid(E{:});

%--------------------------------------------------------------------------
function G = CombinedGrid(varargin)
% If N >= 2, G = COMBINEDGRID(x1,x2,...,xN) is a memory-efficient
% equivalent to:
%
%   G = cell(N,1);
%   [G{:}] = ndgrid(x1,x2,...,xN);
%   G = cat(N + 1, G{:});
%
% (where x{i} = varargin{i} and D(i) is the number of
% elements in x{}).
%
% If N == 1, COMBINEDGRID returns x1 as a column vector. (The code
% with NDGRID replicates x1 across N columns, returning an N-by-N
% matrix -- which is inappropriate for our application.)
%
% N == 0 is not allowed.

N = length(varargin);
if N == 0
    eid = sprintf('Images:%s:combinegridCalledNoArgs',mfilename);
    error(eid,'%s','Internal Error: COMBINEDGRID called with no arguments.');
end

for i=1:N
    D(1,i) = numel(varargin{i});
end

if N == 1    % Special handling required to avoid calling RESHAPE with a single-element size vector.
    G = varargin{1}(:);
else    % Pre-allocate the output array.
    G = zeros([D N]);
    % Prepare to generate a comma-separated list of N colons.
    colons = repmat({':'},[1 N]);
    
    for i=1:N
        % Extract the i-th vector, reshape it to run along the
        % i-th dimension (adding singleton dimensions as needed),
        % replicate it across all the other dimensions, and
        % copy it to G(:,:,...,:,i) -- with N colons.
        x = varargin{i}(:);
        x = reshape( x, [ones(1,i-1) D(i) ones(1,N-i)] );
        x = repmat(  x, [ D(1:i-1)    1     D(i+1:N) ] );
        G(colons{:},i) = x;
    end
end

%--------------------------------------------------------------------------
function blockingfactors = GetBlocking( tsize, P )

% With large input arrays, the memory used by the grid and/or map
% (G and M) in the main function can be substantial and may exceed
% the memory used by the input and output images -- depending on the
% number of input transform dimensions, the number of non-transform
% dimensions, and the storage class of the input and output arrays.
% So the main function may compute the map block-by-block, resample
% the input for just that block, and assign the result to a pre-allocated
% output array.  We define the blocks by slicing the output transform space
% along lines, planes, or hyperplanes normal to each of the subscript
% axes.  We call the number of regions created by the slicing along a
% given axis the 'blocking factor' for that dimension.  (If no slices
% are made, the blocking factor is 1.  If all the blocking factors are
% 1 (which can be tested by checking if prod(blockingfactors) == 1),
% then the output array is small enough to process as a single block.)
% The blockingfactors array created by this function has the same
% size as TSIZE, and its dimensions are ordered as in TSIZE.
%
% This particular implementation has the following properties:
% - All blocking factors are powers of 2
% - After excluding very narrow dimensions, it slices across the
%   remaining dimensions to defined blocks whose dimensions tend
%   toward equality. (Of course, for large aspect ratios and relatively
%   few blocks, the result may not actually get very close to equality.)

% Control parameters
Nt = 20;  % Defines target block size
Nm =  2;  % Defines cut-off to avoid blocking on narrow dimensions
% The largest block will be on the order of 2^Nt doubles (but may be
% somewhat larger).  Any dimensions with size less than 2^Nm (= 4)
% will have blocking factors of 1.  They are tagged by zeros in the
% qDivisible array below.

% Separate the dimensions that are too narrow to divide into blocks.
qDivisible = tsize > 2^(Nm+1);
L = sum(qDivisible);

% The number of blocks will be 2^N.  As a goal, each block will
% contain on the order of 2^Nt values, but larger blocks may be
% needed to avoid subdividing narrow dimensions too finely.
%
% If all dimensions are large, then we'd like something like
%
%   (1) N = floor(log2(prod(tsize)/(2^Nt/P)));
%
% That's because we'd like prod(size(M)) / 2^N to be on the order of 2^Nt,
% and prod(size(M)) = prod(tsize) * P.  The following slightly more complex
% formula is equivalent to (1), but divides out the non-divisible dimensions
%
%   (2) N = floor(log2(prod(tsize(qDivisible)) / ...
%                             (2^Nt/(P*prod(tsize(~qDivisible)))))); 
%
% However, it is possible that we might not be able to block an image as
% finely as we'd like if its size is due to having many small dimensions
% rather than a few large ones.  That's why the actual formula for N in
% the next line of code replaces the numerator in (2) with 2^((Nm+1)*L)
% if this quantity is larger. In such a case, we'd have
%
%   (3) N = floor(log2(prod(tsize(qDivisible)) / 2^((Nm+1)*L)));
%
% and would have to make do with fewer blocks in order to ensure a minimum
% block size greater than or equal to 2^Nm.  The fact that our block sizes
% always satisfy this constraint is proved in the comments at the end of
% this function.

N = floor(log2(prod(tsize(qDivisible)) / ...
               max( (2^Nt)/(P*prod(tsize(~qDivisible))), 2^((Nm+1)*L) )));

% Initialize the blocking factor for the each divisible dimensions
% as unity.  Iterate N times, each time multiplying one of the
% blocking factors by 2.  The choice of which to multiply is driven
% by the goal of making the final set of average divisible block
% dimensions (stored in D at the end of the iteration) as uniform
% as possible.
B = ones(1,L);
D = tsize(qDivisible);
blockingfactors = zeros(size(tsize));

for i = 1:N
    k = find(D == max(D));
    k = k(1);  % Take the first if there is more than one maximum
    B(k) = B(k) * 2;
    D(k) = D(k) / 2;
end
blockingfactors( qDivisible) = B;
blockingfactors(~qDivisible) = 1;

% Assertion: After the loop is complete, all(D >= 2^Nm).
if any(D < 2^Nm)
    eid = sprintf('Images:%s:blockDimsTooSmall',mfilename);
    error(eid,'%s','Internal Error: Block dimension below minimum.');
end

% Let Dmin = min(D) after the loop is complete.  Dmin is the smallest
% average block size among the 'divisible' dimensions.  The following
% is a proof that Dmin >= 2^Nm.
%
% There are two possibilities: Either the blocking factor corresponding
% to Dmin is 1 or it is not.  If the blocking factor is 1, the proof
% is trival: Dmin > 2^(Nm+1) > 2^Nm (otherwise this dimension would
% not be included by virtue of qDivisible being false). Otherwise,
% because the algorithm continually divides the largest
% element of D by 2, it is clear that when the loop is complete
%
%   (1)  Dmin >= (1/2) * max(D).
%
% max(D) must equal or exceed the geometric mean of D,
%
%   (2)  max(D) >= (prod(D))^(1/L).
%
% From the formula for N,
%
%   (3)  N <= log2(prod(tsize(qDivisible)) / 2^((Nm+1)*L)).
%
% Because exactly N divisions by 2 occur in the loop,
%
%   (4)  prod(tsize(qDivisible) = (2^N) * prod(D);
%
% Combining (3) and (4),
%
%   (5)  prod(D) >= 2^((Nm+1)*L).
%
% Combining (1), (2), and (5) completes the proof,
%
%   (6)  Dmin >= (1/2) * max(D)
%             >= (1/2)*(prod(D))^(1/L)
%             >= (1/2)*(2^((Nm+1)*L))^(1/L) = 2^Nm.

%--------------------------------------------------------------------------
function [lo, hi] = GetBlockLimits( tsize_B, blockingfactors )
% Construct matrices LO and HI containing the lower and upper limits for
% each block, making sure to enlarge the upper limit of the right-most
% block in each dimension as needed to reach the edge of the image.
% LO and HI are each nBlocks-by-N.

N = length(tsize_B);
blocksizes = floor(tsize_B ./ blockingfactors);
blocksubs  = [0 ones( 1,N-1)];
delta      = [1 zeros(1,N-1)];
nBlocks = prod(blockingfactors);
lo = zeros(nBlocks,N);
hi = zeros(nBlocks,N);
for i = 1:nBlocks;
    blocksubs = blocksubs + delta;
    while any(blocksubs > blockingfactors)        % 'Carry' to higher dimensions as required.
        k = find(blocksubs > blockingfactors);
        blocksubs(k) = 1;
        blocksubs(k+1) = blocksubs(k+1) + 1;
    end
    lo(i,:) = 1 + (blocksubs - 1) .* blocksizes;
    hi(i,:) = blocksubs .* blocksizes;
    qUpperLim = blocksubs == blockingfactors;
    hi(i,qUpperLim) = tsize_B(qUpperLim);
end
if any(blocksubs ~= blockingfactors)
    eid = sprintf('Images:%s:failedToIterate',mfilename);
    error(eid,'%s','Internal Error: Failed to iterate over all blocks.');
end

%--------------------------------------------------------------------------
function [fsize_A, fsize_B, osize] = fullsizes(size_A, tsize_B, tdims_A, tdims_B)
% Constructs row vectors indicating the full sizes of A and B (FSIZE_A and
% FSIZE_B), including trailing singleton dimensions. Also constructs a row
% vector (OSIZE) listing (in order) the sizes of the non-transform
% dimensions of both A and B.
%
% There are two ways for trailing singletons to arise: (1) values in
% TDIMS_A exceed the length of SIZE_A or values in TDIMS_B exceed the
% length of TSIZE_B plus the number on non-transform dimensions of A and
% (2) all dimensions of A are transform dimensions (e.g., A is a grayscale
% image) -- in this case a trailing singleton is added to A and then
% transferred to B.
%
% Example:
%
%   [fsize_A, fsize_B] = ...
%      fullsizes( [7 512 512 512 3 20 ], [200 300], [2 3 4], [1 2] );
%
%   returns fsize_A = [  7   512   512   512     3    20]
%   and     fsize_B = [200   300     7     3    20]
%   and     osize   = [  7     3    20].

% Actual dimensionality of input array
ndims_A = length(size_A);

% 'Full' dimensionality of input array --
% Increase ndims(A) as needed, to allow for the largest
% transform dimensions as specified in tdims_A, then
% make sure there is at least one non-transform dimension.
fdims_A = max([ndims_A max(tdims_A)]); 
if fdims_A == length(tdims_A)
    fdims_A = fdims_A + 1;
end

% 'Full' size of input array --
% Pad size_A with ones (as needed) to allow for values
% in tdims_A that are higher than ndims(A):
fsize_A = [size_A ones(1,fdims_A - ndims_A)];

% The non-transform sizes of A and B:
osize = fsize_A(~ismember(1:fdims_A,tdims_A));

% The minimum ndims of B:
ndims_B = length(osize) + length(tdims_B);

% Increase ndims_B as needed, to allow for the largest
% transform dimensions as specified in tdims_B:
fdims_B = max(ndims_B, max(tdims_B)); 

% The full size of B, including possible padding:
isT_B = ismember(1:fdims_B, tdims_B);
padding = ones(1,fdims_B - (length(osize) + length(tsize_B)));
[sdims_B index] = sort(tdims_B);
fsize_B = zeros(size(isT_B));
fsize_B( isT_B) = tsize_B(index);
fsize_B(~isT_B) = [osize padding];

%---------------------------------------------------------------------
function K = PredictNDimsG( N, L )

if (N == 1),                K = 2;
elseif (N > 1 && L == 1)    K = N;
elseif (N > 1 && L > 1)     K = N + 1;
else
    eid = sprintf('Images:%s:NLMustBePositive',mfilename);
    error(eid,'%s','Internal Error: N and L must both be positive.');
end

%---------------------------------------------------------------------
function tsize_B = CheckDimsAndSizes( tdims_A, tdims_B, tsize_B, tmap_B, T )

% Check dimensions
if ~IsFinitePositiveIntegerRowVector( tdims_A )
    eid = sprintf('Images:%s:invalidTDIMS_A',mfilename);
    error(eid,'%s','TDIMS_A must be a row vector of finite, positive integers.');
end

if ~IsFinitePositiveIntegerRowVector( tdims_B )
    eid = sprintf('Images:%s:invalidTDIMS_B',mfilename);
    error(eid,'%s','TDIMS_B must be a row vector of finite, positive integers.');
end

P = length(tdims_A);
N = length(tdims_B);

if length(unique(tdims_A)) ~= P
    eid = sprintf('Images:%s:nonUniqueTDIMS_A',mfilename);
    error(eid,'%s','All values in TDIMS_A must be unique.');
end

if length(unique(tdims_B)) ~= N
    eid = sprintf('Images:%s:nonUniqueTDIMS_B',mfilename);
    error(eid,'%s','All values in TDIMS_B must be unique.');
end

% If tmap_B is supplied, ignore the input value of tsize_B and 
% construct tsize_B from tmap_B instead. Allow for the possibility
% that the last dimension of tmap_B is a singleton by copying no
% more than N values from size(tmap_B).

if isempty(tmap_B)
    L = N;
else
    if ~isempty(tsize_B)
        wid = sprintf('Images:%s:ignoringTSIZE_B',mfilename);       
        msg1 = 'Both TMAP_B and TSIZE_B are non-empty; ';
        warning(wid,'%s %s',msg1,'TSIZE_B will be ignored.');              
    end
    
    if ~isa(tmap_B,'double') || ~isreal(tmap_B) || issparse(tmap_B) || ...
          any(isinf(tmap_B(:)))
        eid = sprintf('Images:%s:invalidTMAP_B',mfilename);       
        msg1 = 'TMAP_B must be a non-sparse, finite real-valued array';
        error(eid,'%s %s',msg1, ' of class double.');
    end
    
    L = size(tmap_B,N+1);
    
    if isempty(T) && L ~= P
        eid = sprintf('Images:%s:invalidSizeTMAP_B',mfilename);        
        msg1 = 'SIZE(TMAP_B) is inconsistent with LENGTH(TDIMS_A) ';
        error(eid,'%s %s',msg1, 'or LENGTH(TDIMS_B).');
    end
    
    if ndims(tmap_B) ~= PredictNDimsG(N,L)
        eid = sprintf('Images:%s:invalidDimsTMAP_B',mfilename);        
        error(eid,'%s','NDIMS(TMAP_B) is inconsistent with TDIMS_B and SIZE(TMAP_B).');
    end
    
    size_G = size(tmap_B);
    tsize_B(1:N) = size_G(1:N);
end

% Check T
if ~isempty(T)
    if ~istform(T)
        eid = sprintf('Images:%s:invalidTFORMStruct',mfilename);        
        error(eid,'%s','T is not a valid TFORM struct.');
    end
    
    if length(T) > 1
        eid = sprintf('Images:%s:multipleTFORMStructs',mfilename);        
        error(eid,'%s','T must be a single TFORM struct.');
    end
    
    if T.ndims_in ~= P
        eid = sprintf('Images:%s:dimsMismatchA',mfilename);        
        error(eid,'%s','T.NDIMS_IN must match LENGTH(TDIMS_A).');
    end
    
    if T.ndims_out ~= L
        eid = sprintf('Images:%s:dimsMismatchB',mfilename);        
        error(eid,'%s','T.NDIMS_OUT is inconsistent with LENGTH(TDIMS_B) or SIZE(TMAP_B).');
    end
end

if ~IsFinitePositiveIntegerRowVector( tsize_B )
    eid = sprintf('Images:%s:invalidTSIZE_B',mfilename);        
    error(eid,'%s','TSIZE_B must be a row vector of finite, positive integers.');
end

% tsize_B and tdims_B must have the same length.
if length(tsize_B) ~= N
    eid = sprintf('Images:%s:lengthMismatchB',mfilename);        
    error(eid,'%s','Lengths of TSIZE_B and TDIMS_B must match.');
end

%---------------------------------------------------------------------
function CheckInputArray( A )
if (~isnumeric(A) && ~islogical(A)) || issparse(A)
    eid = sprintf('Images:%s:invalidA',mfilename);        
    error(eid,'%s','A must be a non-sparse numeric or logical array.');
end

%---------------------------------------------------------------------
function CheckResampler( R, tdims_A )

if ~isresampler(R) || (length(R) ~= 1)
    eid = sprintf('Images:%s:invalidResampler',mfilename);        
    error(eid,'%s','R must be a valid resampler struct.');
end

if R.ndims ~= Inf && R.ndims ~= length(tdims_A)
    eid = sprintf('Images:%s:resamplerDimsMismatch',mfilename);        
    error(eid,'%s','R''s ''ndims'' field doesn''t match LENGTH(TDIMS_A).');
end

%---------------------------------------------------------------------
function F = CheckFillArray( F, osize )

if isempty(F)
    F = 0;
else
    if ~isa(F,'double') || issparse(F)
        eid = sprintf('Images:%s:invalidF',mfilename);        
        error(eid,'%s','F must be a non-sparse array of class double.');
    end
    
    % Validate SIZE(F), stripping off trailing singletons.
    size_F = size(F);
    last = max([1 find(size_F ~= 1)]);
    size_F = size_F(1:last);
    
    % SIZE_F can't be longer than OSIZE.
    N = length(osize);
    q = (length(size_F) <= N);
    if q
        % Add (back) enough singletons to make size_F the same length as OSIZE.
        size_F = [size_F ones(1,N-length(size_F))];
        
        % Each value in SIZE_F must be unity (or zero), or must match OSIZE.
        q = all(size_F == 1 | size_F == osize);
    end
    if ~q
        eid = sprintf('Images:%s:sizeMismatchFA',mfilename);        
        error(eid,'%s','SIZE(F) is inconsistent with the non-tranform sizes of A.');
    end
end

%---------------------------------------------------------------------
function q = IsFinitePositiveIntegerRowVector( v )

q = isa(v,'double') & ~issparse(v) & ndims(v) == 2 & size(v,1) == 1 & size(v,2) >= 1;
if q
    q = ~(~isreal(v) | ~all(v == floor(v)) | ~all(isfinite(v)) | any(v < 1));
end

%--------------------------------------------------------------------------
function q = isresampler(R)
%ISRESAMPLER True for valid resampling structure.
%   ISRESAMPLER(R) returns 1 if R is a valid resampler struct, such as
%   one created by MAKERESAMPLER, and 0 otherwise.
% 
%   See also MAKERESAMPLER.
%   $Revision: 1.3.4.2 $ $Date: 2004/08/10 01:45:31 $

q = isa(R,'struct') & isfield(R,'ndims') & isfield(R,'padmethod') ...
    & isfield(R,'resamp_fcn') & isfield(R,'rdata');

if (q),    q = (length(R) == 1);    end
if (q),    q = ~isempty(R.ndims) & ~isempty(R.resamp_fcn) & ~isempty(R.padmethod);  end

if q
    q = isa(R.ndims,'double') & length(R.ndims) == 1    ...
        & isreal(R.ndims) & isa(R.padmethod,'char') ...
        & isa(R.resamp_fcn,'function_handle');
end

if (q && R.ndims ~= Inf)
    q = R.ndims == floor(R.ndims) & R.ndims >= 1;
end

if q
    q = (length(strmatch(R.padmethod,{'fill','bound','replicate','circular','symmetric'},'exact')) == 1);
end

% --------------------------------------------------------------------------------
function r = makeresampler( varargin )
%MAKERESAMPLER Create resampling structure.
%   R = MAKERESAMPLER(INTERPOLANT,PADMETHOD) creates a separable
%   resampler structure for use with TFORMARRAY and IMTRANSFORM.
%   In its simplest form, INTERPOLANT can be one of these strings:
%   'nearest', 'linear', or 'cubic'.  INTERPOLANT specifies the
%   interpolating kernel that the separable resampler uses.  PADMETHOD
%   can be one of these strings: 'replicate', 'symmetric', 'circular',
%   'fill', or 'bound'.  PADMETHOD controls how the resampler to
%   interpolates or assigns values to output elements that map close
%   to or outside the edge of input array.
%
%   PADMETHOD options
%   -----------------
%   In the case of 'fill', 'replicate', 'circular', or 'symmetric',
%   the resampling performed by TFORMARRAY or IMTRANSFORM occurs in
%   two logical steps: (1) pad A infinitely to fill the entire input
%   transform space, then (2) evaluate the convolution of the padded
%   A with the resampling kernel at the output points specified by
%   the geometric map.  Each non-transform dimension is handled
%   separately.  The padding is virtual, (accomplished by remapping
%   array subscripts) for performance and memory efficiency.
%  
%   'circular', 'replicate', and 'symmetric' have the same meanings as
%   in PADARRAY as applied to the transform dimensions of A:
%   
%     'replicate' -- Repeats the outermost elements
%     'circular'  -- Repeats A circularly
%     'symmetric' -- Mirrors A repeatedly.
%   
%   'fill' generates an output array with smooth-looking edges (except
%   when using nearest neighbor interpolation) because for output points
%   that map near the edge of the input array (either inside or outside),
%   it combines input image and fill values .
%  
%   'bound' is like 'fill', but avoids mixing fill values and input image
%   values.  Points that map outside are assigned values from the fill
%   value array.  Points that map inside are treated as with 'replicate'. 
%   'bound' and 'fill' produce identical results when INTERPOLANT is
%   'nearest'.
% 
%   It is up to the user to implement these behaviors in the case of a
%   custom resampler.
% 
%   Advanced options for INTERPOLANT
%   --------------------------------
%   In general, INTERPOLANT can have one of these forms:
%
%       1. One of these strings: 'nearest', 'linear', 'cubic'
%
%       2. A cell array: {HALF_WIDTH, POSITIVE_HALF}
%          HALF_WIDTH is a positive scalar designating the half width of
%          a symmetric interpolating kernel.  POSITIVE_HALF is a vector
%          of values regularly sampling the kernel on the closed interval
%          [0 POSITIVE_HALF].
%
%       3. A cell array: {HALF_WIDTH, INTERP_FCN}
%          INTERP_FCN is a function handle that returns interpolating
%          kernel values given an array of input values in the interval 
%          [0 POSITIVE_HALF].
%
%       4. A cell array whose elements are one of the three forms above.
%
%   Forms 2 and 3 are used to interpolate with a custom interpolating
%   kernel.  Form 4 is used to specify the interpolation method
%   independently along each dimension.  The number of elements in the
%   cell array for form 4 must equal the number of transform dimensions.
%   For example, if INTERPOLANT is {'nearest', 'linear', {2
%   KERNEL_TABLE}}, then the resampler will use nearest-neighbor
%   interpolation along the first transform dimension, linear
%   interpolation along the second, and a custom table-based
%   interpolation along the third.
%
%   Custom resamplers
%   -----------------
%   The syntaxes described above construct a resampler structure that
%   uses the separable resampler function that ships with the Image
%   Processing Toolbox.  It is also possible to create a resampler
%   structure that uses a user-written resampler by using this syntax: 
%   R = MAKERESAMPLER(PropertyName,PropertyValue,...).  PropertyName can
%   be 'Type', 'PadMethod', 'Interpolant', 'NDims', 'ResampleFcn', or
%   'CustomData'.
%
%   'Type' can be either 'separable' or 'custom' and must always be
%   supplied.  If 'Type' is 'separable', the only other properties that can
%   be specified are 'Interpolant' and 'PadMethod', and the result is
%   equivalent to using the MAKERESAMPLER(INTERPOLANT,PADMETHOD) syntax.
%   If 'Type' is 'custom', then 'NDims' and 'ResampleFcn' are required
%   properties, and 'CustomData' is optional.  'NDims' is a positive
%   integer and indicates what dimensionality the custom resampler can
%   handle.  Use a value of Inf to indicate that the custom resampler can
%   handle any dimension.  The value of 'CustomData' is unconstrained.
%
%   'ResampleFcn' is a handle to a function that peforms the resampling.
%   The function will be called with the following interface:
%
%       B = RESAMPLE_FCN(A,M,TDIMS_A,TDIMS_B,FSIZE_A,FSIZE_B,F,R)
%
%   See the help for TFORMARRAY for information on the inputs A, TDIMS_A,
%   TDIMS_B, and F.
%
%   M is an array that maps the transform subscript space of B to the
%   transform subscript space of A.  If If A has N transform dimensions (N =
%   length(TDIMS_A)) and B has P transform dimensions (P = length(TDIMS_B)),
%   then NDIMS(M) = P + 1 if N > 1 and P if N == 1, and SIZE(M, P + 1) =
%   N.  The first P dimensions of M correspond to the output transform
%   space, permuted according to the order in which the the output
%   transform dimensions are listed in TDIMS_B.  (In general TDIMS_A and
%   TDIMS_B need not be sorted in ascending order, although such a
%   limitation may be imposed by specific resamplers.)  Thus the first P
%   elements of SIZE(M) determine the sizes of the transform dimensions of
%   B.  The input transform coordinates to which each point is mapped are
%   arrayed across the final dimension of M, following the order given in
%   TDIMS_A.  M must be double.
%
%   FSIZE_A and FSIZE_B are the full sizes of A and B, padded with 1s as
%   necessary to be consistent with TDIMS_A, TDIMS_B, and SIZE(A).
%
%   Example
%   -------
%   Stretch an image in the y-direction using separable resampler that
%   applies in cubic interpolation in the y-direction and nearest
%   neighbor interpolation in the x-direction. (This is equivalent to,
%   but faster than, applying bicubic interpolation.)
%
%       A = imread('moon.tif');
%       resamp = makeresampler({'nearest','cubic'},'fill');
%       stretch = maketform('affine',[1 0; 0 1.3; 0 0]);
%       B = imtransform(A,stretch,resamp);
%
%   See also IMTRANSFORM, TFORMARRAY.

%   Copyright 1993-2003 The MathWorks, Inc.
%   $Revision: 1.5 $ $Date: 2003/01/17 16:27:45 $

message = nargchk(2,10,nargin);
if ~isempty(message),    error(message);    end

if mod(nargin,2) ~= 0
    error('Invalid number of arguments.');
end
npairs = nargin / 2;

property_strings = {'type','padmethod','interpolant','ndims','resamplefcn','customdata'};

% Check for the shorthand syntax for separable resamplers.
if npairs == 1
    if ischar(varargin{1})
        if isempty(FindValue('type', property_strings, varargin{:}))
            r = MakeSeparable( varargin{1}, varargin{2} );
            return;
        end
    else
        r = MakeSeparable( varargin{1}, varargin{2} );
        return;
    end
end

% Parse property name/property value syntax.
type = FindValue('type', property_strings, varargin{:});
if isempty(type)
    error('Missing ''Type'' property.');
end
switch GetCanonicalString(type,'Type',{'separable','custom'})
  case 'separable'
    [interpolant, padmethod] = ParseSeparable( property_strings, varargin{:} );   
    r = MakeSeparable( interpolant, padmethod );
  case 'custom'
    r = MakeCustom( property_strings, varargin{:} );
  otherwise
    error('Internal Error: ''Type'' must be ''separable'' or ''custom''.');
end

%--------------------------------------------------------------------------
function r = MakeCustom( property_strings, varargin )

CheckPropertyNames('Custom',...
  {'type','padmethod','ndims','resamplefcn','customdata'}, varargin{:});

padmethod = FindValue('padmethod', property_strings, varargin{:});
if isempty(padmethod)
    padmethod = 'replicate';
end
padmethod = GetCanonicalString( padmethod, 'PadMethod', ...
                  {'fill','bound','replicate','circular','symmetric'});

n_dimensions = FindValue('ndims', property_strings, varargin{:});
if isempty(n_dimensions)
    error('The ''NDims'' property must be specified.');
end
if length(n_dimensions) ~= 1
    error('''NDims'' must have a scalar value.');
end
if ~isa(n_dimensions,'double') || ~isreal(n_dimensions)
    error('''NDims'' must have numeric, real value.');
end
if n_dimensions ~= floor(n_dimensions) || n_dimensions < 1
    error('''NDims'' must have a positive integer value.');
end

resample_fcn = FindValue('resamplefcn', property_strings, varargin{:});
if isempty(resample_fcn)
    error('The ''ResampleFcn'' property must be specified.');
end
if length(resample_fcn) ~= 1
    error('''ResampleFcn'' must have a scalar value.');
end
if ~isa(resample_fcn,'function_handle')
    error('''ResampleFcn'' must have function handle value.');
end

rdata = FindValue('customdata', property_strings, varargin{:});
    
r = AssignResampler(n_dimensions, padmethod, resample_fcn, rdata);

%--------------------------------------------------------------------------
function [interpolant, padmethod] = ParseSeparable(property_strings, varargin)

CheckPropertyNames('Separable',...
  {'type','interpolant','padmethod'}, varargin{:});

interpolant = FindValue('interpolant', property_strings, varargin{:});
if isempty(interpolant)
    interpolant = 'linear';
end 
           
padmethod = FindValue('padmethod', property_strings, varargin{:});
if isempty(padmethod)
    padmethod = 'replicate';
end

%--------------------------------------------------------------------------
function r = MakeSeparable( interpolant, padmethod )

standardFrequency = 1000;  % Standard number of samples per unit 
n_dimensions = inf;
              
if isa(interpolant,'cell')
    if HasCustomTable(interpolant)
        rdata.K = interpolant;
    elseif HasCustomFunction(interpolant)
        rdata.K = CustomKernel(interpolant, standardFrequency);
    else
        n_dimensions = length(interpolant);
        rdata.K = MultipleKernels(interpolant, standardFrequency);
    end
else
    rdata.K = StandardKernel(interpolant, standardFrequency);
end

padmethod = GetCanonicalString( padmethod, 'PadMethod', ...
                  {'fill','bound','replicate','circular','symmetric'});

r = AssignResampler(n_dimensions, padmethod, @resampsep, rdata);

%--------------------------------------------------------------------------
function q = HasCustomTable(interpolant)

q = isa(interpolant,'cell');

if q
    q = length(interpolant) == 2;
end

if q
    q = isa(interpolant{1},'double') ...
        & length(interpolant{1}) == 1 ...
        & isreal(interpolant{1}) ...
        & isa(interpolant{2},'double') ...
        & isreal(interpolant{2});
end

if q
    q = interpolant{1} > 0;
end

%--------------------------------------------------------------------------
function q = HasCustomFunction(interpolant)

q = isa(interpolant,'cell');

if q
    q = length(interpolant) == 2;
end

if q
    q = isa(interpolant{1},'double') ...
        & length(interpolant{1}) == 1 ...
        & isreal(interpolant{1}) ...
        & isa(interpolant{2},'function_handle');
end

if q
    q = interpolant{1} > 0;
end

%--------------------------------------------------------------------------
function K = MultipleKernels( interpolant, frequency )

for i = 1:length(interpolant)
    if HasCustomTable(interpolant{i})
        K{i} = interpolant{i};
    elseif HasCustomFunction(interpolant{i})
        K{i} = CustomKernel(interpolant{i}, frequency);
    else
        K{i} = StandardKernel(interpolant{i}, frequency);
    end
end

%--------------------------------------------------------------------------
function K = CustomKernel( interpolant, frequency )

halfwidth  = interpolant{1};
kernel_fcn = interpolant{2};
positiveHalf = SampleKernel(kernel_fcn, halfwidth, frequency);
K = {halfwidth, positiveHalf};

%--------------------------------------------------------------------------
function K = StandardKernel( interpolant, frequency )

interpolant = GetCanonicalString( interpolant, 'Interpolant', ...
                                  {'nearest','linear','cubic'});
switch interpolant
  case 'nearest'
    K = [];
  case 'linear' 
    halfwidth = 1.0;
    positiveHalf = SampleKernel(@LinearKernel, halfwidth, frequency);
    K = {halfwidth, positiveHalf};
  case 'cubic'
    halfwidth = 2.0;
    positiveHalf = SampleKernel(@CubicKernel, halfwidth, frequency);
    K = {halfwidth, positiveHalf};
  otherwise
    error('Internal Error: INTERPOLANT must be ''nearest'', ''linear'', or ''cubic''.');
end

%--------------------------------------------------------------------------
function positiveHalf = SampleKernel( kernel, halfwidth, frequency )

if length(kernel) ~= 1 || ~isa(kernel,'function_handle')
    error('kernel must be a function handle.');
end
n = floor(halfwidth * frequency);
positiveHalf = feval( kernel, (halfwidth / n) * (0:n) ); 

%--------------------------------------------------------------------------
function y = LinearKernel( x )

y = zeros(1,length(x));
reshape(y,size(x));
x(x < 0) = -x(x < 0);
q = (x <= 1);
y(q) = 1 - x(q);

%--------------------------------------------------------------------------
function y = CubicKernel( x )

% There is a whole family of "cubic" interpolation kernels. The 
% particular kernel used here is described in the article Keys,
% "Cubic Convolution Interpolation for Digital Image Processing,"
% IEEE Transactions on Acoustics, Speech, and Signal Processing,
% Vol. ASSP-29, No. 6, December 1981, p. 1155.

y = zeros(1,length(x));
reshape(y,size(x));
x(x < 0.0) = -x(x < 0.0);

q = (x <= 1);            % Coefficients: 1.5, -2.5, 0.0, 1.0
y(q) = ((1.5 * x(q) - 2.5) .* x(q)) .* x(q) + 1.0;

q = (1 < x & x <= 2);    % Coefficients: -0.5, 2.5, -4.0, 2.0
y(q) = ((-0.5 * x(q) + 2.5) .* x(q) - 4.0) .* x(q) + 2.0;

%--------------------------------------------------------------------------
function r = AssignResampler(n_dimensions, padmethod, resamp_fcn, rdata)

% Use this function to ensure consistency in the way we assign
% the fields of each resampling struct. Note that r.ndims = Inf
% is used to denote that the resampler supports an arbitrary
% number of dimensions.

r.ndims      = n_dimensions;
r.padmethod  = padmethod;
r.resamp_fcn = resamp_fcn;
r.rdata      = rdata;

%--------------------------------------------------------------------------
function value = FindValue( property_name, property_strings, varargin )

value = [];
for i = 1:((nargin-2)/2)
    current_name = varargin{2*i-1};
    if ischar(current_name)
        imatch = strmatch(lower(current_name),property_strings);
        nmatch = length(imatch);
        if nmatch > 1
                error(['Ambiguous property name ' current_name '.']);
        end
        if nmatch == 1
            canonical_name = property_strings{imatch};
            if strcmp(canonical_name, property_name)
                if isempty(value)
                    if isempty(varargin{2*i})
                        error(['Empty value for ' property_name '.']);
                    end
                    value = varargin{2*i};
                else
                    error(['Property ' property_name ' is specified more than once.']);
                end
            end
        end
    end
end

%--------------------------------------------------------------------------
function canonical_string = GetCanonicalString(...
               input_string, property_name, canonical_strings)

if ~ischar(input_string)
    error([property_name ' must be a string.']);
end

imatch = strmatch(lower(input_string),canonical_strings);
nmatch = length(imatch);

if nmatch == 0
    error(['Unrecognized ' property_name ' ' input_string '.']);
end

if nmatch > 1
    error(['Ambiguous ' property_name ' ' input_string '.']);
end

canonical_string = canonical_strings{imatch};

%--------------------------------------------------------------------------
function CheckPropertyNames( type, valid_property_names, varargin )

for i = 1:((nargin-2)/2)
    current_name = varargin{2*i-1};
    if ischar(current_name)
        nmatch = length(strmatch(lower(current_name),valid_property_names));
        if nmatch == 0
            warning(['Property name ' current_name ' ignored with type ' type '.']);
        end
        if nmatch > 1
            error(['Ambiguous property name ' current_name '.']);
        end
    else
        error('Non-character property name.');
    end
end

% ---------------------------------------------------------------------------
function tflip = fliptform( t )
%FLIPTFORM Flip the input and output roles of a TFORM structure.
%   TFLIP = FLIPTFORM(T) creates a new spatial transformation structure (a
%   "TFORM struct") by flipping the roles of the inputs and outputs in an
%   existing TFORM struct.
%
%   Example
%   -------
%       T = maketform('affine',[.5 0 0; .5 2 0; 0 0 1]);
%       T2 = fliptform(T);
%
%   The following are equivalent:
%       x = tformfwd([-3 7],T)
%       x = tforminv([-3 7],T2)

%   Copyright 1993-2003 The MathWorks, Inc.
%   $Revision: 1.9 $ $Date: 2003/01/17 16:27:29 $

checknargin(1,1,nargin,mfilename);

if ~istform(t) || (length(t) ~= 1)
    eid = sprintf('Images:%s:tMustBeSingleTformStruct',mfilename);
    error(eid, 'T must be a single TFORM struct.');
end

tflip = maketform('custom', t.ndims_out, t.ndims_in, t.inverse_fcn, t.forward_fcn, t.tdata);

% ---------------------------------------------------------------------------
function varargout = tforminv(varargin)
%TFORMINV Apply inverse spatial transformation.
%   TFORMINV applies an inverse spatial transformation based on a TFORM
%   structure created with MAKETFORM, FLIPTFORM, or CP2TFORM.
%   
%   [U,V] = TFORMINV(T,X,Y) applies the 2D-to-2D inverse transformation
%   defined in TFORM structure T to coordinate arrays X and Y, mapping
%   the point [X(k) Y(k)] to the point [U(k) V(k)].  Both T.ndims_in
%   and T.ndims_out must equal 2.  X and Y will typically be column
%   vectors matching in length.  In general, X and Y can have any
%   dimensionality, but must have the same size.  In any case, U and V
%   will have the same size as X and Y.
%
%   [U1,U2,U3,...] = TFORMINV(T,X1,X2,X3,...) applies the NDIMS_OUT-to-
%   NDIMS_IN inverse transformation defined in TFORM structure T to the
%   coordinate arrays X1,X2,...,XNDIMS_OUT (where NDIMS_IN = T.ndims_in
%   and NDIMS_OUT = T.ndims_out).  The number of output arguments
%   must equal NDIMS_IN.  The transformation maps the point
%              [X1(k) X2(k) ... XNDIMS_OUT(k)]
%   to the point
%              [U1(k) U2(k) ... UNDIMS_IN(k)].
%   X1,X2,X3,... can have any dimensionality, but must be the same size.
%   U1,U2,U3,... will have this size also.
%
%   U = TFORMINV(T,X) applies the NDIMS_OUT-to-NDIMS_IN inverse
%   transformation defined in TFORM structure T to each row of X, where
%   X is an M-by-NDIMS_OUT matrix.  It maps the point X(k,:) to the
%   point U(k,:).  U will be an M-by-NDIMS_IN matrix.
%
%   U = TFORMINV(T,X), where X is an (N+1)-dimensional array, maps
%   the point X(k1,k2,...,kN,:) to the point U(k1,k2,...,kN,:).
%   SIZE(X,N+1) must equal NDIMS_OUT.  U will be an (N+1)-dimensional
%   array, with SIZE(U,I) equal to SIZE(X,I) for I = 1,...,N and
%   SIZE(U,N+1) equal to NDIMS_IN.
%
%   [U1,U2,U3,...] = TFORMINV(T,X) maps an (N+1)-dimensional array
%   to NDIMS_IN equally-sized N-dimensional arrays.
%
%   U = TFORMINV(T,X1,X2,X3,...) maps NDIMS_OUT N-dimensional arrays
%   to one (N+1)-dimensional array.
%
%   Note
%   ----
%   U = TFORMINV(X,T) is an older form of the two-argument syntax
%   that remains supported for backward compatibility.
% 
%   Example
%   -------
%   Create an affine transformation that maps the triangle with vertices
%   (0,0), (6,3), (-2,5) to the triangle with vertices (-1,-1), (0,-10),
%   (4,4):
%
%       u = [ 0   6  -2]';
%       v = [ 0   3   5]';
%       x = [-1   0   4]';
%       y = [-1 -10   4]';
%       tform = maketform('affine',[u v],[x y]);
%   
%   Validate the mapping by applying TFORMINV:
%
%       [um, vm] = tforminv(tform, x, y)  % Results should equal [u, v]

%   Copyright 1993-2003 The MathWorks, Inc.
%   $Revision: 1.11 $ $Date: 2003/02/12 20:24:06 $

varargout = tform('inv', nargout, varargin{:});

% ---------------------------------------------------------------------------
function varargout = tformfwd(varargin)
%TFORMFWD Apply forward spatial transformation.
%   TFORMFWD applies a forward spatial transformation based on a TFORM
%   structure created with MAKETFORM, FLIPTFORM, or CP2TFORM.

varargout = tform('fwd', nargout, varargin{:});

% --------------------------------------------------------------------------
function outputs = tform(direction, numout, varargin)
% Perform either a forward or inverse spatial transformation,
% in support of TFORMFWD and TFORMINV.

%   Copyright 1993-2003 The MathWorks, Inc.
%   $Revision: 1.3 $  $Date: 2003/07/17 21:14:30 $

% Set up a structure to control which direction things go,
% and to provide strings to use in error calls.

if strcmp(direction,'fwd')
    f.name      = 'TFORMFWD';
    f.ndims_in  = 'ndims_in';
    f.ndims_out = 'ndims_out';
    f.fwd_fcn   = 'forward_fcn';
    f.argname   = 'U';
    f.arglist   = 'U1,U2,U3,...';
else
    f.name      = 'TFORMINV';
    f.ndims_in  = 'ndims_out';
    f.ndims_out = 'ndims_in';
    f.fwd_fcn   = 'inverse_fcn';
    f.argname   = 'X';
    f.arglist   = 'X1,X2,X3,...';
end   

checknargin(2,inf,numel(varargin),f.name);

% Get the TFORM struct, t, and the remaining input arguments, A.
[t, A] = checkTform(f, numout, varargin{:});

P = t.(f.ndims_in);   % Dimensionality of input space.
L = t.(f.ndims_out);  % Dimensionality of output space.

% Validate and organize the input coordinate data so that:
% U contains the input coordinates in an (N + 1)-dimensional array.
% D contains the sizes of the first N dimensions of U, after correctly
% accounting for trailing singleton dimensions that may have been
% dropped or added.

if length(A) == 1
    [U, D] = checkCoordinates( f, P, A{1} );
else
    [U, D] = concatenateCoordinates( f, P, A );
end

% Reshape U to collapse the first N dimensions into a single column.
% Each row of U is now a coordinate vector (of length NDIMS_IN) in
% t's input space.

U = reshape( U, [prod(D) P] );

% Apply the transformation, mapping U from the input space to the
% output space.

X = feval( t.(f.fwd_fcn), U, t );

% X is an array in which each row is a coordinate vector in t's output
% space. Reshape X to be D(1)-by-D(2)-by-...-by-D(N)-by-NDIMS_OUT.

X = reshape( X, [D L] );

% Construct output.
if (numout <= 1)
    outputs{1} = X;
else
    outputs = separateCoordinates( X, numel(D), L );
end

%----------------------------------------------------------------------
function [t, A] = checkTform(f, numout, varargin)
% Look at varargin{1} and varargin{2} to locate the TFORM structure.
% Return the remaining arguments (coordinate arrays) in cell array A.

if isstruct(varargin{1})
    t = varargin{1};
    A = varargin(2:end);
elseif isstruct(varargin{2})
    % Old syntax: TFORMFWD(U,T) or TFORMINV(X,T).
    if numel(varargin) == 2
        t = varargin{2};
        A = varargin(1);
    else
        error(sprintf('Images:%s:TooManyInputs',f.name),...
              'Function %s accepts only two arguments when T is the second argument (older syntax).',f.name);
    end
else
    error(sprintf('Images:%s:MissingTform',f.name),...
          'Function %s expected a TFORM struct as its first or second argument.',f.name);
end

if ~istform(t) || (numel(t) ~= 1)
    error(sprintf('Images:%s:InvalidTform',f.name),...
         'Function %s expected T to be a single TFORM struct.', f.name);
end

if isempty(t.(f.fwd_fcn))
    error(sprintf('Images:%s:TformMissingFcn',f.name),...,
         'Function %s expected T''s %s field to be non-empty.',f.name, f.fwd_fcn);
end

if (numout > 1) && (numout ~= t.(f.ndims_out))
    error(sprintf('Images:%s:OutputCountTformMismatch',f.name),...
          'Function %s expected the number of its output arguments to be consistent with T.%s.',...
          f.name, f.ndims_out);
end

%----------------------------------------------------------------------
function [U, D] = checkCoordinates(f, P, U)
% Let U have (N + 1) dimensions after allowing for dropped or extraneous
% trailing singletons.  Determine the sizes of the first N dimensions of
% U and return them in row vector D.  P is the dimensionality of the input space.

iptcheckinput(U, {'double'}, {'real','finite'}, f.name, 'U', 2);

M = ndims(U);
S = size(U);

if S(M) > 1
    if P == 1;
        N = M;         % PROD(S) points in 1-D
    elseif P == S(M)
        N = M - 1;     % PROD(S(1:M-1)) points in P-D
    else
        error(sprintf('Images:%s:ArraySizeTformMismatch1',f.name),...
             'Function %s: SIZE(%s) is inconsistent with T.%s.',...
              f.name, f.argname, f.ndims_in);
    end
    D = S(1:N);
else % S == [S(1) 1]
    if P == 1
        D = S(1);      % S(1) points in 1-D
    elseif P == S(1)
        D = 1;         % 1 point in P-D
    else
        error(sprintf('Images:%s:ArraySizeTformMismatch2',f.name),...
             'Function %s: SIZE(%s) is inconsistent with T.%s.',...
              f.name, f.argname, f.ndims_in);
    end        
end

%----------------------------------------------------------------------
function [U, D] = concatenateCoordinates(f, P, A)
% A is a cell array containing the coordinate arguments:
%
%                A = {U1, U2, U3, ...}
%
% Validate arguments U1,U2,U3,..., then concatenate them into U.
% If the Uk are column vectors, then U is a matrix (2D array) and
% D = length(Uk).  Otherwise, U is an (N + 1)-dimensional array,
% where N = ndims(Uk), including implicit trailing singletons.
% D is a row vector containing the sizes of the first N dimensions
% of U.  P is the dimensionality of the input space.

% Check argument count.
if length(A) ~= P
    error(sprintf('Images:%s:InputCountTformMismatch',f.name),...
          'Function %s expected either 2 or (1 + T.%s) input arguments.',f.name, f.ndims_in);
end

% Check argument class, properties, consistency.
ndims1 = ndims(A{1});
size1 = size(A{1});
for k = 1:P
    desc = sprintf('%s%d', f.argname, k);
    iptcheckinput(A{k}, {'double'}, {'real','finite'}, f.name, desc, k+1);
    if any(ndims(A{k}) ~= ndims1) || any(size(A{k}) ~= size1)
        error(sprintf('Images:%s:ArraySizeMismatch',f.name),...
              'Function %s expected the sizes of %s to match.', f.name, f.arglist);
    end
end

% Determine the size vector, D.
D = size1;
if (numel(D) == 2) && (D(2) == 1)
    % U1,U2,... are column vectors.  They must be specified as
    % 1-D because MATLAB does not support explicit 1-D arrays.    
    D = D(1);
end

% Concatenate the coordinate arguments.
N = numel(D);
U = cat(N + 1, A{:});

%----------------------------------------------------------------------
function outputs = separateCoordinates(X, N, L)
% Distribute X to output arguments along its (N+1)-th dimension, which has size L.
subs = repmat( {':'}, [1 N] );
for k = 1:L
    outputs{k} = X(subs{:},k);
end

%----------------------------------------------------------------------
function varargout = imrotate(varargin)
%IMROTATE Rotate image.
%   B = IMROTATE(A,ANGLE) rotates image A by ANGLE degrees in a 
%   counterclockwise direction around its center point. To rotate the image
%   clockwise, specify a negative value for ANGLE. IMROTATE makes the output
%   image B large enough to contain the entire rotated image. IMROTATE uses
%   nearest neighbor interpolation, setting the values of pixels in B that 
%   are outside the rotated image to 0 (zero).
%
%   B = IMROTATE(A,ANGLE,METHOD) rotates image A, using the interpolation
%   method specified by METHOD. METHOD is a string that can have one of the
%   following values. The default value is enclosed in braces ({}).
%
%        {'nearest'}  Nearest neighbor interpolation
%
%        'bilinear'   Bilinear interpolation
%
%        'bicubic'    Bicubic interpolation. Note: This interpolation method
%                     method can produce pixel values outside the original range
%
%   B = IMROTATE(A,ANGLE,METHOD,BBOX) rotates image A, where BBOX specifies 
%   the size of the output image B. BBOX is a text string that can have 
%   either of the following values. The default value is enclosed in braces
%   ({}).
%
%        {'loose'}    Make output image B large enough to contain the
%                     entire rotated image. B is generally larger than A.
%
%        'crop'       Make output image B the same size as the input image
%                     A, cropping the rotated image to fit. 
%
%   Class Support
%   -------------
%   The input image can be numeric or logical.  The output image is of the
%   same class as the input image.
%
%   Example
%   -------
%        % This example brings image I into horizontal alignment by
%        % rotating the image by -1 degree.
%        
%        I = fitsread('solarspectra.fts');
%        I = mat2gray(I);
%        J = imrotate(I,-1,'bilinear','crop');
%        figure, imshow(I), figure, imshow(J)
%
%   See also IMCROP, IMRESIZE, IMTRANSFORM, TFORMARRAY.

%   Copyright 1992-2004 The MathWorks, Inc.
%   $Revision: 5.25.4.7 $  $Date: 2004/12/18 07:35:58 $

[A,ang,method,bbox] = parse_inputsImrotate(varargin{:});

so = size(A);
twod_size = so(1:2);

if (rem(ang,90) == 0)       % Catch and speed up 90 degree rotations
    % determine if angle is +- 90 degrees or 0,180 degrees.
    multiple_of_ninety = mod(floor(ang/90), 4);
    v = repmat({':'},[1 ndims(A)]);       % initialize array of subscripts

    switch multiple_of_ninety 
        case 0                   % 0 rotation;
        B = A;
    case {1,3}               % +- 90 deg rotation 
        thirdD = prod(so(3:end));
        A = reshape(A,[twod_size thirdD]);
        not_square = twod_size(1) ~= twod_size(2);
        if (bbox(1) == 'c') && not_square       % center rotated image and preserve size
            imbegin = (max(twod_size) == so)*abs(diff(floor(twod_size/2)));
            vec = 1:min(twod_size);
            v(1) = {imbegin(1)+vec};
            v(2) = {imbegin(2)+vec};
            new_size = [twod_size thirdD];
        else       % don't preserve original size
            new_size = [fliplr(twod_size) thirdD];
        end
     
        % pre-allocate array
        if islogical(A),    B = false(new_size);
        else                B = alloc_mex(new_size,class(A));       % B = zeros(new_size,class(A));
        end
     
        for k = 1:thirdD
            B(v{1},v{2},k) = rot90(A(v{1},v{2},k),multiple_of_ninety);
        end
     
        B = reshape(B,[new_size(1) new_size(2) so(3:end)]);
    
    case 2    % 180 rotation
        v(1) = {twod_size(1):-1:1};
        v(2) = {twod_size(2):-1:1};
        B = A(v{:});
    end

else        % use tformarray
    phi = ang*pi/180; % Convert to radians
    rotate = maketform('affine',[ cos(phi)  sin(phi)  0; ...
                               -sin(phi)  cos(phi)  0; ...
                                    0       0       1 ]);
    % Coordinates from center of A
    hiA = (twod_size-1)/2;
    loA = -hiA;
    if bbox(1) == 'l'  % Determine limits for rotated image
        hiB = ceil(max(abs(tformfwd([loA(1) hiA(2); hiA(1) hiA(2)],rotate)))/2)*2;
        loB = -hiB;
        sn = hiB - loB + 1;
    else % Cropped image
        hiB = hiA;
        loB = loA;
        sn = twod_size;
    end
  
    boxA = maketform('box',twod_size,loA,hiA);
    boxB = maketform('box',sn,loB,hiB);
    T = maketform('composite',[fliptform(boxB),rotate,boxA]);
  
    if strcmp(method,'bicubic')
        R = makeresampler('cubic','fill');
    elseif strcmp(method,'bilinear')
        R = makeresampler('linear','fill');
    else
        R = makeresampler('nearest','fill');
    end

    F = 0;
    if (isa(A,'single') || isa(A,'double'))         % JL
        F = NaN;
    elseif (isa(A,'int32'))
        F = -2147483648;
    elseif (isa(A,'int16'))
        F = -32768;
    end
    B = tformarray(A, T, R, [1 2], [1 2], sn, [], F);
    %B = tformarray(A, T, R, [1 2], [1 2], sn, [], 0);
end
   
% Output
switch nargout,
	case 1,
        varargout{1} = B;
	case 3,
        wid = 'Images:imrotate:obsoleteSyntax';    
        warning(wid, '%s', ['[R,G,B] = IMROTATE(RGB) is an obsolete output syntax. ',...
        'Use one output argument the receive the 3-D output RGB image.']);
        for k=1:3,
            varargout{k} = B(:,:,k);
        end
	otherwise
        error('Images:imrotate:tooManyOutputs', '%s', 'Invalid number of output arguments.');
end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Function: parse_inputs
function [A,ang,method,bbox] = parse_inputsImrotate(varargin)
% Outputs:  A       the input image
%           ang     the angle by which to rotate the input image
%           method  interpolation method (nearest,bilinear,bicubic)
%           bbox    bounding box option 'loose' or 'crop'

% Defaults:
method = 'n';   bbox = 'l';

error(nargchk(2,4,nargin));
switch nargin
	case 2,             % imrotate(A,ang)        
        A = varargin{1};
        ang=varargin{2};
	case 3,             % imrotate(A,ang,method) or
        A = varargin{1};  % imrotate(A,ang,box)
        ang=varargin{2};
        method=varargin{3};
	case 4,             % imrotate(A,ang,method,box) 
        A = varargin{1};
        ang=varargin{2};
        method=varargin{3};
        bbox=varargin{4};
	otherwise,
        error('Images:imrotate:invalidInputs', '%s', 'Invalid input arguments.');
end

% Check validity of the input parameters 
if ischar(method) && ischar(bbox),
  strings = {'nearest','bilinear','bicubic','crop','loose'};
  idx = strmatch(lower(method),strings);
  if isempty(idx),
    eid = 'Images:imrotate:unrecognizedInterpolationMethod';
    error(eid, 'Unknown interpolation method: %s', method);
  elseif length(idx)>1,
    eid = 'Images:imrotate:ambiguousInterpolationMethod';
    error(eid, 'Ambiguous interpolation method: %s', method);
  else
    if idx==4,bbox=strings{4};method=strings{1};
    elseif idx==5,bbox = strings{5};method=strings{1};
    else method = strings{idx};
    end
  end  
  idx = strmatch(lower(bbox),strings(4:5));
  if isempty(idx),
    error('Images:imrotate:unrecognizedBBox', 'Unknown BBOX parameter: %s', bbox);
  elseif length(idx)>1,
    error('Images:imrotate:ambiguousBBox', 'Ambiguous BBOX string: %s', bbox);
  else
    bbox = strings{3+idx};
  end 
else
  eid = 'Images:imrotate:expectedString';
  error(eid, '%s', 'Interpolation method and BBOX have to be a string.');  
end

% ----------------------------------------------------------------------------------
function outbounds = findbounds(varargin)
%FINDBOUNDS Find output bounds for spatial transformation.
%   OUTBOUNDS = FINDBOUNDS(TFORM,INBOUNDS) estimates the output bounds
%   corresponding to a given spatial transformation and a set of input
%   bounds.  TFORM is a spatial transformation structure as returned by
%   MAKETFORM or CP2TFORM.  INBOUNDS is 2-by-NUM_DIMS matrix.  The first row
%   of INBOUNDS specifies the lower bounds for each dimension, and the
%   second row specifies the upper bounds. NUM_DIMS has to be consistent
%   with the ndims_in field of TFORM.
%
%   OUTBOUNDS has the same form as INBOUNDS.  It is an estimate of the
%   smallest rectangular region completely containing the transformed
%   rectangle represented by the input bounds.  Since OUTBOUNDS is only an
%   estimate, it may not completely contain the transformed input rectangle.
%
%   Notes
%   -----
%   IMTRANSFORM uses FINDBOUNDS to compute the 'OutputBounds' parameter
%   if the user does not provide it.
%
%   If TFORM contains a forward transformation (a nonempty forward_fcn
%   field), then FINDBOUNDS works by transforming the vertices of the input
%   bounds rectangle and then taking minimum and maximum values of the
%   result.
%
%   If TFORM does not contain a forward transformation, then FINDBOUNDS
%   estimates the output bounds using the Nelder-Mead optimization
%   function FMINSEARCH.  If the optimization procedure fails, FINDBOUNDS
%   issues a warning and returns OUTBOUNDS=INBOUNDS.
%
%   Example
%   -------
%       inbounds = [0 0; 1 1]
%       tform = maketform('affine',[2 0 0; .5 3 0; 0 0 1])
%       outbounds = findbounds(tform, inbounds)
%
%   See also CP2TFORM, IMTRANSFORM, MAKETFORM, TFORMARRAY, TFORMFWD, TFORMINV.

%   Copyright 1993-2003 The MathWorks, Inc.
%   $Revision: 1.7 $  $Date: 2003/01/17 16:27:29 $

% I/O details
% -----------
% tform     - valid TFORM structure; checked using private/istform.
%
% inbounds  - 2-by-NUM_DIMS real double matrix.  NUM_DIMS must be equal to
%             tform.ndims_in.  It may not contain NaN's or Inf's.
%
% outbounds - 2-by-NUM_DIMS_OUT real double matrix.  NUM_DIMS_OUT is
%             equal to tform.ndims_out.

[tform,inbounds] = parse_inputs_findbounds(varargin{:});

if isempty(tform.forward_fcn)
    outbounds = find_bounds_using_search(tform, inbounds);
else
    outbounds = find_bounds_using_forward_fcn(tform, inbounds);
end

%--------------------------------------------------
function out_bounds = find_bounds_using_forward_fcn(tform, in_bounds)

in_vertices = bounds_to_vertices(in_bounds);
in_points = add_in_between_points(in_vertices);
out_points = tformfwd(in_points, tform);
out_bounds = points_to_bounds(out_points);

%--------------------------------------------------
function out_bounds = find_bounds_using_search(tform, in_bounds)

% Strategy
% --------
% For each point u_k in a set of points on the boundary or inside of the
% input bounds, find the corresponding output location by minimizing this
% objective function:
%
%    norm(u_k - tforminv(x, tform))
%
% It seems reasonable to use the u_k values as starting points for the
% optimization routine, FMINSEARCH.

if isempty(tform.inverse_fcn)
    msg = sprintf('%s: forward_fcn and inverse_fcn fields of TFORM cannot both be empty.',...
          upper(mfilename));
    error(sprintf('Images:%s:fwdAndInvFieldsCantBothBeEmpty',mfilename),msg);
end

in_vertices = bounds_to_vertices(in_bounds);
in_points = add_in_between_points(in_vertices);
out_points = zeros(size(in_points));
success = 1;
options = optimset('Display','off');
for k = 1:size(in_points,1)
    [x,fval,exitflag] = fminsearch(@objective_function, in_points(k,:), ...
                                   options, tform, in_points(k,:));
    if exitflag <= 0
        success = 0;
        break;
    else
        out_points(k,:) = x;
    end
end

if success
    out_bounds = points_to_bounds(out_points);
else
    % Optimization failed; the fallback strategy is to make the output
    % bounds the same as the input bounds.  However, if the input
    % transform dimensionality is not the same as the output transform
    % dimensionality, there doesn't seem to be anything reasonable to do.
    if tform.ndims_in == tform.ndims_out
        msg = sprintf('%s: Search procedure failed; returning OUTBOUNDS = INBOUNDS.',upper(mfilename));
        warning(sprintf('Images:%s:searchFailed',mfilename),msg);
        out_bounds = in_bounds;
    else
        msg = sprintf('%s: Search procedure failed with a mixed-dimensionality TFORM.',upper(mfilename));
        error(sprintf('Images:%s:mixedDimensionalityTFORM',mfilename),msg);
    end
end

%--------------------------------------------------
function s = objective_function(x, tform, u0)
% This is the function to be minimized by FMINSEARCH.
s = norm(u0 - tforminv(x, tform));

%--------------------------------------------------
function vertices = bounds_to_vertices(bounds)
% Convert a 2-by-num_dims bounds matrix to a 2^num_dims-by-num_dims
% matrix containing each of the vertices of the region corresponding to
% BOUNDS.
%
% Strategy: the k-th coordinate of each vertex bound can be either
% bounds(k,1) or bounds(k,2).  One way to enumerate all the possibilities
% is to count in binary from 0 to (2^num_dims - 1).

num_dims = size(bounds,2);
num_vertices = 2^num_dims;

binary = repmat('0',[num_vertices,num_dims]);
for k = 1:num_vertices
    binary(k,:) = dec2bin(k-1,num_dims);
end

mask = binary ~= '0';

low = repmat(bounds(1,:),[num_vertices 1]);
high = repmat(bounds(2,:),[num_vertices 1]);
vertices = low;
vertices(mask) = high(mask);

%--------------------------------------------------
function points = add_in_between_points(vertices)
% POINTS contains all of the input vertices, plus all the unique points
% that are in between each pair of vertices.

[num_vertices,num_dims] = size(vertices);
ndx = nchoosek(1:num_vertices,2);
new_points = (vertices(ndx(:,1),:) + vertices(ndx(:,2),:))/2;
new_points = unique(new_points, 'rows');
points = [vertices; new_points];

%--------------------------------------------------
function bounds = points_to_bounds(points)
% Find a 2-by-num_dims matrix bounding the set of points in POINTS.
bounds = [min(points,[],1) ; max(points,[],1)];

%--------------------------------------------------
function [tform,inbounds] = parse_inputs_findbounds(varargin)

checknargin(2,2,nargin,mfilename)

tform = varargin{1};
inbounds = varargin{2};

if ~istform(tform)
    msg = sprintf('%s: First input argument must be a TFORM struct.',upper(mfilename));
    error(sprintf('Images:%s:firstInputMustBeTformStruct',mfilename),msg);
end

if numel(tform) ~= 1
    msg = sprintf('%s: First input argument must be a 1-by-1 TFORM struct.',upper(mfilename));
    error(sprintf('Images:%s:firstInputMustBeOneByOneTformStruct',mfilename),msg);
end

if tform.ndims_in ~= tform.ndims_out
    msg = sprintf('%s: Input and output dimensions of TFORM must be the same.',upper(mfilename));
    error(sprintf('Images:%s:inOutDimsOfTformMustBeSame',mfilename),msg);
end

if ~isnumeric(inbounds) || (ndims(inbounds) > 2) || (size(inbounds,1) ~= 2)
    msg = sprintf('%s: INBOUNDS must be a 2-by-NUM_DIMS numeric matrix.',upper(mfilename));
    error(sprintf('Images:%s:inboundsMustBe2byN',mfilename),msg);
end

num_dims = size(inbounds,2);

if num_dims ~= tform.ndims_in
    msg = sprintf('%s: size(INBOUNDS,2) must equal TFORM.ndims_in.',upper(mfilename));
    error(sprintf('Images:%s:secondDimOfInbundsMustEqualTformNdimsIn',mfilename),msg);
end

if any(~isfinite(inbounds(:)))
    msg = sprintf('%s: INBOUNDS must contain only finite values.', upper(mfilename));
    error(sprintf('Images:%s:inboundsMustBeFinite',mfilename),msg);
end

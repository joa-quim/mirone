function [nxout,nyout,nzout] = getnormals(x,y,z,s)
%   [Nx,Ny,Nz] = GETNORMALS(X,Y,Z) returns the components of the 3-D surface normal
%   for the surface with components (X,Y,Z). The normal is normalized to length 1.
%
%   s = GETNORMALS(X,Y,Z,[azim elev]) returns the components of the 3-D surface normal
%   projected along the vector S. S = [azim,elev]. This form is used in one of the
%   illumination algorithm, and the fact that is computed here allows saving a lot of space
%
%   Nz = GETNORMALS(X,Y,Z) returns only the vertical component of the normal vector
%	This saves a little space since the horizontal components are not computed
%
%   The surface normals returned are based on a bicubic fit of the data.

%   This function is based on SURFNORM but, contrary to it that takes about 14 times
%   the grid size in memory, this one "only" takes 7 times. A number of useless operations
%   where also removed, resulting in a much faster routine.
%
%	Coffeeright (c) 2004-2012 by J. Luis

[m,n] = size(z);
if any([m n]<3), error('Z must be at least 3-by-3.'); end
% If we get X & Y meshgrided, get rid of that stupidity
if (size(x,1) > 1)
    x = x(1,:);
    y = y(:,1);
else                % Ensure that X is a row vector and Y a column vector
    x = x(:)';      y = y(:);
end
if ~isequal(length(y),m), error('Y must have the same number of rows as Z.'); end
if ~isequal(length(x),n), error('X must have the same number of columns as Z.'); end

stencil1 = [1 0 -1]/2;    stencil2 =  [-1; 0; 1]/2;

% Expand x,y,z so interpolation is valid at the boundaries.
x = [3*(x(1,1)-x(1,2))+x(1,3),x,3*(x(1,n)-x(1,n-1))+x(1,n-2)];
y = [3*(y(1,1)-y(2,1))+y(3,1);y;3*(y(m,1)-y(m-1,1))+y(m-2,1)];
z = [3*(z(1,:)-z(2,:))+z(3,:);z;3*(z(m,:)-z(m-1,:))+z(m-2,:)];
z = [3*(z(:,1)-z(:,2))+z(:,3),z,3*(z(:,n)-z(:,n-1))+z(:,n-2)];

rows = 2:m+1; cols = 2:n+1;
%ay = ffilter(stencil1,y);      % This component is allways 0
tmp = ffilter(stencil1,x(1,:));     ax = repmat(tmp(cols),m,1);
az = ffilter(stencil1,z);           az = az(rows,cols);

%bx = ffilter(stencil2,x);      % This component is allways 0
tmp = ffilter(stencil2,y(:,1));     by = repmat(tmp(rows),1,n);
bz = ffilter(stencil2,z);           bz = bz(rows,cols);

% Perform cross product to get normals
% nx = -(ay.*bz - az.*by);
% ny = -(az.*bx - ax.*bz);
% nz = -(ax.*by - ay.*bx);
nx = az.*by;    clear az;       % by = cte (for grids)
ny = ax.*bz;    clear bz;       % ax = cte (for grids)
nz = -ax.*by;   clear ax by;    % nz = cte (for grids)

% Normalize the length of the surface normals to 1.
mag = sqrt(nx.*nx + ny.*ny + nz.*nz);
d = find(mag==0); mag(d) = eps*ones(size(d));

if (nargin == 4 && numel(s) == 2)
	D2R = pi / 180;
	azim = s(1);    elev = s(2);
	theta = elev*D2R; phi = azim*D2R; 
	St = sin(theta); Ct = cos(theta); Sp = sin(phi); Cp = cos(phi);
	nxout = (nx*Ct*Sp + ny*Ct*Cp + nz*St) ./ mag;
	if (nargout == 3),  nyout = [];    nzout = [];   end   % Just a precaution
else
	if (nargout == 1)
		nxout = nz ./mag;
	else
		nxout = nx ./mag;   clear nx;
		nyout = ny ./mag;   clear ny;
		nzout = nz ./mag;
	end
end

% --------------------------------------------------------------------
function y = ffilter(b,x)
%see FILTER2 for help (I did this for striping unnecessry code and reducing the .exe version size).

	stencil = rot90(b,2);
	[ms,ns] = size(stencil);

	% 1-D filter
	if (ms == 1)
	  y = conv2(1,stencil,x,'same');
	elseif (ns == 1)
	  y = conv2(stencil,1,x,'same');
	else
		error('The filter must be 1-D')
	end

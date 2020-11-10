function [slope,aspect,gradN,gradE] = gradient_geo(varargin)

% [slope,aspect,gradN,gradE] = gradient_geo(lat,lon,Z) computes the slope, 
% aspect and the components of the gradient for the matrix Z.
%
% Limiting what goes out is also possible.
% slope = gradient_geo(lat,lon,Z,'slope') returns only the slope
% aspect = gradient_geo(lat,lon,Z,'aspect') returns only the aspect
% gradN = gradient_geo(lat,lon,Z,'gradN') returns only the gradient due North
% gradE = gradient_geo(lat,lon,Z,'gradE') returns only the gradient due East
% [gradN gradE] = gradient_geo(lat,lon,Z,'grad') returns only the gradients due North and East
% [...] = gradient_geo(lat,lon,Z,...,'cart') Do as above but Z is taken as a cartesian array 
% [slope,aspect,gradN,gradE] = gradient_geo(lat,lon,Z,'cart') Compute all 4 on a cartesian array
%
%  11-10-03     Updated to do Tiling. Now it hopfully can be used for any grid size
%  03-01-09		Operates also on cartesian arrays

if (nargin < 3 || nargin > 5)
	error('Incorrect number of input arguments.')
end
only_slope = 0;     only_aspect = 0;    only_gradN = 0;     only_gradE = 0;     only_grad = 0;  todos = 0;
geog = true;

% Test  to see if varargin{1} & varargin{2} are vectors or meshgrids
if ( min(size(varargin{1})) == 1 && min(size(varargin{2})) == 1 && ...    % lat,lon vectors as input
        isequal(numel(varargin{2}), size(varargin{3},2)) && isequal(numel(varargin{1}), size(varargin{3},1)) )
    do_mesh = 1;                % Do it later inside the Tiling
    map = varargin{3};
elseif isequal( size(varargin{1}), size(varargin{2}), size(varargin{3}) )
    do_mesh = 0;                % Input already "meshgrided"
else
    errordlg('Eror in map input arguments to gradient_geo','Error');   return
end

[m,n] = size(varargin{3});
if (n < 3 || m < 3)
	errordlg('Matrix must be 4 by 4 or larger to compute gradient','Error'),	return
end

if (nargin == 3)
    todos = true;
elseif (nargin >= 4)
    switch varargin{4}
        case 'slope'
            if (nargout ~= 1),  errordlg('GRADIENT_GEO: Must have only one output (slope).','Error');   return; end
            only_slope = 1;
        case 'aspect'
            if (nargout ~= 1),  errordlg('GRADIENT_GEO: Must have only one output (aspect).','Error');   return; end
            only_aspect = 1;
        case 'gradN'
            if (nargout ~= 1),  errordlg('GRADIENT_GEO: Must have only one output (gradN).','Error');   return; end
            only_gradN = 1;
        case 'gradE'
            if (nargout ~= 1),  errordlg('GRADIENT_GEO: Must have only one output (GradE).','Error');   return; end
            only_gradE = 1;
        case 'grad'
            if (nargout ~= 2),  errordlg('GRADIENT_GEO: Must have two outputs (gradN & gradE).','Error');   return; end
            only_grad = 1;
		case 'cart'
			todos = true;
        otherwise
            errordlg(['GRADIENT_GEO: wrong code ' '"' varargin{4} '"'],'Error'),	return
    end
end
if (nargin == 5 && strncmp(varargin{5},'cart',4)),		geog = false;	end

D2R = pi/180;
geoid = [6378137 0.0818191910428158];       % grs80 ellipsoid

% convert from the geodetic latitude to the rectifying latitude.
% n = (semimajor axis - semiminor axis)/(semimajor axis + semiminor axis)
n = 0.00167922039463;       % Value for grs80 ellipsoide used here
f1 = 3*n / 2 - 9*n^3 / 16;  f2 = 15*n^2 / 16 - 15*n^4 / 32;
f3 = 35*n^3 / 48;           f4 = 315*n^4 / 512;

% Do the tiling
[ind_s,ind] = tile(m,300,4);
if size(ind_s,1) > 1
	gradE = zeros(size(varargin{3}));
	gradN = zeros(size(varargin{3}));
	old_height = 0;
	for i = 1:size(ind_s,1)
		tmp1 = (ind_s(i,1):ind_s(i,2));     % Indexes with overlapping zone
		tmp2 = ind(i,1):ind(i,2);           % Indexes of chunks without the overlaping zone
		if (do_mesh),	[lonmesh,latmesh] = meshgrid(varargin{2},varargin{1}(tmp1));	% varargin{1} contains lat
		else			[latmesh,lonmesh,map] = deal(varargin{1:3});					% ??? tenho de adaptar este tb.
		end
		if (geog)
			% Compute the gradient for the cell spacing in the projected cylindrical coordinates
			latmesh = latmesh * D2R;    lonmesh = lonmesh * D2R;
			latmesh = latmesh - f1*sin(2*latmesh) + f2*sin(4*latmesh) - f3*sin(6*latmesh) + f4*sin(8*latmesh);
	        [tmp_gradE,tmp_gradN] = mtxgradient(map(tmp1,:), geoid(1) * lonmesh, geoid(1) * latmesh, geog);
		else
	        [tmp_gradE,tmp_gradN] = mtxgradient(map(tmp1,:), lonmesh, latmesh, geog);
		end
        % Adjust the longitude gradient for the convergence of the meridians
		if (geog)
			convfactor = departure(zeros(size(latmesh)), ones(size(latmesh)),...
				latmesh,geoid) / departure(0,1,0,geoid);
			convfactor(convfactor == 0) = NaN;		% avoid divisions by zero     
			tmp_gradE = tmp_gradE ./ convfactor;    clear convfactor;
		end
		gradE(old_height+1:old_height+numel(tmp2), :) = tmp_gradE(tmp2,:);
		gradN(old_height+1:old_height+numel(tmp2), :) = tmp_gradN(tmp2,:);
		old_height = old_height + numel(tmp2);
	end
else
	if (do_mesh),		[lonmesh,latmesh] = meshgrid(varargin{2},varargin{1});		% varargin{1} contains lat
	else				[latmesh,lonmesh,map] = deal(varargin{1:3});				% ??? tenho de adaptar este tb.
	end
	if (geog)
		% Compute the gradient for the cell spacing in the projected cylindrical coordinates
		latmesh = latmesh * D2R;    lonmesh = lonmesh * D2R;
		latmesh = latmesh - f1*sin(2*latmesh) + f2*sin(4*latmesh) - f3*sin(6*latmesh) + f4*sin(8*latmesh);
		[gradE,gradN] = mtxgradient(map, geoid(1) * lonmesh, geoid(1) * latmesh, geog);
	else
		[gradE,gradN] = mtxgradient(map,lonmesh, latmesh, geog);
	end
	if (geog)
        % Adjust the longitude gradient for the convergence of the meridians
        convfactor = departure(zeros(size(latmesh)), ones(size(latmesh)),...
			latmesh,geoid) / departure(0,1,0,geoid);
        convfactor(convfactor == 0) = NaN;			% avoid divisions by zero
        gradE = gradE ./ convfactor;
	end
end
clear convfactor latmesh lonmesh;

% Now, compute only the strictly necessary
if (only_aspect)
    aspect  = cart2pol(gradE,gradN);
    aspect(gradN == 0 & gradE == 0) = NaN;
    aspect  = zero_to_2pi(-aspect-pi/2);    % convert back to degrees
    aspect  = aspect / D2R;
elseif (only_slope)
    [aspect,mag]  = cart2pol(gradE,gradN);      clear aspect;
    slope = atan(mag) / D2R;
elseif (todos)
    [aspect,mag] = cart2pol(gradE,gradN);
    aspect(gradN == 0 & gradE == 0) = NaN;
    aspect  = zero_to_2pi(-aspect-pi/2) / D2R;    % convert back to degrees
    slope = atan(mag) / D2R;
end

if (only_slope == 1)    % nothing to swap
elseif (only_aspect)    % swap slope with aspect
    slope = aspect;
elseif (only_gradN)     % swap slope with gradN
    slope = gradN;
elseif (only_gradE)     % swap slope with gradE
    slope = gradE;
elseif (only_grad)      % swap slope and aspect with gradN & gradE
    slope = gradN;  aspect = gradE;
end                     % output all four quantities

%-----------------------------------------------------------------------------------
function [dfdx,dfdy] = mtxgradient(f,x,y, geog)
%gradient for matrix with variable cell spacing
% New thing - If ~GEOG operate on an cartesian array

% Derivatives of function with respect to rows and columns
dfdc = zeros(size(f));		dfdr = zeros(size(f));

% Take forward differences on left and right edges
dfdr(1,:) = f(2,:) - f(1,:);		dfdr(end,:) = f(end,:) - f(end-1,:);
dfdc(:,1) = f(:,2) - f(:,1);		dfdc(:,end) = f(:,end) - f(:,end-1);

% Take centered differences on interior points
dfdr(2:end-1,:) = (f(3:end,:)-f(1:end-2,:)) / 2;
dfdc(:,2:end-1) = (f(:,3:end)-f(:,1:end-2)) / 2;

if (geog)
	% Differences of x and y with respect to row and column numbers
	dxdr = zeros(size(x));      dxdc = zeros(size(x));
	dxdr(1,:) = x(2,:)-x(1,:);  dxdr(end,:) = x(end,:)-x(end-1,:);  % Take forward differences on left and right edges
	dxdc(:,1) = x(:,2)-x(:,1);  dxdc(:,end) = x(:,end)-x(:,end-1);
	dxdr(2:end-1,:) = (x(3:end,:)-x(1:end-2,:))/2;                  % Take centered differences on interior points
	dxdc(:,2:end-1) = (x(:,3:end)-x(:,1:end-2))/2;

	dydr = zeros(size(y));      dydc = zeros(size(y));
	dydr(1,:) = y(2,:)-y(1,:);  dydr(end,:) = y(end,:)-y(end-1,:);  % Take forward differences on left and right edges
	dydc(:,1) = y(:,2)-y(:,1);  dydc(:,end) = y(:,end)-y(:,end-1);
	dydr(2:end-1,:) = (y(3:end,:)-y(1:end-2,:))/2;                  % Take centered differences on interior points
	dydc(:,2:end-1) = (y(:,3:end)-y(:,1:end-2))/2;

	colang = atan2(dydc,dxdc);				% Angles of mesh columns
	coldist = sqrt(dxdc.^2 + dydc.^2);		% distances between elements along mesh columns
	dfdc = dfdc ./ coldist;
	clear dxdc dydc coldist;				% free memory

	rowang = atan2(dydr,dxdr);				% Angles of mesh rows
	rowdist = sqrt(dxdr.^2 + dydr.^2);		% distances between elements along mesh rows
	dfdr = dfdr ./ rowdist;
	clear dxdr dydr rowdist;				% free memory
end

% derivatives in the x and y directions
if (geog)
	dfdx = dfdc .* cos(colang) + dfdr .* cos(rowang);
	dfdy = dfdr .* sin(rowang) + dfdc .* sin(colang);
else
	rowdist = y(2) - y(1);
	coldist = x(1,2) - x(1,1);
	dfdx = dfdc / coldist;
	dfdy = dfdr / rowdist;
end

%-----------------------------------------------------------------------------------
function dist = departure(lon1,lon2,lat,geoid)
D2R = pi/180;
% Ensure that longitudes are in the [0 2pi] range since they will be treated as distances.
lon1 = lon1 * D2R;      lon1 = zero_to_2pi(lon1);
lon2 = lon2 * D2R;      lon2 = zero_to_2pi(lon2);
r = parallel_rad(geoid,lat);  dist = r .* abs(lon1-lon2);

%-----------------------------------------------------------------------------------
function angout = zero_to_2pi(angin)
angout = pi*((abs(angin)/pi) - 2*ceil(((abs(angin)/pi)-1)/2)) .* sign(angin);
epsilon = -1e-8;    indx = find(angout<epsilon);

%  Shift the points in the -pi to 0 range to the pi to 2pi range
if ~isempty(indx);  angout(indx) = angout(indx) + 2*pi;  end;
indx = find(angout < 0);		%  Reset near zero points
if ~isempty(indx);  angout(indx) = zeros(size(indx));  end

%-----------------------------------------------------------------------------------
function r = parallel_rad(geoid,lat)
%  r = PARALLEL_RAD(geoid,lat) computes the parallel radius of curvature for the ellipsoid.
semimajor = geoid(1);		eccent = geoid(2);
% Compute the distance from the center of the geoid to the specified point
num = 1-eccent^2;			den = 1 - (eccent * cos(lat)).^2;
rho = semimajor * sqrt(num ./ den);  r = rho .* cos(lat);

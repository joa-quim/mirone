function [xout,yout] = spline_interp(x,y,BC)
% SPLINE_INTERP - Interpolate data with cardinal spline
%   SPLINE_INTERP(P,BC), where P is nx3 matrix of coordinates and BC is a
%   2x3 array containing the tangencies at the beginning and end of the curve
%
% Author unknown. I found this function inside GSPLINE that resides in FEX.
% Made several changes to make it more (much more) memory friend/efficient
% namely the sparse array calculations.
% J. Luis

NDIM = 2;			% Currently this function is accepting only 2D points. It should
					% not be dificult to change to 3D (to do whenever it will be needed)
RES  = 0.1;			% Interpolatiom resolution. 0.1 means that we'll interpolate
					% 10 intervals between each input data point
% P = [x(:) y(:) zeros(numel(x),1)];	% If we had 3D points
P = [x(:) y(:)];
n = size(P,1);		% Determine number of points

if (n <= 2)
	xout = x;	yout = y;   
	return
end

m = n - 2;			% Determine number of tangencies to account for

% Our known vector Reduces to -3*(Pn - Pn+2)
% Pa = P(1:m,:);	Pb = P(3:n,:);
PP = -3*(P(1:m,:) - P(3:n,:)); %;-[0.5*(1-t)*(1+b)*(1-c)*(Pa-Pb)];

% Create Try diagonal matrix (1,4,1)
% TM = (diag(4*ones(n,1)) + diag(ones(n-1,1),1) + diag(ones(n-1,1),-1));
e = ones(n,1);
TM = local_spdiags([e 4*e e], -1:1, n, n);		clear e

% Check for boundary conditions
if (nargin == 2)		% If no BC specified, impose zero curvature at endpoints
	PP = [6*(P(2,:) - P(1,:)); PP; 6*(P(n,:) - P(n-1,:))];
	TM(1,1:3) = [4 2 0];
	TM(n,n-2:n) = [0 2 4];
else					% Use BC as tangencies at endpoints
	PP = [BC(1,:); PP; BC(2,:)];
	TM(1,1:2) = [1 0];
	TM(n,n-1:n) = [0 1];
end

% Solve for uknown tangencies 
P_dot = TM\PP;
clear TM PP

% Set up matricies for solving
CMx = [P(1:n-1,1) P(2:n,1) P_dot(1:n-1,1) P_dot(2:n,1)];
CMy = [P(1:n-1,2) P(2:n,2) P_dot(1:n-1,2) P_dot(2:n,2)];
% CMz = [P(1:n-1,3) P(2:n,3) P_dot(1:n-1,3) P_dot(2:n,3)];		% If 3D points
clear P P_dot

% Sort our data points so that it makes sense for plotting
uu = 0:RES:1;			% To calculate interpolated points
MM = zeros(numel(uu),NDIM,n-1);
j = 1;
for u = uu
	% Calculate the Hermite Basis Functions
	H_u = [(1-3*u^2) + 2*u^3; 3*u^2 - 2*u^3; u - 2*u^2 + u^3; u^3 - u^2;];
% 	M = [CMx * H_u CMy * H_u CMz * H_u];		% If we had 3D points
	M = [CMx * H_u CMy * H_u];
	for q = 1:n-1
		MM(j,:,q) = M(q,:);        
	end
	j = j+1;
end
%for (i = 1:n-1),	out = [out; MM(:,:,i)];		end
% Cat along the 3rth dimension. E.G. MxNxK = K*MxN
out = reshape(permute(MM,[1 3 2]),(j-1)*(n-1),NDIM);

% Remove repeated points
out(diff(sum(out.^2,2)) == 0, :) = [];

xout = out(:,1);		yout = out(:,2);

if (size(x,1) == 1)		% Return result as row vectors, like the input
	xout = xout';		yout = yout';
end

% --------------------------------------------------------------------------
function res = local_spdiags(arg1,arg2,arg3,arg4)
% SPDIAGS Sparse matrix formed from diagonals.
% This is a minimalist excerpt of spdiags.m. Just the enough to solve the issue in question
	B = arg1;
	d = arg2(:);
	p = length(d);
	A = sparse(arg3,arg4);
	% Process A in compact form
	a = [];		i = [];
	[m,n] = size(A);
	for k = 1:p
		% Append new d(k)-th diagonal to compact form
		i = (max(1,1-d(k)):min(m,n-d(k)))';
		a = [a; i i+d(k) B(i+(m>=n)*d(k),k)];
	end
	res = sparse(a(:,1),a(:,2),full(a(:,3)),m,n);

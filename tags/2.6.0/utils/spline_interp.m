function [xout,yout] = spline_interp(x,y,varargin)
% SPLINE_INTERP - Interpolate data with cardinal spline
%   SPLINE_INTERP(X,Y), reinterpolate the X,Y vectors at 10 subdivisions between consecutive points
%   SPLINE_INTERP(X,Y,'Z',Z), does the same as above but along a 3D path
%   SPLINE_INTERP(...,'N',N), markes N subdivisions between consecutive data points.
%		For example, if N = 2 the interval between the data pairs X(k),Y(k) and 
%		X(k+1),Y(k+1) is split into N = 2 intervals by adding one data point in between them.
%		Default is N = 10.
%   SPLINE_INTERP(...,'BC',BC), where BC is a 2x3 array containing the tangencies
%		at the beginning and end of the curve
%
% Author unknown. I found this function inside GSPLINE that resides in FEX.
% Made several changes to make it more (much more) memory friend/efficient
% namely the sparse array calculations.
% J. Luis

% $Id$

	if (nargin < 2),		error('spline_interp: wrong number of arguments'),	end
	P = [x(:) y(:)];	NDIM = 2;	N = 10;		BC = [];
	if (nargin > 2)
		for (k = 1:2:numel(varargin))
			switch lower(varargin{k})
				case 'n'
					N = round(varargin{k+1});
					if (N <= 1)
						xout = x;	yout = y;  	return		% $€!?#
					end
				case 'z'
					z = varargin{k+1};
					P = [P z(:)];		NDIM = 3;
				case 'bc'
					BC = varargin{k+1};
			end
		end
	end

	% P = [x(:) y(:) zeros(numel(x),1)];	% If we had 3D points
	% P = [x(:) y(:)];
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
	if (isempty(BC))		% If no BC specified, impose zero curvature at endpoints
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
	if (NDIM == 3)
		CMz = [P(1:n-1,3) P(2:n,3) P_dot(1:n-1,3) P_dot(2:n,3)];
	end
	clear P P_dot

	% Sort our data points so that it makes sense for plotting
	uu = 0:1/N:1;			% To calculate interpolated points
	MM = zeros(numel(uu),NDIM,n-1);
	j = 1;
	for u = uu
		% Calculate the Hermite Basis Functions
		H_u = [(1-3*u^2) + 2*u^3; 3*u^2 - 2*u^3; u - 2*u^2 + u^3; u^3 - u^2;];
		if (NDIM == 2)
			M = [CMx * H_u CMy * H_u];
		else
			M = [CMx * H_u CMy * H_u CMz * H_u];
		end
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
	a = [];
	[m,n] = size(A);
	for k = 1:p
		% Append new d(k)-th diagonal to compact form
		i = (max(1,1-d(k)):min(m,n-d(k)))';
		a = [a; i i+d(k) B(i+(m>=n)*d(k),k)];
	end
	res = sparse(a(:,1),a(:,2),full(a(:,3)),m,n);

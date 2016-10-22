function [F, row, col] = bi_linear(arg1,arg2,arg3,arg4,arg5)
%LINEAR 2-D bilinear data interpolation.
%   ZI = BI_LINEAR(X,Y,Z,XI,YI) uses bilinear interpolation to
%   find ZI, the values of the underlying 2-D function in Z at the points
%   in matrices XI and YI.  Matrices X and Y specify the points at which 
%   the data Z is given.  X and Y can also be vectors specifying the 
%   abscissae for the matrix Z as for MESHGRID. In both cases, X
%   and Y must be equally spaced and monotonic.
%
%   Values of NaN are returned in ZI for values of XI and YI that are 
%   outside of the range of X and Y.
%
%   If XI and YI are vectors, BI_LINEAR returns vector ZI containing
%   the interpolated values at the corresponding points (XI,YI).
%
%   Clay M. Thompson 3-22-93.

	[nrows,ncols,dumb] = size(arg3);
	mx = numel(arg1);   my = numel(arg2);
	if any([mx my] ~= [ncols nrows]) && ~isequal(size(arg1),size(arg2),size(arg3))
		error('The lengths of the X and Y vectors must match Z.');
	end
	s = 1 + (arg4-arg1(1))/(arg1(mx)-arg1(1))*(ncols-1);
	t = 1 + (arg5-arg2(1))/(arg2(my)-arg2(1))*(nrows-1);
	row = round(t);
	col = round(s);

	if (~isequal(size(s),size(t)))
		error('XI and YI must be the same size.');
	end

	% Check for out of range values of s and set to 1
	sout = find((s < 1) | (s > ncols));
	if ~isempty(sout), s(sout) = ones(size(sout)); end

	% Check for out of range values of t and set to 1
	tout = find((t < 1) | (t > nrows));
	if ~isempty(tout), t(tout) = ones(size(tout)); end

	% Matrix element indexing
	ndx = floor(t)+floor(s-1)*nrows;

	% Compute intepolation parameters, check for boundary value.
	if isempty(s),	d = s;
	else			d = find(s == ncols);
	end
	s(:) = (s - floor(s));
	if ~isempty(d), s(d) = s(d)+1; ndx(d) = ndx(d)-nrows; end

	% Compute intepolation parameters, check for boundary value.
	if isempty(t),	d = t;
	else			d = find(t == nrows);
	end
	t(:) = (t - floor(t));
	if (~isempty(d))
		t(d) = t(d)+1;		ndx(d) = ndx(d)-1;
	end

	% Make sure arg3_? is of double type
	if (~isa(arg3,'double'))
		arg3_1 = double(arg3(ndx));
		arg3_2 = double(arg3(ndx+1));
		arg3_3 = double(arg3(ndx+nrows));
		arg3_4 = double(arg3(ndx+(nrows+1)));
	else
		arg3_1 = arg3(ndx);
		arg3_2 = arg3(ndx+1);
		arg3_3 = arg3(ndx+nrows);
		arg3_4 = arg3(ndx+(nrows+1));
	end

	% Now interpolate
	F = (arg3_1.*(1-t) + arg3_2.*t).*(1-s) + (arg3_3.*(1-t) + arg3_4.*t).*s;

	% Now set out of range values to NaN.
	if ~isempty(sout), F(sout) = NaN; end
	if ~isempty(tout), F(tout) = NaN; end

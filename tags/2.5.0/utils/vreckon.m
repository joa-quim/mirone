function [lat2,lon2,a21] = vreckon(lat1, lon1, rng, azim, npts, ellipsoid)
% VRECKON -  Computes points at a specified azimuth and range in an ellipsoidal earth
%
%		- Using the WGS-84 Earth ellipsoid, travel a given distance along
%		  a given azimuth starting at a given initial point, and return the
%		  endpoint within a few millimeters of accuracy, using Vincenty's algorithm.
%
% USAGE:
% [lat2,lon2] = vreckon(lat1, lon1, rng, azim)
% ... = vreckon(...,ellipsoid)
%  Transmits ellipsoid definition (either as [a,b] or [a,f]) as fifth argument ELLIPSOIDE
%
% VARIABLES:
% lat1 = inital latitude (degrees)
% lon1 = initial longitude (degrees)
% rng  = distance (meters)
%		It can be a scalar or a vector. Latter case computes a series of
%		circles (or arc circles, see azim) centered on X,Y (which are scalars)
% azim = intial azimuth (degrees)
% 		"azim" is a one or two-column vector. For single column, returns
%		the arc between 0 and azim. For two columns, returns the arc
%		between azim(1) and azim(2).
% ellipsoid = two-element ellipsoid vector. Either [a b] or [a f]
%		If omitted, defaults to WGS-84
% lat2, lon2 = second point (degrees)
% a21  = reverse azimuth (degrees), at final point facing back toward the
%        intial point
%
% Original algorithm source:
% T. Vincenty, "Direct and Inverse Solutions of Geodesics on the Ellipsoid
% with Application of Nested Equations", Survey Review, vol. 23, no. 176,
% April 1975, pp 88-93.
% Available at: http://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf
%
% Notes: 
% (1) The Vincenty reckoning algorithm was transcribed verbatim into
%     JavaScript by Chris Veness. It was modified and translated to Matlab
%     by Michael Kleder. Mr. Veness's website is:
%     http://www.movable-type.co.uk/scripts/latlong-vincenty-direct.html
% (2) Error correcting code, polar error corrections, WGS84 ellipsoid
%     parameters, testing, and comments by Michael Kleder.
% (3) By convention, when starting at a pole, the longitude of the initial
%     point (otherwise meaningless) determines the longitude line along
%     which to traverse, and hence the longitude of the final point.
% (4) The convention noted in (3) above creates a discrepancy with VDIST
%     when the the intial or final point is at a pole. In the VDIST
%     function, when traversing from a pole, the azimuth is  0 when
%     heading away from the south pole and 180 when heading away from the
%     north pole. In contrast, this VRECKON function uses the azimuth as
%     noted in (3) above when traversing away form a pole.
% (5) In testing, where the traversal subtends no more than 178 degrees,
%     this function correctly inverts the VDIST function to within 0.2
%     millimeters of distance, 5e-10 degrees of forward azimuth,
%     and 5e-10 degrees of reverse azimuth. Precision reduces as test
%     points approach antipodal because the precision of VDIST is reduced
%     for nearly antipodal points. (A warning is given by VDIST.)
% (6) Tested but no warranty. Use at your own risk.
% (7) Ver 1.0, Michael Kleder, November 2007
% (8) Ver 2.0, Joaquim Luis, September 2008
%
% Added 5th argin (ellipsoid) and vectorized whenever possible.
% It can now also compute circles or circle arcs. 
% Also, lon2 is always converted to the [-180 180] interval
% Joaquim Luis

% $Id: $

	% Input check:
	if (abs(lat1) > 90)
		error('VRECKON: Input latitude must be between -90 and 90 degrees, inclusive.')
	end
	if (numel(lat1) == 1 && numel(rng) > 1)
		error('VRECKON: Variable ranges are only allowed for a single point.')
	end

	if (nargin == 4),		npts = 180;		end
	if (nargin == 6)			% An ellipsoid vector (with a & b OR a & f) was supplyied
		a = ellipsoid(1);		% b = ellipsoid(2);
		if (ellipsoid(2) < 1)	% Second ellipsoid argument contains flattening instead of minor axis
			f = ellipsoid(2);	b = a * (1 - f);
		else					% Second ellipsoid argument contains minor axis
			f = (a - ellipsoid(2)) / a;
		end
	else
		% Supply WGS84 earth ellipsoid axis lengths in meters:
		a = 6378137; 			% semimajor axis
		b = 6356752.31424518; 	% computed from WGS84 earth flattening coefficient definition
		f = (a-b)/a;
	end

	lat1 = lat1(:) * .1745329251994329577e-1;		% intial latitude in radians
	lon1 = lon1(:) * .1745329251994329577e-1;		% intial longitude in radians

	% correct for errors at exact poles by adjusting 0.6 millimeters:
	kidx = abs(pi/2-abs(lat1)) < 1e-10;
	if any(kidx)
		lat1(kidx) = sign(lat1(kidx))*(pi/2-(1e-10));
	end
        
	%  Allow for multiple circles starting from the same point
	if (numel(lat1) == 1 && numel(lon1) == 1 && numel(rng) > 1)
		lat1 = lat1(ones(size(rng)));		lon1 = lon1(ones(size(rng)));
	end
	
	if isempty(azim),
		azim = [0 360];			azim = azim(ones([size(lat1,1) 1]), :);
	elseif (isequal(size(azim), [1 2]))
		azim = azim(ones([size(lat1,1) 1]), :);
	end
	if (numel(rng) == 1),		rng = ones(size(azim,1), 1) * rng;
	else						rng = rng(:);
	end
	if (numel(rng) ~= size(azim ,1))
		error('VRECKON: Range must be a scalar or a vector with the same size as azim(:,1).')
	end

	%  Expand the azimuth inputs
	if size(azim,2) == 1				% Single column azimuth inputs
		az_neg = zeros(size(azim));		az_pos = azim;
	else								% Two column azimuth inputs
		az_neg = azim(:,1);				az_pos = azim(:,2);
	end
	
	az = zeros([size(az_neg,1) npts]);
	for i = 1:size(az,1)
		% Handle case where limits give more than half of the circle. Go clockwise from start to end.
		if (az_neg(i) > az_pos(i)),		az_pos(i) = az_pos(i)+2*pi;		end
		az(i,:) = linspace(az_neg(i),az_pos(i),real(npts));	
	end
	
	% Each circle occupies a row of the output matrices.
	lat1 = lat1(:,ones([1,size(az,2)]) );
	lon1 = lon1(:,ones([1,size(az,2)]) );
	rng  = rng(:,ones([1,size(az,2)]) );
	
	alpha1 = az * 0.1745329251994329577e-1;	% inital azimuth in radians
	sinAlpha1 = sin(alpha1);	cosAlpha1 = cos(alpha1);
	tanU1 = (1-f) * tan(lat1);
	cosU1 = 1 ./ sqrt(1 + tanU1 .* tanU1);
	sinU1 = tanU1 .* cosU1;
	sigma1 = atan2(tanU1, cosAlpha1);
	sinAlpha = cosU1 .* sinAlpha1;
	cosSqAlpha = 1 - sinAlpha .* sinAlpha;
	uSq = cosSqAlpha * (a*a - b*b) / (b*b);
	A = 1 + uSq/16384 .* (4096+uSq .* (-768+uSq .* (320-175 * uSq)));
	B = uSq/1024 .* (256+uSq .* (-128+uSq .* (74-47 * uSq)));
	sigma = rng ./ (b*A);
	sigmaP = 2*pi;

	if (numel(sigma) == 1)
		while (abs(sigma-sigmaP) > 1e-12)
			cos2SigmaM = cos(2*sigma1 + sigma);
			sinSigma = sin(sigma);
			cosSigma = cos(sigma);
			deltaSigma = B .* sinSigma*(cos2SigmaM + B/4 .* (cosSigma .* (-1 + 2*cos2SigmaM .* cos2SigmaM)-...
				B/6 .* cos2SigmaM .* (-3+4*sinSigma .* sinSigma) .* (-3 + 4*cos2SigmaM .* cos2SigmaM)));
			sigmaP = sigma;
			sigma = rng ./ (b*A) + deltaSigma;
		end
	else
		% This part is not vectorized
		cos2SigmaM = zeros(size(sigma));		sinSigma = zeros(size(sigma));		cosSigma = zeros(size(sigma));
		for (k = 1:numel(sigma))
			while (abs(sigma(k)-sigmaP) > 1e-12)
				cos2SigmaM(k) = cos(2*sigma1(k) + sigma(k));
				sinSigma(k)   = sin(sigma(k));
				cosSigma(k)   = cos(sigma(k));
				tmp = 2*cos2SigmaM(k) * cos2SigmaM(k);
				deltaSigma = B(k) * sinSigma(k) * (cos2SigmaM(k) + B(k)/4 * (cosSigma(k) * (-1 + tmp)-...
					B(k)/6 * cos2SigmaM(k) * (-3+4*sinSigma(k) * sinSigma(k)) * (-3 + 2*tmp)));
				sigmaP = sigma(k);
				sigma(k) = rng(k) ./ (b*A(k)) + deltaSigma;
			end
		end
	end

	tmp = sinU1  .* sinSigma - cosU1 .* cosSigma .* cosAlpha1;
	lat2 = atan2(sinU1 .* cosSigma + cosU1 .* sinSigma .* cosAlpha1, (1-f)*sqrt(sinAlpha .* sinAlpha + tmp .* tmp));
	lambda = atan2(sinSigma .* sinAlpha1, cosU1 .* cosSigma - sinU1 .* sinSigma .* cosAlpha1);
	C = f/16*cosSqAlpha .* (4+f*(4-3*cosSqAlpha));
	L = lambda - f * (1-C) .* sinAlpha .* (sigma + C .* sinSigma .* (cos2SigmaM + C .* cosSigma .* (-1+2*cos2SigmaM .* cos2SigmaM)));
	lon2 = lon1 + L;
	
	% Truncates angles into the [-pi pi] range
	if (any(lon2 > pi))
		lon2 = pi*((abs(lon2)/pi) - 2*ceil(((abs(lon2)/pi)-1)/2)) .* sign(lon2);
	end

	% output degrees
	lat2 = lat2 * 57.295779513082322865;
	lon2 = lon2 * 57.295779513082322865;
	% lon2 = mod(lon2,360); % follow [0,360] convention -- SORRY, but don't like that J.L.
	if nargout > 2
		a21 = atan2(sinAlpha, -tmp); 
		a21  = 180 + a21  * 57.295779513082322865; % note direction reversal
		a21=mod(a21,360);
	end

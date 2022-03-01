function [fout, k] = rtp3d(f3d, incl_fld, decl_fld, incl_mag, decl_mag, component) 
% RTP3D  Reduce a magnetic field anomaly map to the pole using Fourier transform
%        and specifying inclination and declination of the field and magnetization
%	f3d		input array
%	fout	output array
%	k		Wavenumber array 
%
% usage: fout=rtp3d(f3d,incl_fld,decl_fld,incl_mag,decl_mag)
%
%  Maurice A. Tivey March 16 1992/ Oct 1994
%   Retouched by Joaquim Luis

	if (nargin < 6),	component = 0;	end			% I.e. default to RTP

	i = sqrt(-1);  % complex i
	D2R = pi/180;  % conversion radians to degrees

	incl_fld = incl_fld * D2R;    decl_fld = decl_fld * D2R;
	incl_mag = incl_mag * D2R;    decl_mag = decl_mag * D2R;

	[ny,nx] = size(f3d);
	ni = 1/nx;          nx2 = nx/2;
	x = -.5:ni:.5-ni;
	ni = 1/ny;          ny2 = ny/2;
	y = -.5:ni:.5-ni;

	X = ones(size(y))'*x;
	Y = y'*ones(size(x));
	if (nargout == 2)
		k = 2*pi*sqrt(X.^2+Y.^2);  % wavenumber array
	end

	aux = atan2(Y,X);		clear X	Y	% Auxiliary variable to not compute this more than once

	%------ calculate geometric and amplitude factors
	Ob = sin(incl_fld) + i*cos(incl_fld) * sin(aux + decl_fld);
	if ((abs(incl_fld - incl_mag) > 0.03) || (abs(decl_fld - decl_mag) > 0.03) || component) % 0.03 is < 2 deg
		Om = sin(incl_mag) + i*cos(incl_mag) * sin(aux + decl_mag);
		O  = Ob .* Om;
	else
		O  = Ob .* Ob;
	end
	clear Ob
	%O  = (Ob.*Om);   O = fftshift(O);    %clear Ob Om

	if (component)		% Compute one of the X, Y or Z components intead of a RTP
		if     (component == 1),	new_inc = 0;	new_dec = 0;		% X or North component
		elseif (component == 2),	new_inc = 0;	new_dec = 90;		% Y or East component
		elseif (component == 3),	new_inc = 90;	new_dec = 0;		% Z or Up component
		end
		new_inc = new_inc * D2R;	new_dec = new_dec * D2R;
		O = O ./ ((sin(new_inc) +i*(cos(new_inc) * sin(aux + new_dec)) + eps) .* Om);	% So that at the end <=> F .* ON ./ O
	end
	O = fftshift(O);

% % calculate northpole phase
%  ra1 = 90 * D2R;
%  rb1 = 0 * D2R;
%  ONP = (sin(ra1)+i*cos(ra1)*sin(aux + rb1));
%  ONP = fftshift(ONP.*Om);
	
	% save mean of input magnetic field
	mfin = mean(mean(f3d));
 	F = fft2(f3d - mfin);
	%------------------ INVERSE fft ----------------
	fout = ifft2(F ./ O);
%	fout = ifft2((F .* ONP ./ O));
	% add back in mean of input field?
	fout = real(fout); 

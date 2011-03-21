function [fout, k] = rtp3d(f3d,incl_fld,decl_fld,incl_mag,decl_mag) 
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

	i = sqrt(-1);  % complex i
	D2R = pi/180;  % conversion radians to degrees

	incl_fld = (incl_fld)*D2R;    decl_fld = (decl_fld)*D2R;
	incl_mag = (incl_mag)*D2R;    decl_mag = (decl_mag)*D2R;

	[ny,nx] = size(f3d);
	ni = 1/nx;          nx2 = nx/2;
	nx2plus = nx2+1;    x = -.5:ni:.5-ni;
	ni = 1/ny;          ny2 = ny/2;
	ny2plus = ny2+1;    y = -.5:ni:.5-ni;

	X = ones(size(y))'*x;
	Y = y'*ones(size(x));
	if (nargout == 2)
		k = 2*pi*sqrt(X.^2+Y.^2);  % wavenumber array
	end

	%------ calculate geometric and amplitude factors
	Ob = (sin(incl_fld)+i*cos(incl_fld)*sin(atan2(Y,X)+decl_fld));
	Om = (sin(incl_mag)+i*cos(incl_mag)*sin(atan2(Y,X)+decl_mag));
	O  = (Ob.*Om);   O = fftshift(O);    clear Ob Om X Y;

	% save mean of input magnetic field
	mfin = mean(mean(f3d));
	F = fft2(f3d);
	%------------------ INVERSE fft ----------------
	fout = ifft2((F./O));
	% add back in mean of input field?
	fout = real(fout)+mfin; 

function f3d = syn3d(m3d,h,rlat,rlon,yr,zobs,thick,slin,dx,dy,sdip,sdec)
% SYN3D Calculate magnetic field given a magnetization and bathymetry
% map using Parker's [1973] Fourier series summation approach
% Input arrays:
%    m3d 	magnetization (A/m)
%    h 		bathymetry (km +ve up)
%    rlat 	latitude of survey area (dec. deg.)
%    rlon 	longitude of survey area (dec. deg.)
%    yr 	year of survey (dec. year)
%    slin 	azimuth of lineations (deg) set to 0
%    zobs 	observation level (+km up)
%    thick 	thickness of source layer (km)
%    dx 	x grid spacing  (km) 
%    dy		y grid spacing  (km)
%    sdec	declination of magnetization (optional)
%    sdip	inclination of magnetization (optional)
% Output array:
%    f3d	magnetic field (nT)
%
%
% Usage: f3d=syn3d(m3d,h,rlat,rlon,yr,zobs,thick,slin,dx,dy,sdip,sdec)
%   or geocentric dipole: 
%        f3d=syn3d(m3d,h,rlat,rlon,yr,zobs,thick,slin,dx,dy);
%
% Maurice A. Tivey MATLAB Version 5 August 1992
%                                   March  1996
%  Joaquim Luis         May   2004
%       Cleaned the code for saving RAM and replaced the igrf routine
%-----------------------------------------------------------

	if (nargin == 0)
% 		fprintf('\n\n             DEMO OF SYN3D mfile\n\n');
% 		fprintf('GENERATE A PRISM AT DEPTH\n');
% 		m=zeros(32,32);  m(16:20,16:20)=m(16:20,16:20)+10;
% 		fprintf('SET PARAMETERS FOR CALCULATING THE FIELD\n');
% 		fprintf('ASSUME TOPOGRAPHY IS FLAT\\nn');
% 		h=ones(32,32).*(-2);
% 		dx=1; dy=1; sdip=0; sdec=0; thick=1;
% 		rlat=26;rlon=-45;yr=1990;slin=0;zobs=0;
% 		f3d=syn3d(m,h,rlat,rlon,yr,zobs,thick,slin,dx,dy,sdip,sdec);
% 		return
	end
	if (nargin <= 10) % geocentric dipole hypothesis assumed
		sdip = 0;    sdec = 0;
	end  
	% parameters defined
	D2R = pi/180;  % conversion radians to degrees
	mu = 100;      % conversion factor to nT

	% changeable parameters 
	nterms = 8;    tol = 0.1;

	% If a window for verbose exists
	old_show = get(0,'ShowHiddenHandles');
	set(0,'ShowHiddenHandles','on')
	h = get(0,'Children'); 
	figdmsg = findobj(h,'flat','tag','Wdmsgfig');
	set(0,'ShowHiddenHandles',old_show)

	if (isempty(figdmsg))
		try		message_win('create','     3D MAGNETIC FIELD FORWARD MODEL');    pause(0.001)
		catch,  fprintf('     3D MAGNETIC FIELD FORWARD MODEL\n')
		end
	else
		try		message_win('add','     3D MAGNETIC FIELD FORWARD MODEL')
		catch,  fprintf('     3D MAGNETIC FIELD FORWARD MODEL\n')
		end
	end

	try	message_win('add','       Constant thickness layer');
	catch,  fprintf('       Constant thickness layer\n')
	end
	try		message_win('add',' M.A.Tivey      Version: May, 2004')
	catch,  fprintf(' M.A.Tivey      Version: May, 2004\n')
	end
	try		message_win('add',sprintf(' Zobs = %12.5f\n Rlat = %12.4f\t Rlon = %12.4f',zobs,rlat,rlon))
	catch,  fprintf(' Zobs= %12.5f\n Rlat = %12.4f\t Rlon = %12.4f\n',zobs,rlat,rlon)
	end
	try		message_win('add',sprintf(' Yr = %12.4f',yr))
	catch,  fprintf(' Yr = %12.4f\n',yr)
	end
	try		message_win('add',sprintf(' Thick = %12.4f',thick))
	catch,  fprintf(' Thick = %12.4f\n',thick)
	end
	try		message_win('add',sprintf(' Slin,Sdec,Sdip = %12.6f %12.6f %12.6f',slin,sdec,sdip))
	catch,  fprintf(' Slin,Sdec,Sdip = %12.6f %12.6f %12.6f\n',slin,sdec,sdip)
	end
	try		message_win('add',sprintf(' Nterms,Tol %6.0f %10.5f',nterms,tol))
	catch,  fprintf(' Nterms,Tol %6.0f %10.5f \n',nterms,tol)
	end
	[ny,nx] = size(m3d);
	try		message_win('add',sprintf(' Number of points in map are : %6.0f x%6.0f',nx,ny));
	catch,  fprintf(' Number of points in map are : %6.0f x%6.0f\n',nx,ny)
	end
	try		message_win('add',sprintf(' Spacing of points : %10.4f X %10.4f',dx,dy));
	catch,  fprintf(' Spacing of points : %10.4f X %10.4f \n',dx,dy)
	end

	%y = igrf(rlat, rlon, zobs, yr, 'igrf_1945_2005.dat');
	y = igrf_m(rlon, rlat, zobs, yr);

	% compute skewness parameter
	decl1 = y(6);   incl1 = y(7);   clear y;
	if (abs(sdec) > 0. || abs(sdip) > 0.)
		%[theta,ampfac] = nskew(yr,rlat,rlon,zobs,slin,sdec,sdip);
	else
		%[theta,ampfac] = nskew(yr,rlat,rlon,zobs,slin);
		sdip = atan2( 2.*sin(rlat*D2R),cos(rlat*D2R) )/D2R;
		sdec = 0;
	end
	%
	slin = 0;           % slin is forced to zero
	ra1 = incl1*D2R;    rb1 = (decl1-slin)*D2R;
	ra2 = sdip*D2R;     rb2 = (sdec-slin)*D2R;

	% calculate wavenumber array
	% nx2 = nx/2;                 ny2 = ny/2;       % And if nx|ny are odd?
	% dkx = pi/(nx*dx);           dky = pi/(ny*dy);
	% kx = (-nx2:nx2-1).*dkx;     ky = (-ny2:ny2-1).*dky;
	% X = ones(size(ky))'*kx;     Y = ky'*ones(size(kx));
	% k = fftshift(2*sqrt(X.^2+Y.^2));  % wavenumber array

	nx2 = fix(nx/2);		ny2 = fix(ny/2);
	if (rem(nx,2) == 0),	sft_x = 1;
	else					sft_x = 0;
	end
	if (rem(ny,2) == 0),	sft_y = 1;
	else					sft_y = 0;
	end
	dkx = 2*pi / (nx*dx);   dky = 2*pi / (ny*dy);
	kx = (-nx2:nx2-sft_x).*dkx;     ky = (-ny2:ny2-sft_y).*dky;
	X = repmat(kx,length(ky),1);    Y = repmat(ky',1,length(kx));
	k = ifftshift(sqrt(X.^2+Y.^2));      % wavenumber array

	%
	i = sqrt(-1);
	aux = atan2(Y,X);       clear X Y;
	Ob = (sin(ra1)+i*cos(ra1)*sin(aux+rb1));
	Om = (sin(ra2)+i*cos(ra2)*sin(aux+rb2));
	clear aux;
	O = fftshift(Ob.*Om);
	%amp = abs(O);   % amplitude factor
	clear Ob Om;
	const = 2*pi*mu;

	% shift zero level of bathy
	hmax = max(max(h));     hmin = min(min(h));
	try		message_win('add',sprintf(' %10.3f %10.3f = MIN, MAX OBSERVED BATHY',hmin,hmax))
	catch,  fprintf(' %10.3f %10.3f = MIN, MAX OBSERVED BATHY\n',hmin,hmax)
	end
	shift = max(max(h));
	hwiggl = abs(hmax-hmin)/2;
	zup = zobs-shift;
	try		message_win('add',sprintf(' SHIFT ZERO OF BATHY WILL BE %8.3f',shift))
	catch,  fprintf(' SHIFT ZERO OF BATHY WILL BE %8.3f\n',shift)
	end
	try		message_win('add',' THIS IS OPTIMUM FOR INVERSION')
	catch,  fprintf(' THIS IS OPTIMUM FOR INVERSION.\n')
	end
	try		message_win('add',sprintf(' NOTE OBSERVATIONS ARE %8.3f KM ABOVE BATHY',zup))
	catch,  fprintf(' NOTE OBSERVATIONS ARE %8.3f KM ABOVE BATHY\n',zup)
	end
	try		message_win('add',sprintf('ZOBS=%8.3f\t ZUP=%8.3f',zobs,zup))
	catch,  fprintf('ZOBS=%8.3f ZUP=%8.3f\n',zobs,zup)
	end
	zup = zup+hwiggl;     
	try		message_win('add',sprintf('%8.3f = HWIGGL, DISTANCE TO MID-LINE OF BATHY',hwiggl))
	catch,  fprintf('%8.3f = HWIGGL, DISTANCE TO MID-LINE OF BATHY\n',hwiggl)
	end
	try		message_win('add',' THIS IS OPTIMUM ZERO LEVEL FOR FORWARD PROBLEM')
	catch,  fprintf(' THIS IS OPTIMUM ZERO LEVEL FOR FORWARD PROBLEM\n')
	end
	h = h-shift;
	h = h+hwiggl;

	eterm = exp(-k.*zup);       % do upcon term
	% now do summing over nterms
	MH = fft2(m3d);
	msum1 = eterm.*MH;
	last = 0;
	%first = max(abs(msum1));        % NOT USED ?
	for (n = 1:nterms)
		MH = (fft2(m3d.*h.^n));
		msum = eterm.*((k.^n)./factorial(n)).*MH+msum1;
		errmax = max(max(abs(real(msum))));
		try		message_win('add',sprintf(' AT TERM  %6.0f\t\tMAXIMUM PERTURBATION TO SUM %12.6e',n,errmax-last))
		catch,  fprintf(' AT TERM  %6.0f\t\t MAXIMUM PERTURBATION TO SUM %12.6e\n',n,errmax-last)
		end
		last = errmax;
		msum1 = msum;
	end
	clear MH eterm msum1;
	alap = 1 - exp(-k.*thick);  % do thickness term
	f3d = real(ifft2((const.*msum.*alap.*O)));


% --------------------------------------------------------------------------
function [theta,ampfac] = nskew(yr,rlat,rlon,zobs,slin,sdec,sdip)

% NSKEW  Compute skewness parameter and amplitude factor
%  following Schouten (1971)
%  NEW VERSION THAT COMPUTES GEOCENTRIC DIPOLE UNLESS GIVEN
%  DEC AND DIP OF MAGNETIZATION      
% Usage: [theta,ampfac]=skew(yr,rlat,rlon,zobs,slin,sdec,sdip)
%  Input variables:
%   YR: decimal year of survey
%   RLAT: regional decimal latitude degrees
%   RLON: regional decimal longitude degrees
%   ZOBS: level of observation in km above sealevel
%   SLIN: strike of lineations normal to profile
%   SDEC: magnetization declination cw degrees from north
%   SDIP: magnetization inclination degrees
% Output variables
%   THETA: phase angle
%   AMPFAC: amplitude factor
% Calls <magfd>
% Maurice A. Tivey February 3, 1993
%    checked April 1996
%  Joaquim Luis         May   2004
%       Cleaned a bit of the code and replaced the igrf routine
%---------------------------------------------------------
	D2R=pi/180;
	% get unit vectors
	%y = igrf(rlat, rlon, zobs, yr, 'igrf_1945_2005.dat');
	y = igrf_m(rlon, rlat, zobs, yr);
	decl1 = y(6);   incl1 = y(7);   clear y;

	% compute skewness parameter
	% If a window for verbose exists
	old_show = get(0,'ShowHiddenHandles');
	set(0,'ShowHiddenHandles','on')
	h = get(0,'Children'); 
	figdmsg = findobj(h,'flat','tag','Wdmsgfig');
	set(0,'ShowHiddenHandles',old_show)

	if (isempty(figdmsg))
		message_win('create',' EARTH''S MAGNETIC FIELD DIRECTION:');    pause(0.001)
	else
		message_win('add',' EARTH''S MAGNETIC FIELD DIRECTION:');
	end
	message_win('add',sprintf(' %10.3f = MAGNETIC DECLINATION ( STRIKE, CW FROM N )',decl1));
	message_win('add',sprintf(' %10.4f = MAGNETIC INCLINATION ( DIP, POS DOWN )',incl1));
	if (nargin > 5)
		%  NOTE FOR GEOCENTRIC DIPOLE TAN(INC)=2*TAN(LAT)
		message_win('add',' NON-GEOCENTRIC MAGNETIZATION VECTOR SPECIFIED:');
		message_win('add',sprintf(' %10.4f = DESIRED MAGNETIZATION DECLINATION (+CW FROM N)',sdec));
		message_win('add',sprintf(' %10.4f = DESIRED MAGNETIZATION INCLINATION (+DN)', sdip));
	else %
		sdip = atan2( 2.*sin(rlat*D2R),cos(rlat*D2R) )/D2R;
		sdec = 0;
		message_win('add',' GEOCENTRIC MAGNETIZATION VECTOR SPECIFIED:');
		message_win('add',sprintf(' %10.4f = GEOCENTRIC DIPOLE INCLINATION',sdip));
		message_win('add',sprintf(' %10.3f = GEOCENTRIC DECLINATION ASSUMED\n',sdec));
	end

	% compute phase and amplitude factors
	ra1 = incl1*D2R;    rb1 = (decl1-slin)*D2R;
	ra2 = sdip*D2R;     rb2 = (sdec-slin)*D2R;
	% compute phase and amplitude factors
	inclm = atan2(tan(ra2),sin(rb2));
	inclf = atan2(tan(ra1),sin(rb1));
	ampfac= ((sin(ra2))*(sin(ra1)))/((sin(inclm))*(sin(inclf)));
	theta = (inclm/D2R)+(inclf/D2R)-180;
	% compute unit vectors for a check
	hatm(1) = cos(sdip*D2R)*sin((sdec-slin)*D2R);
	hatm(2) = cos(sdip*D2R)*cos((sdec-slin)*D2R);
	hatm(3) = -sin(sdip*D2R);
	hatb(1) = cos(incl1*D2R)*sin((decl1-slin)*D2R);
	hatb(2) = cos(incl1*D2R)*cos((decl1-slin)*D2R);
	hatb(3) = -sin(incl1*D2R);

	message_win('add',sprintf('  %10.6f %10.6f %10.6f = MAGNETIZATION UNIT VECTOR',hatm(1),hatm(2),hatm(3)));
	message_win('add',sprintf('  %10.6f %10.6f %10.6f = AMBIENT FIELD UNIT VECTOR',hatb(1),hatb(2),hatb(3)));
	message_win('add','  COMPONENTS ARE (X,Y,Z=ALONG, ACROSS PROFILE, AND UP');

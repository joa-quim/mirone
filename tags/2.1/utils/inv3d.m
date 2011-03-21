function m3d=inv3d(f3d,h,wl,ws,rlat,rlon,yr,zobs,thick,slin,dx,dy,sdec,sdip);
% INV3D Calculate magnetization from magnetic field and bathymetry for a map.
%       Assumes constant thickness source layer whose upper bound is bathymetry
%       Use the Parker & Huestis [1974] Fourier inversion approach.  
% Input arrays:
%    f3d 	magnetic field (nT)
%    h 		bathymetry (km +ve up)
%    wl		filter long wavelength cutoff (km)
%    ws		filter short wavelength cutoff (km)
%    rlat 	latitude of survey area dec. deg.
%    rlon 	longitude of survey area dec. deg.
%    yr 	year of survey (dec. year)
%    slin 	azimuth of lineations (deg)
%    zobs 	observation level (+km up)
%    thick 	thickness of source layer (km)
%    azim	azimuth of grid (degrees) hard wired to 0
%    dx 	x grid spacing  (km)
%    dy 	y grid spacing  (km)
%    sdec	declination of magnetization (optional)
%    sdip	inclination of magnetization (optional)
% Output array:
%    m3d	magnetization (A/m)
%
% Usage: m3d=inv3d(f3d,h,wl,ws,rlat,rlon,yr,zobs,thick,azim,dx,dy,sdec,sdip);
%   or for geocentric dipole
%     m3d=inv3d(f3d,h,wl,ws,rlat,rlon,yr,zobs,thick,azim,dx,dy);
%
% Maurice Tivey  Mar 27 1996
% Joaquim Luis   May   2004
%   Cleaned the code for saving RAM and replaced the igrf routine (the original was to damn slow)
%   May 2005 -> Added an option to deal with case of the field is already RTP
%               In order to not change the nargins, this option is activated with:
%               sdec = rlon = 0 & sdip = rlat = 90
%-----------------------------------------------------------

if (nargin == 0)
%     try,    message_win('create','DEMO OF INV3D');  pause(0.001);
%             message_win('add','FIRST GENERATE A PRISM AT DEPTH');
%             message_win('add','SET PARAMETERS FOR CALCULATING THE FIELD');
%             message_win('add','ASSUME TOPOGRAPHY IS FLAT');
%     catch,  fprintf('\n            DEMO OF INV3D mfile\n');
%             fprintf('FIRST GENERATE A PRISM AT DEPTH\n');
%             fprintf('SET PARAMETERS FOR CALCULATING THE FIELD\n');
%             fprintf('ASSUME TOPOGRAPHY IS FLAT\n');
%     end
%     m=zeros(64,64);  m(32:40,32:40)=m(32:40,32:40)+10;
%     h=ones(64,64).*(-2);
%     dx=.5; dy=0.5; sdip=0;; sdec=0; thick=1;
%     rlat=26;rlon=-45;yr=1990;slin=0;zobs=0;
%     f3d=syn3d(m,h,rlat,rlon,yr,zobs,thick,slin,dx,dy,sdip,sdec);
%     try,    message_win('add','NOW DO INVERSION.....');
%     catch,  fprintf('NOW DO INVERSION.....\n');  end
%     wl = 0; ws = 0;
%     m3d=inv3d(f3d,h,wl,ws,rlat,rlon,yr,zobs,thick,slin,dx,dy);
%     surf(m3d);view(-30,65);
    return
end

if (nargin > 12), % user defined sdip sdec
    pflag = 1;
else            % geocentric dipole hypothesis assumed
    sdip = 0;    sdec = 0;    pflag = 0;
end  
% parameters defined
i = sqrt(-1);
D2R = pi/180;  % conversion D2Rians to degrees
mu = 100;      % conversion factor to nT
% changeable parameters
nterms = 8;     nitrs = 10;
tol = 0.0001;   tolmag = 0.0001;
flag = 0;       xmin = 0;

[ny,nx] = size(f3d);
if size(h) ~= size(f3d), 
    fprintf(' bathy and field arrays must e of the same length\n');
    return; 
end
% remove mean from input field
mnf3d = mean(mean(f3d));
f3d = f3d - mnf3d;

%y = igrf(rlat, rlon, zobs, yr, 'igrf_1945_2005.dat');
y = igrf_m(rlon, rlat, zobs, yr);

% compute phase and amplitude factors from 2D method
%bx = y(3);      by = y(4);      bz = y(5);      bh = y(2);
decl = y(6);    incl = y(7);    clear y;
%if (abs(sdec) > 0 | abs(sdip) > 0)
if (abs(sdec) == 0 & abs(sdip) == 90)
    %[theta,ampfac] = nskew(yr,rlat,rlon,zobs,slin,sdec,sdip);
    if (abs(rlat) == 90 & abs(rlon) == 0)       % Trick used to inform that the Field is RTP
        incl = 90;      decl = 0;
    end
else
    %[theta,ampfac] = nskew(yr,rlat,rlon,zobs,slin);        % Not used
    sdip = atan2( 2.*sin(rlat*D2R),cos(rlat*D2R) )/D2R;
    sdec = 0;
end
%
slin= 0;            % slin is forced to zero
ra1 = incl*D2R;     rb1 = (decl-slin)*D2R;
ra2 = sdip*D2R;     rb2 = (sdec-slin)*D2R;

% make wave number array
% nx2 = nx/2;                 ny2 = ny/2;       % And if nx|ny are odd?
% dkx = pi/(nx*dx);           dky = pi/(ny*dy);
% kx = (-nx2:nx2-1).*dkx;     ky = (-ny2:ny2-1).*dky;
% X = ones(size(ky))'*kx;     Y = ky'*ones(size(kx));
% k = fftshift(2*sqrt(X.^2+Y.^2));  % wavenumber array

nx2 = fix(nx/2);        ny2 = fix(ny/2);
if (rem(nx,2) == 0)     sft_x = 1;
else                    sft_x = 0;     end
if (rem(ny,2) == 0)     sft_y = 1;
else                    sft_y = 0;     end
dkx = 2*pi / (nx*dx);   dky = 2*pi / (ny*dy);
kx = (-nx2:nx2-sft_x).*dkx;     ky = (-ny2:ny2-sft_y).*dky;
X = repmat(kx,length(ky),1);    Y = repmat(ky',1,length(kx));
k = ifftshift(sqrt(X.^2+Y.^2));      % wavenumber array
%
%------ calculate geometric and amplitude factors
aux = atan2(Y,X);       clear X Y;
Ob = (sin(ra1)+i*cos(ra1)*sin(aux+rb1));
Om = (sin(ra2)+i*cos(ra2)*sin(aux+rb2));
clear aux;
amp = abs(fftshift(Ob.*Om));
phase = fftshift( exp(i*(angle(Ob)+angle(Om))));    % phase angle
clear Ob Om;
const = 2*pi*mu;

% shift zero level of bathy
hmax = max(max(h)); hmin = min(min(h));
conv = 1;
shift = hmax;
hwiggl = abs(hmax-hmin)/2;
zup = zobs - shift;

% bathy zero placed halfway between extremes this is optimum for summation
% but not for iteration which needs zero at highest point of bathy
h = h - shift + hwiggl;

% set up bandpass filter
wts = bpass3d(nx,ny,dx,dy,wl,ws);
% do eterm
dexpz = exp(k.*zup);    dexpw = exp(-k.*hwiggl);
alap = (1-exp(-k.*thick));      alap(1,1) = eps;      % do thickness term
% take fft of observed magnetic field and initial m3d
m3d = zeros(ny,nx); % make an initial guess of 0 for m3d
%sum1 = fft2(m3d);  % Not used
F = (fft2(f3d));
% now do summing over nterms
intsum = 0;
mlast = zeros(ny,nx);
lastm3d = zeros(ny,nx);
B = (F.*dexpz) ./ (const.*alap.*amp.*phase);
B(1,1) = 0;
clear F dexpz alap amp phase;

for iter=1:nitrs,       % summation loop
    sum = zeros(ny,nx);
    for nkount=1:nterms,
        n = nkount-1;
        MH = (fft2(m3d.*(h).^n));
        dsum = dexpw.*((k.^n)./factorial(n)).*MH;
        sum = dsum+sum;
        errmax = max(max( abs(real(sum)+imag(sum)) ));
    end
    % transform to get new solution
    M = (B-(sum))+mlast;
    % filter before transforming to ensure no blow ups
    M(1,1) = 0;
    mlast = M.*wts;
    m3d = ifft2(mlast);
    % do convergence test
    errmax = 0;
    dif_max = max(max(abs(lastm3d-m3d)));
    dif_avg = mean(mean(abs(lastm3d-m3d)));
    if (errmax-dif_max < 0)
        errmax = dif_max;
    end
    lastm3d = m3d;
    if (iter == 1) 
        first1 = errmax+1e-10; 
        erpast = errmax;
    end
    if (errmax > erpast)
        flag = 1;  % set the flag to show diverging solution
        break
    end
    erpast = errmax;
    if (errmax < tolmag)    % test for errmax less than tolerance
        flag = 0;        break
    end
    try,    message_win('add',sprintf('%3.0f, %10.4e, %6.0f %10.4e\n',iter,errmax,nkount,dif_avg));
    catch,  fprintf('%3.0f, %10.4e, %6.0f %10.4e\n',iter,errmax,nkount,dif_avg);    end
end  % end of iteration loop

if (flag == 1)
    warndlg('I would be quitting now error < tolerance ','Warning');
end

m3d = real(m3d);

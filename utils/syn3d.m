function f3d=syn3d(m3d,h,rlat,rlon,yr,zobs,thick,slin,dx,dy,sdip,sdec)
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
%     fprintf('\n\n             DEMO OF SYN3D mfile\n\n');
%     fprintf('GENERATE A PRISM AT DEPTH\n');
%     m=zeros(32,32);  m(16:20,16:20)=m(16:20,16:20)+10;
%     fprintf('SET PARAMETERS FOR CALCULATING THE FIELD\n');
%     fprintf('ASSUME TOPOGRAPHY IS FLAT\\nn');
%     h=ones(32,32).*(-2);
%     dx=1; dy=1; sdip=0; sdec=0; thick=1;
%     rlat=26;rlon=-45;yr=1990;slin=0;zobs=0;
%     f3d=syn3d(m,h,rlat,rlon,yr,zobs,thick,slin,dx,dy,sdip,sdec);
    return
end
if (nargin > 10) % user defined sdip sdec
    pflag=1;
else            % geocentric dipole hypothesis assumed
    sdip = 0;    sdec = 0;    pflag = 0;
end  
% parameters defined
i = sqrt(-1);
D2R = pi/180;  % conversion radians to degrees
mu = 100;      % conversion factor to nT

% changeable parameters 
nterms = 8;    tol = 0.1;

% If a window for verbose exists
figdmsg = wfindobj('figure','tag','Wdmsgfig');
%figdmsg = findobj(get(0,'Children'),'flat','tag','Wdmsgfig');
if (isempty(figdmsg))
    try,    message_win('create','     3D MAGNETIC FIELD FORWARD MODEL');    pause(0.001)
    catch,  fprintf('     3D MAGNETIC FIELD FORWARD MODEL\n');  end
else
    try,    message_win('add','     3D MAGNETIC FIELD FORWARD MODEL');
    catch,  fprintf('     3D MAGNETIC FIELD FORWARD MODEL\n');  end
end

try,    message_win('add','       Constant thickness layer');
catch,  fprintf('       Constant thickness layer\n');   end
try,    message_win('add',' M.A.Tivey      Version: May, 2004');
catch,  fprintf(' M.A.Tivey      Version: May, 2004\n');    end
try,    message_win('add',sprintf(' Zobs = %12.5f\n Rlat = %12.4f\t Rlon = %12.4f',zobs,rlat,rlon));
catch,  fprintf(' Zobs= %12.5f\n Rlat = %12.4f\t Rlon = %12.4f\n',zobs,rlat,rlon);    end
try,    message_win('add',sprintf(' Yr = %12.4f',yr));
catch,  fprintf(' Yr = %12.4f\n',yr);   end
try,    message_win('add',sprintf(' Thick = %12.4f',thick));
catch,  fprintf(' Thick = %12.4f\n',thick); end
try,    message_win('add',sprintf(' Slin,Sdec,Sdip = %12.6f %12.6f %12.6f',slin,sdec,sdip));
catch,  fprintf(' Slin,Sdec,Sdip = %12.6f %12.6f %12.6f\n',slin,sdec,sdip); end
try,    message_win('add',sprintf(' Nterms,Tol %6.0f %10.5f',nterms,tol));
catch,  fprintf(' Nterms,Tol %6.0f %10.5f \n',nterms,tol);  end
[ny,nx] = size(m3d);
try,    message_win('add',sprintf(' Number of points in map are : %6.0f x%6.0f',nx,ny));
catch,  fprintf(' Number of points in map are : %6.0f x%6.0f\n',nx,ny); end
try,    message_win('add',sprintf(' Spacing of points : %10.4f X %10.4f',dx,dy));
catch,  fprintf(' Spacing of points : %10.4f X %10.4f \n',dx,dy);   end

%y = igrf(rlat, rlon, zobs, yr, 'igrf_1945_2005.dat');
y = igrf_m(rlon, rlat, zobs, yr);

% compute skewness parameter
bx = y(3);      by = y(4);      bz = y(5);      bh = y(2);
decl1 = y(6);   incl1 = y(7);   clear y;
if abs(sdec) > 0. | abs(sdip) > 0.
    [theta,ampfac] = nskew(yr,rlat,rlon,zobs,slin,sdec,sdip);
else
    [theta,ampfac] = nskew(yr,rlat,rlon,zobs,slin);
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
try,    message_win('add',sprintf(' %10.3f %10.3f = MIN, MAX OBSERVED BATHY',hmin,hmax));
catch,  fprintf(' %10.3f %10.3f = MIN, MAX OBSERVED BATHY\n',hmin,hmax);    end
shift = max(max(h));
hwiggl = abs(hmax-hmin)/2;
zup = zobs-shift;
try,    message_win('add',sprintf(' SHIFT ZERO OF BATHY WILL BE %8.3f',shift));
catch,  fprintf(' SHIFT ZERO OF BATHY WILL BE %8.3f\n',shift);  end
try,    message_win('add',' THIS IS OPTIMUM FOR INVERSION');
catch,  fprintf(' THIS IS OPTIMUM FOR INVERSION.\n');   end
try,    message_win('add',sprintf(' NOTE OBSERVATIONS ARE %8.3f KM ABOVE BATHY',zup));
catch,  fprintf(' NOTE OBSERVATIONS ARE %8.3f KM ABOVE BATHY\n',zup);   end
try,    message_win('add',sprintf('ZOBS=%8.3f\t ZUP=%8.3f',zobs,zup));
catch,  fprintf('ZOBS=%8.3f ZUP=%8.3f\n',zobs,zup); end
zup = zup+hwiggl;     
try,    message_win('add',sprintf('%8.3f = HWIGGL, DISTANCE TO MID-LINE OF BATHY',hwiggl));
catch,  fprintf('%8.3f = HWIGGL, DISTANCE TO MID-LINE OF BATHY\n',hwiggl);  end
try,    message_win('add',' THIS IS OPTIMUM ZERO LEVEL FOR FORWARD PROBLEM');
catch,  fprintf(' THIS IS OPTIMUM ZERO LEVEL FOR FORWARD PROBLEM\n');   end
h = h-shift;
h = h+hwiggl;

eterm = exp(-k.*zup);       % do upcon term
% now do summing over nterms
MH = fft2(m3d);
msum1 = eterm.*MH;
last = 0;
%first = max(abs(msum1));        % NOT USED ?
for n=1:nterms,
    MH = (fft2(m3d.*h.^n));
    msum = eterm.*((k.^n)./factorial(n)).*MH+msum1;
    errmax = max(max(abs(real(msum))));
    try,    message_win('add',sprintf(' AT TERM  %6.0f\t\tMAXIMUM PERTURBATION TO SUM %12.6e',n,errmax-last));
    catch,  fprintf(' AT TERM  %6.0f\t\t MAXIMUM PERTURBATION TO SUM %12.6e\n',n,errmax-last);  end
    last = errmax;
    msum1 = msum;
end
clear MH eterm msum1;
alap = 1 - exp(-k.*thick);  % do thickness term
f3d = real(ifft2((const.*msum.*alap.*O)));

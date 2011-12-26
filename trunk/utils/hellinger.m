function [along, alat, rho] = hellinger(along, alat, rho, isow, isoe, DP_tol)
% c simplified September 9, 2000, July 2001
% c
% c program implementing section 1 of 'on reconstructing tectonic plate
% c   motion from ship track crossings'
% c
% c input data file: file code 10
% c   first line is number of sections
% c   data lines follow: each line to represent one crossing and contain
% c      iside--side number (either 1--moving or 2--fixed)
% c      isect--section number
% c      crossing latitude
% c      crossing longitude
% c      estimated error in crossing location (expressed in kilometers)
% c
% c this program estimates the rotation to rotate side 1 into side 2.  
% c
% c starting with an initial guess of the optimal rotation, the program
% c does a grid search of all rotations within .2 radians of the initial
% c guess to refine the guess.
% c it then uses subroutine amoeba, a downhill simplex method, for
% c minimization.  it is recommended that amoeba be called at least twice
% c and at least enough times to produce a minimum which is different from
% c the initial guess.
% c
% c***********************************************************************

	global eta etai

% calculating matrices sigma and sigma-tilde

	if (nargin == 0)
		alat = 64;		along = 135;	rho = 0.66;
	end

% 	isow = load('c:\j\chang\a2w.dat');
% 	isoe = load('c:\j\chang\a2e.dat');
	if (nargin < 6),	DP_tol = 0.05;		end
	%[lat, lon] = reducem_m(isow(:,2), isow(:,1), DP_tol);
	DP_tol = DP_tol * 111;			% Crude conversion of TOL from degrees to km
	B = cvlib_mex('dp',isow,DP_tol,'GEOG');
	lat = B(:,2);	lon = B(:,1);
	
	isow_dp = [lat lon];
	[c,ind] = intersect(isow(:,2),lat);
	ind = ind(end:-1:1);
	flags_w = zeros(size(isow,1),1);
	for (k = 1:numel(ind)-1)
		flags_w(ind(k):ind(k+1)-1) = k;
	end
	flags_w(end) = numel(ind) - 1;

	[r_lon,r_lat] = rot_euler(isow_dp(:,2),isow_dp(:,1), along, alat, rho, -1);		% Rotate DP moving isoc
	ndata = size(isoe,1);
	n_pts_dp = size(isow_dp,1);
	flags_e = zeros(ndata,1);
	for (k = 1:ndata)						% Loop over points of the static isoc
		P  = isoe(k,1:2);
		Dsts = sqrt((P(1)-r_lon).^2 + (P(2)-r_lat).^2);
		[D,ind] = min(Dsts);
		if ~(ind == 1 || ind == n_pts_dp)
			Q1 = [r_lon(ind-1) r_lat(ind-1)];		Q2 = [r_lon(ind) r_lat(ind)];
			D1 = abs(det([Q2-Q1; P-Q1])) / norm(Q2-Q1);
			Q1 = [r_lon(ind) r_lat(ind)];			Q2 = [r_lon(ind+1) r_lat(ind+1)];
			D2 = abs(det([Q2-Q1; P-Q1])) / norm(Q2-Q1);
		else
			D1 = 1;		D2 = 2;		% dumb values for D1 < D2
		end
		if (D1 < D2)
			flags_e(k) = ind;
		else
			flags_e(k) = ind + 1;
		end
	end
	data = [ones(size(isow,1),1) flags_w isow(:,2:-1:1) ones(size(isow,1),1); ones(size(isoe,1),1)*2 flags_e isoe(:,2:-1:1) ones(ndata,1)];

% 	fname = 'c:\j\chang\a2hell_lin.dat';
% 	data = load(fname);
	
	eps_ = 1e-4;
	nsig = 4;
	maxfn = 1000;
	ndata = size(data,1);
	nsect = max(data(:,2));
	sigma = zeros(2,nsect,3,3);
	msig  = 2*nsect*nsect+7*nsect+6;
	msig2 = 2*nsect+3;
	plev = .95;				% Confidence level
	eta = zeros(nsect,3,3);
	etai= zeros(nsect,3,2);
	  
	[x,y,z] = trans1(data(:,3), data(:,4));
	axis = [x y z];
	sd = (data(:,5) / 6366.197724) .^2;

	for (k = 1:ndata)
		for (j1 = 1:3)
			for (j2 = 1:3)
				sigma(data(k,1), data(k,2), j1, j2) = ...
					sigma(data(k,1), data(k,2), j1, j2) + axis(k,j1) * axis(k,j2) / sd(k);
			end
		end
	end	

% minimization section--grid search
	qhati = trans3(alat,along,rho);
	h = zeros(3,1);
	%[qhati, eps_, rmin] = grds(eps_, sigma, qhati, nsect);

% 250   write(6,*) 'Initial guess: alat, along, rho? '
%       read(5,*) alat,along,rho
%       write(6,*) 'Search radius (radians)? '
%       read(5,*) eps
%       call trans3(alat,along,rho,qhati)
%       write(6,*) 'qhat: ',qhati
%       do 207 i=1,3
% 207   h(i)=0.
%       rmin=r1(h)
% 251   write(6,*) 'do you want a grid search, eps = ',eps
%       read(5,'(a)') ans
%       if ((ans.ne.'Y').and.(ans.ne.'y')) go to 260
%       write(6,*) 'calling grid search--r = ',rmin
%       call grds(qhati,eps,rmin)
%       if (rmin.gt.1.e30) then
%          write(6,*) 'error return from grid--try another initial guess'
%          go to 250
%       endif
%       write(6,*) 'return from grid search--r = ',rmin
%       write(6,*) 'qhat: ',qhati
%       go to 251

% minimization section--calling amoeba

	count = 1;	rprev = 1e30;		rmin = 0;	r_fac = rmin / rprev;
	while (r_fac < 0.99 && count <= 5)
		hp = [zeros(3,1) eye(3,3) * eps_];
		for (k = 1:4)
			yp(k) = r1(hp(:,k), sigma, qhati, nsect);
		end
	
		rmin = yp(1);
		ftol = (.1)^(2*nsig+1);
		iter = maxfn;
	
		rprev = rmin;
		[hp, yp] = amoeba(hp,yp,3,3,ftol,iter, 'r1', sigma, qhati, nsect);
		rmin = r1(hp(:,1), sigma, qhati, nsect);
		qhat = trans2(hp(:,1),qhati);
		[alat, along, rho] = trans5(qhat);
		eps_ = eps_ / 20;
		qhati = qhat;
		r_fac = rmin / rprev;
		count = count + 1;
	end
	ahat = trans4(qhat);
	if (nargout == 0)
		disp(['Fitted rotation: ' num2str(alat) '  ' num2str(along) '  ' num2str(rho)])
	end
	
% 260   do 215 j=1,4
%       	do 210 i=1,3
% 210   		hp(i,j)=0.
%       	if (j.ne.1) hp(j-1,j)=eps
% 215   	yp(j)=r1(hp(1,j))
%       rmin=yp(1)
%       ftol=(.1)**(2*nsig+1)
%       iter=maxfn
%       write(6,*)
%       write(6,*) 'calling amoeba--r = ',rmin
%       rprev=rmin
%       call amoeba(hp,yp,3,3,ftol,r1,iter,work)
%       rmin=r1(hp(1,1))
%       call trans2(hp(1,1),qhati,qhat)
%       write(6,*) 'return from amoeba--r = ',rmin
%       write(6,*) 'r improvement: ',rprev-rmin
%       write(6,*) 'h: ',(hp(i,1),i=1,3)
%       write(6,*) 'qhat: ',qhat
%       call trans5(qhat,alat,along,rho)
%       write(6,*) 'Fitted rotation: ',alat,along,rho
%       write(6,*) 'ftol, iter: ',ftol,iter
%       write(6,*) 'do you want to call amoeba again? '
%       read(5,'(a)') ans
%       if ((ans.eq.'Y').or.(ans.eq.'y')) then
%           eps=eps/20.
%           do 212 i=1,4
% 212       qhati(i)=qhat(i)
%           go to 260
%       endif
%       call trans1(alat,along,axis)
%       call trans4(qhat,ahat)
%       write(6,*) 'Fitted rotation--alat,along,rho: '
%       write(6,*) alat,along,rho
%       write(6,*) 'ahat: '
%       do 211 i=1,3
% 211   write(6,*) (ahat(i,j),j=1,3)

% replace sigma-tilde with (ahat**t)*sigma-tilde*ahat

	temp = zeros(3,3);
	for (isect = 1:nsect)
		for (i = 1:3)
			for (j = 1:3)
				for (k1 = 1:3)
					for (k2 = 1:3)
						temp(i,j) = temp(i,j) + ahat(k1,i)*sigma(2,isect,k1,k2)*ahat(k2,j);
					end
				end
			end
		end
		sigma(2,isect,1:3,1:3) = temp;
	end

% calculation of (x**t)*x matrix--stored in sigma2 in symmetric mode
% eta is the matrix m(eta) of paper; etai is the matrix o-sub i.
% eta and etai are set by r1.

	sigma2 = zeros(msig,1);
	k = 0;
	for (i = 1:3)
		for (j = 1:i)
			k = k + 1;
			for (isect = 1:nsect)
				for (k1 = 1:3)
					for (k2 = 1:3)
						sigma2(k) = sigma2(k) + eta(isect,i,k1)*sigma(2,isect,k1,k2)*eta(isect,j,k2);
					end
				end
			end
		end
	end
	for (isect = 1:nsect)
		for (i = 1:2)
			for (j = 1:3)
				k = k + 1;
				for (k1 = 1:3)
					for (k2 = 1:3)
						sigma2(k) = sigma2(k) + etai(isect,k1,i)*sigma(2,isect,k1,k2)*eta(isect,j,k2);
					end
				end
			end
			k = k+2*(isect-1);
			for (j = 1:i)
				k = k + 1;
				for (k1 = 1:3)
					for (k2 = 1:3)
						sigma2(k) = sigma2(k) + etai(isect,k1,i)*(sigma(1,isect,k1,k2)+sigma(2,isect,k1,k2))*etai(isect,k2,j);
					end
				end
			end
		end
	end
     
	ndim = 2*nsect + 3;
	k = 0;
	for (i = 1:ndim)
		for (j = 1:i)
			k = k+1;
			sigma3(i,j) = sigma2(k);
			sigma3(j,i) = sigma2(k);
		end
	end
	

%       do 300 i=1,msig
% 300   sigma2(i)=0.
%       k=0
%       do 301 i=1,3
%       do 301 j=1,i
%       k=k+1
%       do 301 isect=1,nsect
%       do 301 k1=1,3
%       do 301 k2=1,3
% 301   sigma2(k)=sigma2(k)+
%      & eta(isect,i,k1)*sigma(2,isect,k1,k2)*eta(isect,j,k2)
%       do 302 isect=1,nsect
%       do 302 i=1,2
%       do 303 j=1,3
%       k=k+1
%       do 303 k1=1,3
%       do 303 k2=1,3
% 303   sigma2(k)=sigma2(k)+
%      & etai(isect,k1,i)*sigma(2,isect,k1,k2)*eta(isect,j,k2)
%       k=k+2*(isect-1)
%       do 302 j=1,i
%       k=k+1
%       do 302 k1=1,3
%       do 302 k2=1,3
% 302   sigma2(k)=sigma2(k)+etai(isect,k1,i)*(sigma(1,isect,k1,k2)+
%      & sigma(2,isect,k1,k2))*etai(isect,k2,j)
% 
%       ndim=2*nsect+3
%       k=0
%       do 304 i=1,ndim
%       do 304 j=1,i
%       k=k+1
%       sigma3(i,j)=sigma2(k)
% 304   sigma3(j,i)=sigma2(k)

% sigma3 is the full variance-covariance matrix of the rotation and
% normal section parameters (except possibly for division by hatkap).
% rotation parameters are stored in upper 3 by 3 corner.

	[V,D] = eig(sigma3);
	[D, ind] = sort(diag(D)');		% Make sure that eigs are in ascending order
	V = V(:,ind);
	sigma3 = zeros(ndim);
	for (i = 1:ndim)
		for (j = 1:ndim)
			for (k = 1:ndim)
				sigma3(i,j) = sigma3(i,j) + V(i,k)*V(j,k) / D(k);
			end
		end
	end

% calculation of critical point
	df = ndata - 2*nsect - 3;
	plev1 = (1-plev) / 2;
	[xchi1,ier] = xidch(plev1,df);
	plev1 = 1 - plev1;
	[xchi,ier] = xidch(plev1,df);
	fkap1 = xchi / rmin;
	fkap2 = xchi1 / rmin;
	hatkap = df / rmin;
	[xchi,ier] = xidf(plev, 3, df);
	xchi = xchi * rmin * 3 / df;

%       write(6,*)
%       write(6,*) 'Enter confidence level: '
%       read(5,*) plev
%       write(6,*)
%       write(6,*) 'Do you want an estimated kappa? '
%       read(5,1000) ans
%       if ((ans.eq.'n').or.(ans.eq.'n')) then
%          call xidch(plev,3.,xchi,ier)
%          if (ier.ne.0) write(6,*) 'subroutine xidch--ier ',ier
%       else
%          df=ndat-2*nsect-3
%          plev1=(1.-plev)/2
%          call xidch(plev1,df,xchi1,ier)
%          if (ier.ne.0) write(6,*) 'subroutine xidch--ier: ',ier
%          plev1=1-plev1
%          call xidch(plev1,df,xchi,ier)
%          if (ier.ne.0) write(6,*) 'subroutine xidch--ier: ',ier
%          fkap1=xchi/rmin
%          fkap2=xchi1/rmin
%          write(6,*) plev,' confidence interval for kappa'
%          write(6,*) fkap2,fkap1
%          hatkap=df/rmin
%          write(6,*) 'kappahat, degrees of freedom: ',hatkap, df
%          call xidf(plev,3.,df,xchi,ier)
%          if (ier.ne.0) write(6,*) 'subroutine xidf--ier: ',ier
%          xchi=xchi*rmin*3./df
%       endif

% calculation of h11.2 matrix--stored sigma4
	cov = sigma3(1:3,1:3);
	[V,D] = eig(sigma3);
	[D, ind] = sort(diag(D)');		% Make sure that eigs are in ascending order
	V = V(:,ind);
	sigma4 = zeros(3);
	for (i = 1:3)
		for (j = 1:3)
			for (k = 1:3)
				sigma4(i,j) = sigma4(i,j) + V(i,k)*V(j,k) / D(k);
			end
		end
	end
	D = 1 ./ D;
	[alat,along] = trans6(V(:,j));
	%angle = (D(j) / xchi)^(-.5);
	%angled = angle*45 / atan(1);

%       do 317 i=1,3
%       do 317 j=1,3
% 317   cov(i,j)=sigma3(i,j)
%       call jacobi(sigma3,3,msig2,d,z,msig2,work,nrot)
%       if (nrot.lt.0) write(6,*) 'subroutine jacobi(2)--nrot ',nrot
%       do 311 i=1,3
%       do 311 j=1,3
%       sigma4(i,j)=0.
%       do 311 k=1,3
% 311   sigma4(i,j)=sigma4(i,j)+z(i,k)*z(j,k)/d(k)
%       do 316 i=1,3
% 316   d(i)=1./d(i)
%       write(6,*)
%       write(6,*) 'h11.2: ',((sigma4(i,j),j=1,i),i=1,3)
%       write(6,*) 'eigenvalues: ',(d(j),j=1,3)
%       do 312 j=1,3
% 312   write(6,*) 'eigenvector: ',(z(i,j),i=1,3)
%       write(6,*)
%       write(6,*) '"confidence region endmembers"'
%       do 313 j=1,3
%       write(6,*) 'u-side eigenvector: ',(sngl(z(i,j)),i=1,3)
%       call trans6(z(1,j),alat,along)
%       write(6,*) 'latitude, longitude: ',sngl(alat),sngl(along)
%       angle=(d(j)/xchi)**(-.5)
%       angled=angle*45./atan(dble(1.))
%       j1=2*j
%       j2=j1+1
%       do 314 i=1,3
%       h(i)=0.
%       do 314 k=1,3
%       h(i)=h(i)+ahat(i,k)*z(k,j)
% 314   continue
%       write(6,*) 'v-side eigenvector: ',(sngl(h(i)),i=1,3)
%       call trans6(h,alat,along)
%       write(6,*) 'latitude, longitude: ',sngl(alat),sngl(along)
% 313   write(6,*) 'angle of rotation: radians, degrees ', sngl(angle),sngl(angled)

	[alat,along,rho] = trans5(qhat);

%       write(15,1460) filnam
%  1460 format ('Results from Hellinger1 using ',a,/)
%       write(15,*) 'Fitted rotation--alat,along,rho: '
%       call trans5(qhat,alat,along,rho)
%       write(15,*) alat,along,rho
%       write(15,*) 'conf. level, conf. interval for kappa: '
%       write(15,*) plev,fkap2,fkap1
%       if(hatkap.eq.0.) hatkap=1.
%       if (df.eq.0.) df=10000.
% c  flag is used by addplus.f to indicate origin of estimates
%       write(15,*) 'kappahat, degrees of freedom,xchi,flag'
%       write(15,*) hatkap,df,xchi,0.0
%       write(15,*) 'Number of points, sections'
%       write(15,*) ndat,nsect
%       write(15,*) 'ahat: '
%       do 460 i=1,3
% 460   write(15,*) (ahat(i,j),j=1,3)
%       write(15,*) 'covariance matrix'
%       do 465 i=1,3
% 465   write(15,*) (cov(i,j),j=1,3)
%       write(15,*) 'H11.2 matrix: '
%       do 470 i=1,3
% 470   write(15,*) (sigma4(i,j),j=1,3)

% ---------------------------------------------------------------------------
function [r, eta, etai] = r1(h, sigma, qhati, nsect)

	global eta etai
	qhat = trans2(h,qhati);
	ahat = trans4(qhat);
	r = 0;
	sig = zeros(3,3);
	% We will create 3D variables because the 4D sigma does not JIT
	sigma_1 = sigma(1,:,:,:);
	sigma_2 = sigma(2,:,:,:);
	siz = size(sigma_2);
	sigma_1 = reshape(sigma_1, siz(2:4));
	sigma_2 = reshape(sigma_2, siz(2:4));

	for (i = 1:nsect)
		for (j = 1:3)
			for (k = 1:3)
				sig(j,k) = sigma_1(i, j, k);
				for (k1 = 1:3)
					for (k2 = 1:3)
						sig(j,k) = sig(j,k) + ahat(k1,j)*sigma_2(i,k1,k2)*ahat(k2,k);
					end
				end
			end
		end
		[V,D] = eig(sig);
		[D, ind] = sort(diag(D));
		V = V(:,ind);
		eta(i,2,1) = V(3,1);
		eta(i,3,1) = -V(2,1);
		eta(i,1,2) = -V(3,1);
		eta(i,3,2) = V(1,1);
		eta(i,1,3) = V(2,1);
		eta(i,2,3) = -V(1,1);
		if (abs(V(3,1)) > .2)
			etai(i,1:3,1) = eta(i,1:3,1);
			etai(i,1:3,2) = eta(i,1:3,2);
		elseif (abs(V(1,1)) > .2)
			etai(i,1:3,1) = eta(i,1:3,2);
			etai(i,1:3,2) = eta(i,1:3,3);
		else
			etai(i,1:3,1) = eta(i,1:3,1);
			etai(i,1:3,2) = eta(i,1:3,3);
		end
		r = r + D(1);
	end

% -------------------------------------------------------------------
function [qhati, eps_, rmin] = grds(eps_, sigma, qhati, nsect)
	rmin = 1.1e30;
	eps5 = eps_ / 5;
	for (i = -5:5)
		for (j = -5:5)
			for (k = -5:5)
				h = [i; j; k] * eps5;
				r = r1(h, sigma, qhati, nsect);
				if (r < rmin)
					imin = i;
					jmin = j;
					kmin = k;
					rmin = r;
				end
			end
		end
	end

	h = [imin; jmin; kmin] * eps5;
	qhati = trans2(h,qhati);
	eps_ = eps5;	

% -------------------------------------------------------------------
function ahat = trans3(alat,along,rho)
% TRANSLATES AXIS LATITUDE AND LONGITUDE, ANGLE OF ROTATION TO A ROTATION EXPRESSED AS A QUARTERNION
	ahat = zeros(4,1);
	D2R = pi / 180;
	ahat(1)=cos(rho*D2R/2.);
	fact=sin(rho*D2R/2.);
	ahat(4)=fact*sin(alat*D2R);
	fact=fact*cos(alat*D2R);
	ahat(2)=fact*cos(along*D2R);
	ahat(3)=fact*sin(along*D2R);
	
% -------------------------------------------------------------------
function xhat = trans2(x,ahat)
% SUBROUTINE TO EXPRESS XHAT=AHAT*EXP(X)--XHAT AND AHAT EXPRESSED AS QUARTERNIONS.
	theta = sqrt(x(1)^2 + x(2)^2 + x(3)^2);
	if (theta <= 0)
		xhat = ahat(1:4);
		return
	end
	expx = zeros(4,1);
	expx(1) = cos(theta/2);
	fact = sin(theta/2)/theta;
	expx(2:4) = fact * x(1:3);
	xhat = qmult(ahat,expx);

% -------------------------------------------------------------------
function [U1,U2,U3] = trans1(ULAT,ULONG)
% SUBROUTINE TO TRANSLATE LATITUDE AND LONGITUDE TO EUCLIDEAN COORDINATES
	ULAT = ULAT *pi/180;
	ULONG = ULONG *pi/180;
	U3=sin(ULAT);
	U1=cos(ULAT) .* cos(ULONG);
	U2=cos(ULAT) .* sin(ULONG);
% -------------------------------------------------------------------
function amhat = trans4(ahat)
% SUBROUTINE TO TRANSLATE A QUARTERNION INTO A ROTATION MATRIX
	amhat = zeros(3,3);
	A0=ahat(1);
	A1=ahat(2);
	A2=ahat(3);
	A3=ahat(4);
	amhat(1,1)=A0*A0+A1*A1-A2*A2-A3*A3;
	amhat(2,1)=2.*(A0*A3+A1*A2);
	amhat(3,1)=2.*(A1*A3-A0*A2);
	amhat(1,2)=2.*(A1*A2-A0*A3);
	amhat(2,2)=A0*A0-A1*A1+A2*A2-A3*A3;
	amhat(3,2)=2.*(A0*A1+A2*A3);
	amhat(1,3)=2.*(A0*A2+A1*A3);
	amhat(2,3)=2.*(A2*A3-A0*A1);
	amhat(3,3)=A0*A0-A1*A1-A2*A2+A3*A3;

% -------------------------------------------------------------------
function C = qmult(A,B)
	C = zeros(4,1);
	C(1)=A(1)*B(1)-A(2)*B(2)-A(3)*B(3)-A(4)*B(4);
	C(2)=A(1)*B(2)+A(2)*B(1)+A(3)*B(4)-A(4)*B(3);
	C(3)=A(1)*B(3)+A(3)*B(1)+A(4)*B(2)-A(2)*B(4);
	C(4)=A(1)*B(4)+A(4)*B(1)+A(2)*B(3)-A(3)*B(2);

% -------------------------------------------------------------------
function [alat, along, rho] = trans5(ahat)
% SUBROUTINE TO TRANSLATE A QUARTERNION TO A AXIS LATITUDE, LONGITUDE AND ANGLE
	R2D = 180 / pi;
	tmp = ahat(1);
	if (tmp > 1),		tmp = 1;	end
	if (tmp < -1),		tmp =-1;	end
	fact = acos(tmp);
	rho = fact * 2 * R2D;
	tmp=ahat(4)/ sin(fact);
	if (tmp > 1),		tmp=1;		end
	if (tmp < -1),		tmp=-1;		end
	alat  = asin(tmp) * R2D;
	along = atan2(ahat(3),ahat(2)) * R2D;

% -------------------------------------------------------------------
function [ULAT,ULONG] = trans6(U)
% SUBROUTINE TO TRANSLATE EUCLIDEAN COORDINATES INTO LATITUDE AND LONGITUDE
	R2D = 180 / pi;
	TEMP=U(3);
	if (TEMP > 1),		TEMP=1;		end
	if (TEMP < (-1)),	TEMP=-1;	end
	ULAT=asin(TEMP)* R2D;
	if ((U(1)^2 + U(2)^2) <= 0)
		ULONG=0;
		return
	end
	ULONG=atan2(U(2),U(1))* R2D;

% ---------------------------------------------------------------------
function [P, Y] = amoeba(P,Y,MP,NDIM,FTOL,ITER, funk, varargin)
% MODIFICATION OF PRESS ET AL: NUMBERICAL RECIPES, PAGES 292-293
% SIMPLEX METHOD OF NELDER AND MEAD
%
% CALLING SEQUENCE:
%     P--MP BY (NDIM+1) MATRIX WHOSE COLUMNS ARE LINEARLY INDEPENDENT
%        AND CLOSE TO INITAL GUESS.  ON OUTPUT, THE SIMPLEX SPANNED BY
%        THE COLUMNS OF P CONTAINS THE MINIMUM FOUND BY AMOEBA.
%     Y--VECTOR OF LENGTH (NDIM+1) CONTAINING THE VALUES OF THE FUNCTION
%        AT THE COLUMNS OF P.
%     MP--ROW DIMENSION OF P AS SPECIFIED IN CALLING PROGRAM'S DIMENSION
%         STATEMENT
%     NDIM--NUMBER OF PARAMETERS TO BE MINIMIZED
%     FTOL--INPUT: FUNCTION VALUES AT THE COLUMNS OF P SHOULD BE WITHIN
%                  FTOL OF MINIMUM VALUE
%           OUTPUT: ESTIMATED ACHIEVED VALUE OF FTOL
%     FUNK--NAME OF FUNCTION TO BE EVALUATED
%     ITER--INPUT: MAXIMUM NUMBER OF ITERATIONS.
%           OUTPUT: ACTUAL NUMBER OF ITERATIONS.
%     WORK--WORK VECTOR OF LENGTH 3*NDIM
%
	MPTS=NDIM+1;
	ITMAX=ITER;
	ALPHA=1;		BETA=0.5;	GAMMA=2.0;
	WORK = zeros(NDIM,3);

	for i=1:MPTS
		Y(i) = feval(funk,P(:,i), varargin{:});
	end
	ITER=0;
	while (true)
		ILO=1;
		if (Y(1) > Y(2))
			IHI=1;		INHI=2;
		else
			IHI=2;		INHI=1;
		end
		for (i=1:MPTS)
			if (Y(i) < Y(ILO)),		ILO=i;	end
			if (Y(i) > Y(IHI))
				INHI=IHI;
				IHI=i;
			elseif (Y(i) > Y(INHI))
				if (i ~= IHI),	INHI=i;		end
			end
		end
		RTOL = 2 * abs(Y(IHI)-Y(ILO))/(abs(Y(IHI))+abs(Y(ILO)));
		if ((RTOL < FTOL) || (ITER >= ITMAX))
			FTOL=RTOL;
			if (ILO == 1),	return,		end
			for (i = 1:NDIM)
				TEMP=P(i,1);
				P(i,1)=P(i,ILO);
				P(i,ILO)=TEMP;
			end
			TEMP=Y(1);
			Y(1)=Y(ILO);
			Y(ILO)=TEMP;
			return
		end
	
		ITER=ITER+1;
		WORK(1:NDIM,1) = 0;
	
		for (i = 1:MPTS)
			if (i ~= IHI)
				WORK(1:NDIM,1) = WORK(1:NDIM,1) + P(1:NDIM,i);
			end
		end
	
		WORK(1:NDIM,1) = WORK(1:NDIM,1)/NDIM;
		WORK(1:NDIM,2) = (1+ALPHA)*WORK(1:NDIM,1) - ALPHA*P(1:NDIM,IHI);
	
		YPR = feval(funk,WORK(:,2), varargin{:});
		if (YPR <= Y(ILO))
			WORK(1:NDIM,3) = GAMMA*WORK(1:NDIM,2) + (1-GAMMA)*WORK(1:NDIM,1);
			YPRR = feval(funk, WORK(:,3), varargin{:});
			if (YPRR < Y(ILO))
				P(1:NDIM,IHI)=WORK(1:NDIM,3);
          		Y(IHI)=YPRR;
			else
				P(1:NDIM,IHI)=WORK(1:NDIM,2);
          		Y(IHI)=YPR;
			end
		elseif (YPR >= Y(INHI))
			if (YPR < Y(IHI))
				P(1:NDIM,IHI)=WORK(1:NDIM,2);
          		Y(IHI)=YPR;
			end
			WORK(1:NDIM,3)=BETA*P(1:NDIM,IHI)+(1-BETA)*WORK(1:NDIM,1);
          	YPRR=feval(funk,WORK(:,3), varargin{:});
			if (YPRR < Y(IHI))
				P(1:NDIM,IHI)=WORK(1:NDIM,3);
				Y(IHI)=YPRR;
			else
				for (i = 1:MPTS)
					if (i ~= ILO)
						WORK(1:NDIM,2)=0.5*(P(1:NDIM,i)+P(1:NDIM,ILO));
						P(1:NDIM,i)=WORK(1:NDIM,2);
						Y(i)=feval(funk,WORK(:,2), varargin{:});
					end
				end
			end
		else
			P(1:NDIM,IHI)=WORK(1:NDIM,2);
			Y(IHI)=YPR;
		end
	end

% -------------------------------------------------------------------------------------
function PLEV = xdf(D,XF,IER)
	X = D(2)/(D(2)+D(1)*XF);
	PLEV = 1 - betai(D(2)/2.,D(1)/2.,X,IER);

function PLEV = xdn(DUMMY,XN,IER)
	X= XN / 1.4142136;
	PLEV = .5 + .5 * erf0(X,IER);

function PLEV = xdch(DF,XCHI,IER)
	PLEV = gammp(DF(1)/2.,XCHI/2.,IER);

% -------------------------------------------------------------------------------------
function ZBR = zbrent(FUNC,FVAL,DUMMY,X1,X2,TOL,NER)
%       Root finding using Brent algorithms
%       Modification of Press et al, Numerical Recipes, pages 253-254

	ITMAX = 100;	eps_ = 3.E-8;

	A=X1;	B=X2;	NER=0;	IER=0;
	FA = feval(FUNC,DUMMY,A,IER);
	FA=FA-FVAL;
	if (IER ~= 0),	NER=NER+1;		end
	IER=0;
	FB = feval(FUNC,DUMMY,B,IER);
	FB=FB-FVAL;
	if (IER ~= 0),	NER=NER+1;		end
	if (FB*FA > 0)
		NER=-NER-1;
		ZBR = 0;
		return
	end
	FC=FB;
	for (ITER = 1:ITMAX)
		if (FB*FC > 0)
			C=A;	FC=FA;		D=B-A;		E=D;
		end
		if (abs(FC) < abs(FB))
			A = B;		B = C;		C = A;		FA = FB;	FB = FC;	FC = FA;
		end
		TOL1 = 2*eps_*abs(B)+0.5*TOL;
		XM = .5*(C-B);
		if (abs(XM) <= TOL1 || FB == 0)
            ZBR = B;
            return
		end
		if (abs(E) >= TOL1 && abs(FA) > abs(FB))
			S=FB/FA;
			if (A == C)
				P=2.*XM*S;
				Q=1.-S;
			else
				Q=FA/FC;
				R=FB/FC;
				P=S*(2.*XM*Q*(Q-R)-(B-A)*(R-1.));
				Q=(Q-1.)*(R-1.)*(S-1.);
			end
			if (P > 0),	Q=-Q;	end
			P=abs(P);
			if (2.*P < min(3.*XM*Q-abs(TOL1*Q),abs(E*Q)))
				E=D;
				D=P/Q;
			else
				D=XM;
				E=D;
			end
		else
			D=XM;
			E=D;
		end
		A=B;
		FA=FB;
		if (abs(D) > TOL1)
			B=B+D;
		else
			if (XM >= 0),	B = B + abs(TOL1);
			else			B = B - abs(TOL1);
			end
			%B=B+SIGN(TOL1,XM);
		end
		IER=0;
		FB = feval(FUNC,DUMMY,B,IER);
		FB=FB-FVAL;
		if (IER ~= 0),	NER=NER+1;	end
	end
	NER = -1000*NER;
	ZBR = B;

% -------------------------------------------------------------------------------------
% INVERSE DISTRIBUTION ROUTINES

function [XN,IER] = xidn(PLEV)
	TOL = 0.0001;
	P = [.5 .6 .7 .8 .85 .9 .95 .975 .9875 .99 .995 .9975 .999 .9995];
	TABLE = [0. .253 .524 .842 1.036 1.282 1.645 1.960 2.240 2.326 2.576 2.807 3.090 3.291];

	IER = 0;
	QLEV = .5 + abs(PLEV-.5);
	from_a_goto = false;
	for (I=1:14)
		if (abs(QLEV-P(I)) < .00049)
			XN = TABLE(I);  
			from_a_goto = true;
			break
		end
		if (QLEV < P(I))
			XN = zbrent('xdn',QLEV,0,TABLE(I-1),TABLE(I),TOL,IER);
			from_a_goto = true;
			break
		end
	end
	if (~from_a_goto)
		XN = TABLE(14);
		IER = 111111;
	end
	if (PLEV >= 0.5),	return,		end
	XN = -XN;

% -------------------------------------------------------------------------------------
function [XCHI,IER] = xidch(PLEV,DF)
      
	TOL = 0.0001;
	P = [.5 .7 .9 .95 .99];
	IER = 0;
	if (DF <= 0)
		XCHI=0.;
		IER=111111;
	end
      
	D = DF;
	IDF = DF + .001;
	if ((abs(IDF-DF) < .001) && (IDF <= 22))
		for (I=1:5)
			if (abs(PLEV-P(I)) < .00049)
				XCHI = x2tab(PLEV,IDF);
				return
			end
			if (I ~= 1)
				if ((PLEV > P(I-1)) && (PLEV < P(I)))
					A = x2tab(P(I-1),IDF);
					B = x2tab(P(I),IDF);
					XCHI = zbrent('xdch',PLEV,D,A,B,TOL,IER);
					return
				end
			end
		end
	end
	[B,IER] = xidn(PLEV);
	B = DF + sqrt(2.*DF)*B;
	if (IER ~= 0),	B=DF;	end
	B = max(1,B);
	A = B;
	PA = xdch(D,A,IER);
	while (PA >= PLEV)
		A = .9 * A;
		PA = xdch(D,A,IER);
	end
	PB = xdch(D,B,IER);
	while (PB <= PLEV)
		B = 1.1 * B;
		PB = xdch(D,B,IER);
	end
	XCHI = zbrent('xdch',PLEV,D,A,B,TOL,IER);

% -------------------------------------------------------------------------------------
function [XF,IER] = xidf(PLEV,D1,D2)
	TOL = [0.0001];
	D = [D1 D2];
	if ((PLEV <= 0) || (PLEV >= 1) || (D1 <= 0) || (D2 <= 0))
		XF = 1;
		IER=111111;
		return
	end
	[XF,IER] = xidch(PLEV,D(1));
	A = XF;
	PA = xdf(D,A,IER);
	while (PA <= PLEV)
		A = A / 2;
		PA = xdf(D,A,IER);
	end
	B = XF;
	PB = xdf(D,B,IER);
	while (PB <= PLEV)
		B = 2 * B;
		PB = xdf(D,B,IER);
	end
	XF = zbrent('xdf',PLEV,D,A,B,TOL,IER);

% -------------------------------------------------------------------------------------
function X2VAL = x2tab(PC,IDF)
	CENT = [.5 .7 .9 .95 .99];
	TABLE = [.455,1.386,2.366,3.357,4.351,5.348,6.346,7.344,8.343, ...
		9.342,10.341,11.34,12.34,13.339,14.339,15.338,16.338, ...
		17.338,18.338,19.337,20.337,21.337, ...
		1.074,2.408,3.665,4.878,6.064,7.231,8.383, ...
		9.524,10.656,11.781,12.899,14.011,15.119,16.222, ...
		17.322,18.418,19.511,20.601,21.689,22.775,23.858,24.939, ...
		2.706,4.605,6.251,7.779,9.236,10.645,12.017,13.362, ...
		14.684,15.987,17.275,18.549,19.812,21.064,22.307, ...
		23.542,24.769,25.989,27.204,28.412,29.615,30.813, ...
		3.841,5.991,7.815,9.488,11.07,12.592,14.067,15.507, ...
		16.919,18.307,19.675,21.026,22.362,23.685,24.996, ...
		26.296,27.587,28.869,30.144,31.41,32.671,33.924, ...
		6.635,9.21,11.341,13.277,15.086,16.812,18.475, ...
		20.09,21.666,23.209,24.725,26.217,27.688,29.141, ...
		30.578,32.,33.409,34.805,36.191,37.566,38.932,40.289];

	TABLE = reshape(TABLE,22,5);
	for (J=1:5)
		if (abs(CENT(J)-PC) < .0005),	break,		end
	end
	X2VAL = TABLE(fix(IDF),J);

% -------------------------------------------------------------------------------------
% C***********************************************************************
% C                    SPECIAL FUNCTION ROUTINES
% C  MODIFIED FROM PRESS ET AL: NUMERICAL RECIPES
% C***********************************************************************
function g = gammp(A,X,IER)
% INCOMPLETE GAMMA FUNCTION
	if (X < 0 || A <= 0)
         IER=100;
         g = 0;
         return
	end
	GLN = 0;
	if (X < A+1.)
		gammpT = gser(A,X,GLN,IER);
		g = gammpT;
	else
		GAMMCF = gcf(A,X,GLN,IER);
		g = 1 - GAMMCF;
	end

% -------------------------------------------------------------------------------------
function GAMSER = gser(A,X,GLN,IER)
% INCOMPLETE GAMMA FUNCTION USING SERIES REPRESENTATION; GLN=LN GAMMA(A)
	ITMAX = 100;	eps_ = 3.E-7;
	GLN=gammln(A);
	if (X <= 0)
		if (X < 0),	IER=100;	end
		GAMSER=0;
		return
	end
	AP=A;
	SOM=1./A;
	DEL=SOM;
	IER=101;
	for (N=1:ITMAX)
		AP=AP+1;
		DEL=DEL*X/AP;
		SOM=SOM+DEL;
		if (abs(DEL) < abs(SOM)*eps_),	IER = 100;	break,	end
	end
	GAMSER=SOM*exp(-X+A*log(X)-GLN);

% -------------------------------------------------------------------------------------
function GAMMCF = gcf(A,X,GLN,IER)
% CONTINUED FRACTION CALCULATION OF INCOMPLETE GAMMA FUNCTION
	ITMAX = 100;	eps_ = 3.E-7;
	GLN=gammln(A);
	GOLD=0;
	A0=1;
	A1=X;
	B0=0;
	B1=1;
	FAC=1;
	for (N=1:ITMAX)
		AN=N;
		ANA=AN-A;
		A0=(A1+A0*ANA)*FAC;
		B0=(B1+B0*ANA)*FAC;
		ANF=AN*FAC;
		A1=X*A0+ANF*A1;
		B1=X*B0+ANF*B1;
		if (A1 ~= 0)
			FAC=1./A1;
			G=B1*FAC;
			if (abs((G-GOLD)/G) < eps_),		break,	end
			GOLD=G;
		end
	end
	%IER=102;
	GAMMCF = exp(-X+A*log(X)-GLN)*G;

% -------------------------------------------------------------------------------------
function G = gammln(XX)
	COF = [76.18009173 -86.50532033 24.01409822 -1.231739516 .120858003e-2 -.536382D-5];
	STP = 2.50662827465;
	FPF = 5.5;
	ONE = 1;
	HALF = 0.5;
	if (XX == 1)
		G = 0;
		return
	end
	X=abs(XX-ONE);
	Z=X;
	TMP=X+FPF;
	TMP=(X+HALF)*log(TMP)-TMP;
	SER=ONE;
	for (J=1:6)
		X=X+ONE;
		SER=SER+COF(J)/X;
	end
	DGAMM=TMP+log(STP*SER);
	if (XX > 1)
		G = DGAMM;
	else
		G = log(pi*Z/sin(pi*Z))-DGAMM;
	end

% -------------------------------------------------------------------------------------
function E = erf0(X,IER)
	if (X < 0)
		E = -gammp(.5,X^2,IER);
	else
		E = gammp(.5,X^2,IER);
	end

% -------------------------------------------------------------------------------------
function B = betai(A,B,X,IER)
	if (X <= 0)
		B=0;
		if (X < 0),		IER=105;	end
	elseif (X >= 1)
		B=1;
		if (X > 1),		IER=105;	end
	else
		BT = exp(gammln(A+B)-gammln(A)-gammln(B)+A*log(X)+B*log(1-X));
		if (X < ((A+1.)/(A+B+2.)))
			B = BT*betacf(A,B,X,IER)/A;
		else
			B = 1 - BT*betacf(B,A,1.-X,IER)/B;
		end
	end

% -------------------------------------------------------------------------------------
function B = betacf(A,B,X,IER)
	ITMAX = 100;	eps_ = 3.E-7;
	AM=1;
	BM=1;
	AZ=1;
	QAB=A+B;
	QAP=A+1;
	QAM=A-1;
	BZ=1.-QAB*X/QAP;
	for (M=1:ITMAX)
		EM=M;
		TEM=EM+EM;
		D=EM*(B-M)*X/((QAM+TEM)*(A+TEM));
		AP=AZ+D*AM;
		BP=BZ+D*BM;
		D=-(A+EM)*(QAB+EM)*X/((A+TEM)*(QAP+TEM));
		APP=AP+D*AZ;
		BPP=BP+D*BZ;
		AOLD=AZ;
		AM=AP/BPP;
		BM=BP/BPP;
		AZ=APP/BPP;
		BZ=1;
		if (abs(AZ-AOLD) < eps_*abs(AZ)), 	break,		end
	end
	IER = 106;
	B = AZ;

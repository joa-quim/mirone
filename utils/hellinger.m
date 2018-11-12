function [along, alat, rho, vol, t_stats, ellip_long, ellip_lat] = ...
	      hellinger(along, alat, rho, isoc_mov, isoc_fix, DP_tol, force_pole, show_segs, isoc1_props, isoc2_props)
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

% The original fortran code had no License so I'm releasing this translation, 
% that is part of Mirone, to the Public Domain
%
%	Copyright (c) 2004-2017 by J. Luis
%
%	Contact info: w3.ualg.pt/~jluis/mirone
% --------------------------------------------------------------------

% $Id: hellinger.m 11303 2018-05-28 21:39:31Z Joaquim Luis $

	do_geodetic = true;
	plot_segmentation_lines = true;
	std_invented = 4;		% If this is ~= 0 the residues transmitted in the ISOC1_PROPS (non pure Hellinger) are ignored

	plev = 0.95;				% Confidence level
	rfact = 4.0528473e07;
	hlim=1e-10;
	nsig = 4;
	maxfn = 1000;
 	D2R = pi / 180;
	eps_ = 10;					% If == 0, ouput pole == input
	if (force_pole),	eps_ = 0;	end

	% Converting to geocentric, as is correct, but makes little difference
	if (do_geodetic)
		ecc = 0.0818191908426215;			% WGS84
		isoc_mov = isoc_mov * D2R;
		isoc_fix = isoc_fix * D2R;
		isoc_mov(:,2) = atan2((1-ecc^2)*sin(isoc_mov(:,2)), cos(isoc_mov(:,2)));	% Lat da isoc_mov geocentrica
		isoc_fix(:,2) = atan2((1-ecc^2)*sin(isoc_fix(:,2)), cos(isoc_fix(:,2)));	% Lat da isoc_fix geocentrica
		isoc_mov = isoc_mov / D2R;
		isoc_fix = isoc_fix / D2R;
	end

	if (nargin == 10)
		data = [ones(size(isoc_mov,1),1) isoc1_props(:,1) isoc_mov(:,2:-1:1) isoc1_props(:,2); ...
			ones(size(isoc_fix,1),1)*2 isoc2_props(:,1) isoc_fix(:,2:-1:1) isoc2_props(:,2)];
	else
		residues = isoc1_props;		% Unfortunatelly this doesn't work as I thought. Solution gets unstable
 		handles = guidata(gcf);
		
		% ------------------------------- Guess Hellinger segmentation ------------------------------------
		% Apply a line simplification to estimate the Hellinger segmentation
		B = cvlib_mex('dp', isoc_mov, DP_tol, 'GEOG');
		lat_dp = B(:,2);	lon_dp = B(:,1);
		[c,ind] = intersect(isoc_mov(:,2),lat_dp);	% Get only those points that are common to both data-sets
		ind = sort(ind);	% We need them in original order and the 'stable' flag to intersect doesn't exist in R13

% 		breaks = (diff(ind) == 1);				% This tells us whose segments have no middle points
% 		breaks(1) = false;						% Otherwise we would be loosing the first point
% 		for (k = 2:numel(breaks)-1)				% And do not let break two segments in a row (would loose points)
% 			if (breaks(k) && breaks(k+1))
% 				breaks(k) = false;
% 			end
% 		end
		segmented_dp = zeros(3*(numel(ind) - 1), 2);
		n = 1;
		for (k = 1:numel(ind)-1)
			segmented_dp(n,:) = [lon_dp(k) lat_dp(k)];		n = n + 1;
			segmented_dp(n,:) = [lon_dp(k+1) lat_dp(k+1)];	n = n + 1;
			segmented_dp(n,:) = [NaN NaN];
% 			if (breaks(k))
% 				segmented_dp(n-2:n,1) = Inf;	% Use Inf to later be able to find these and delete them
% 			end
			n = n + 1;
		end
		ind = isinf(segmented_dp(:,1));			% Remove the segments that have no middle points in them
		segmented_dp(ind,:) = [];

		% Compute the segment indices to which each point in isoc_mov belongs 
		flags_mov = segmentate(isoc_mov, segmented_dp);

		% Now rotate the segmented_dp line and do the same for isoc_fix
		[r_lon,r_lat] = rot_euler(segmented_dp(:,1),segmented_dp(:,2), along, alat, rho, -1);		% Rotate DP moving isoc
		segmented_dp_rot = [r_lon r_lat];
		flags_fix = segmentate(isoc_fix, segmented_dp_rot);

		% I found cases where covariance was NaNs when one of the isocs had more than 2 segment at the end than the other
		% example flags_fix = [... 24 25]; flags_mov = [... 24 25 26 27 28];
		% This trick avoided that case and we should have lost no valid points because extreme points were unpaired anyway.
		delta = flags_fix(end) - flags_mov(end);
		if (abs(delta) > 1)
			k = 0;
			if (flags_fix(end) > flags_mov(end))
				while (flags_mov(end-k) > flags_fix(end) + 1)		% delta > 0
					flags_fix(end-k) = flags_mov(end) + 1;			k = k + 1;
				end
			else													% delta < 0
				while (flags_mov(end-k) > flags_fix(end) + 1)
					flags_mov(end-k) = flags_fix(end) + 1;			k = k + 1;
				end
			end
		end

		if (plot_segmentation_lines)			% Plot the segmentation guess as well as its rotation
			lat_geod = segmented_dp(:,2);
			lat_geod_rot = r_lat;
			if (do_geodetic)					% Convert back to geodetic latitudes
				lat_geod = lat_geod * D2R;
				lat_geod = atan2(sin(lat_geod), (1-ecc^2)*cos(lat_geod)) / D2R;
				lat_geod_rot = lat_geod_rot * D2R;
				lat_geod_rot = atan2(sin(lat_geod_rot), (1-ecc^2)*cos(lat_geod_rot)) / D2R;
			end
			
			hLine = findobj(get(handles.hCallingFig,'CurrentAxes'), 'Tag','broken_dp');
			if (isempty(hLine))
				hLine = line('parent',get(handles.hCallingFig,'CurrentAxes'),'XData',segmented_dp(:,1),'YData',lat_geod, ...
							 'LineStyle','-','LineWidth',1,'Color','b','Tag','broken_dp');
				draw_funs(hLine,'line_uicontext')
				hLine = line('parent',get(handles.hCallingFig,'CurrentAxes'),'XData',r_lon, 'YData',lat_geod_rot, ...
							 'LineStyle','-','LineWidth',1,'Color','b','Tag','broken_dp_rot');
				draw_funs(hLine,'line_uicontext')
			else
				set(hLine, 'XData', segmented_dp(:,1), 'YData',lat_geod)
				hLine = findobj(get(handles.hCallingFig,'CurrentAxes'), 'Tag','broken_dp_rot');
				set(hLine, 'XData',r_lon, 'YData',lat_geod_rot)
			end
		end
		% -------------------------------------------------------------------------------------------------------------

		if (std_invented)
			sig_mov = ones(size(isoc_mov,1),1)*std_invented;
			sig_fix = ones(size(isoc_fix,1),1)*std_invented;
		else
			residues(residues == 0) = 0.1;		% I simply do not believe them
			sig_mov = residues(1:size(isoc_mov,1));
			sig_fix = residues(size(isoc_mov,1)+1:end);
		end
		data = [ones(size(isoc_mov,1),1) flags_mov isoc_mov(:,2:-1:1)   sig_mov; ...
		        ones(size(isoc_fix,1),1)*2 flags_fix isoc_fix(:,2:-1:1) sig_fix];
		clear sig_fix sig_mov

		mixed = classical_hellinger_order(flags_mov, isoc_mov, flags_fix, isoc_fix);
		lat_geod = mixed(:,2);
		data_geod = data;					% For segmentation patches we need geodetic latitudes
		if (do_geodetic)					% Convert back to geodetic latitudes
			lat_geod = lat_geod * D2R;
			lat_geod = atan2(sin(lat_geod), (1-ecc^2)*cos(lat_geod)) / D2R;
			% Set the data column in the 'data_geod' array to geodetic coords (those we use in plots)
			lat_geod_or = data(:,3);
			lat_geod_or = lat_geod_or * D2R;
			lat_geod_or = atan2(sin(lat_geod_or), (1-ecc^2)*cos(lat_geod_or)) / D2R;
			data_geod(:,3) = lat_geod_or;
		end
		hLine = findobj(get(handles.hCallingFig,'CurrentAxes'), 'Tag','HellingerPicks');
		if (isempty(hLine))
			hLine = line('parent',get(handles.hCallingFig,'CurrentAxes'),'XData',mixed(:,1),'YData',lat_geod, ...
						 'LineStyle','-.','LineWidth',1,'Tag','HellingerPicks');
			draw_funs(hLine,'line_uicontext')			% Set lines's uicontextmenu
			set_extra_uicb_options(hLine, data_geod)
		else		% Just update old one
			set(hLine, 'XData', mixed(:,1),'YData',lat_geod)
		end
		update_patches(hLine, data_geod, show_segs)		% Either update or create first time creation
	end

	ndata = size(data,1);
	nsect = max(data(:,2));
	sigma = zeros(2,nsect,3,3);
	msig  = 2*nsect*nsect+7*nsect+6;
	%msig2 = 2*nsect+3;
	eta = zeros(nsect,3,3);
	etai= zeros(nsect,3,2);
	  
	[x,y,z] = trans1(data(:,3), data(:,4));
	axis = [x y z];
	sd = data(:,5).^2;

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
	h = ones(3,1);
	[rmin, eta, etai] = r1(h, sigma, qhati, nsect, eta, etai);		% Not sure but think it's for initializing eta & etai

	%[qhati, eps_, rmin] = grds(eps_, sigma, qhati, nsect);
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

	eps_ = eps_ * D2R;

	h = zeros(1,3);
	[rmin, eta, etai] = r1(h, sigma, qhati, nsect, eta, etai);
	rmin = rmin * rfact;

	% minimization section--calling amoeba
	hp = zeros(3,4);	yp = zeros(1,4);
	kflag = 0;
	count = 0;
	while (count < 15)
		for (j = 1:4)
			hp(1:3, j) = 0;
			if (j ~= 1)
				 hp(j-1,j) = eps_;
			end
			[yp(j), eta, etai] = r1(hp(1:end,j:end), sigma, qhati, nsect, eta, etai);
		end

		ftol = (0.1)^(2*nsig+1);
		iter = maxfn;
		[hp, yp, eta, etai] = amoeba(hp,yp,3,3,ftol,iter, 'r1', sigma, qhati, nsect, eta, etai);
		[rmin, eta, etai] = r1(hp, sigma, qhati, nsect, eta, etai);
		rmin = rmin * rfact;
		qhat=trans2(hp(:,1),qhati);
		[alat,along,rho] = trans5(qhat);
		eps_=eps_ / 20;
		qhati(1:4)=qhat(1:4);

		% run amoeba as long as hp(1:3,1) > 10E-10,  and one more time
		if (kflag == 1),	break,	end		% Break (get out of) the loop here
		if (abs(hp(1,1)) < hlim && abs(hp(2,1)) < hlim && abs(hp(3,1)) < hlim)
			kflag = fix(1+kflag);
		end
		count = count + 1;
	end

	[alat,along,rho] = trans5(qhat);
	ahat = trans4(qhat);

	% replace sigma-tilde with (ahat**t)*sigma-tilde*ahat
	for (isect = 1:nsect)
		temp = zeros(3,3);
		for (i = 1:3)
			for (j = 1:3)
				for (k1 = 1:3)
					for (k2 = 1:3)
						temp(i,j) = temp(i,j) + ahat(k1,i) * sigma(2,isect,k1,k2) * ahat(k2,j);
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
	sigma3 = zeros(ndim, ndim);
	k = 0;
	for (i = 1:ndim)
		for (j = 1:i)
			k = k+1;
			sigma3(i,j) = sigma2(k);
			sigma3(j,i) = sigma2(k);
		end
	end

	% sigma3 is the full variance-covariance matrix of the rotation and
	% normal section parameters (except possibly for division by hatkap).
	% rotation parameters are stored in upper 3 by 3 corner.

	% [D,V,nrot]=jacobi(sigma3,ndim,size(sigma3,1));
	[V,D] = eig(sigma3);
	[D, ind] = sort(diag(D));		% Make sure that eigs are in ascending order
	V = V(:,ind);
	sigma3 = zeros(ndim);
	for (i = 1:ndim)
		for (j = 1:ndim)
			for (k = 1:ndim)
				sigma3(i,j) = sigma3(i,j) + V(i,k)*V(j,k) / (D(k) * rfact);
			end
		end
	end

	% calculation of critical point
	df = ndata - 2*nsect - 3;
	plev1 = (1-plev) / 2;
	[xchi1,ier] = xidch(plev1,df);
	if (ier ~= 0),	fprintf('subroutine xidch--ier: %g\n', ier),	end
	plev1 = 1 - plev1;
	[xchi,ier] = xidch(plev1,df);
	if (ier ~= 0),	fprintf('subroutine xidch--ier: %g\n', ier),	end
	fkap1 = xchi / rmin;
	fkap2 = xchi1 / rmin;
	hatkap = df / rmin;
	[xchi,ier] = xidf(plev, 3, df);
	if(ier ~= 0),	fprintf('subroutine xidf--ier: %g\n', ier),	end
	xchi = xchi * rmin * 3 / df;

	% calculation of h11.2 matrix--stored sigma4
	cov = sigma3(1:3,1:3);
	if (any(isnan(cov)))
		errordlg('Boom: Covariance is NaNs.', 'Error')
		vol = NaN;	ellip_long = NaN;	ellip_lat = NaN;	t_stats = {'Boom'};
		return
	end
	[V,D]    = eig(cov);
	[D, ind] = sort(diag(D));		% Make sure that eigs are in ascending order
	V = V(:,ind);
	sigma4 = zeros(3);
	for (i = 1:3)
		for (j = 1:3)
			for (k = 1:3)
				sigma4(i,j) = sigma4(i,j) + V(i,k)*V(j,k) / D(k);
			end
		end
	end

	qbing = zeros(1,10);	ahatm = zeros(3,4);
	ahatm(1,1)=-qhat(2);	ahatm(1,2)=qhat(1);		ahatm(1,3)=qhat(4);		ahatm(1,4)=-qhat(3);
	ahatm(2,1)=-qhat(3);	ahatm(2,2)=-qhat(4);	ahatm(2,3)=qhat(1);		ahatm(2,4)=qhat(2);
	ahatm(3,1)=-qhat(4);	ahatm(3,2)=qhat(3);		ahatm(3,3)=-qhat(2);	ahatm(3,4)=qhat(1);
	k = 0;
	for (i=1:4)
		for (j=1:i)
			k=k+1;
			for (k1=1:3)
				for (k2=1:3)
					qbing(k) = qbing(k) + 4 * ahatm(k1,i) * sigma4(k1,k2) * ahatm(k2,j);
				end
			end
		end
	end
	
	[x,y,z] = trans1(alat, along);		% geo2cart
	axis = [x y z];
 	[ellip_long, ellip_lat,jer] = bingham(qbing,0,xchi,axis);

	% compute volume of confidence  region: vol=4/3*pi*a*b*c following the conreg.f program from J-Y Royer
	[V,D]  = eig(sigma4);
	D      = diag(D);
	angled = sqrt(xchi ./ D) / D2R;
	vol = 4/3 * pi * prod(angled) * 111.111^3;

	t_stats = cell(13,1);
	t_stats{1}  = sprintf('Results from Hellinger1 using in memory segmentation');
	t_stats{2}  = 'Junk';
	t_stats{3}  = sprintf('Fitted rotation--alat,along,rho:');
	t_stats{4}  = sprintf('%f\t%f\t%f', alat,along,rho);
	t_stats{5}  = sprintf('conf. level, conf. interval for kappa:');
	t_stats{5}  = sprintf('%f\t%f\t%f', plev,fkap2,fkap1);
	t_stats{6}  = sprintf('kappahat, degrees of freedom,xchi:');
	t_stats{7}  = sprintf('%f\t%d\t%f', hatkap,df,xchi);
	t_stats{8}  = sprintf('Number of points, sections, misfit, reduced misfit:');
	t_stats{9}  = sprintf('%d\t%d\t%f\t%f', ndata, nsect, rmin, 1/sqrt(hatkap));
	t_stats{10}  = sprintf('ahat:\n%f\t%f\t%f\n%f\t%f\t%f\n%f\t%f\t%f', ahat');
	t_stats{11} = sprintf('covariance matrix:\n%g\t%g\t%g\n%g\t%g\t%g\n%g\t%g\t%g', cov');
	t_stats{12} = sprintf('H11.2 matrix:\n%f\t%f\t%f\n%f\t%f\t%f\n%f\t%f\t%f', sigma4');

% ---------------------------------------------------------------------------
function [r, eta, etai] = r1(h, sigma, qhati, nsect, eta, etai)

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
		eta(i,1,1) = 0;
		eta(i,2,1) = V(3,1);
		eta(i,3,1) = -V(2,1);
		eta(i,1,2) = -V(3,1);
		eta(i,2,2) = 0;
		eta(i,3,2) = V(1,1);
		eta(i,1,3) = V(2,1);
		eta(i,2,3) = -V(1,1);
		eta(i,3,3) = 0;
		if (abs(V(3,1)) > 0.2)
			etai(i,1:3,1) = eta(i,1:3,1);
			etai(i,1:3,2) = eta(i,1:3,2);
		elseif (abs(V(1,1)) > 0.2)
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
	amhat(2,1)=2*(A0*A3+A1*A2);
	amhat(3,1)=2*(A1*A3-A0*A2);
	amhat(1,2)=2*(A1*A2-A0*A3);
	amhat(2,2)=A0*A0-A1*A1+A2*A2-A3*A3;
	amhat(3,2)=2*(A0*A1+A2*A3);
	amhat(1,3)=2*(A0*A2+A1*A3);
	amhat(2,3)=2*(A2*A3-A0*A1);
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
function [P, Y, eta, etai, FTOL] = amoeba(P,Y,MP,NDIM,FTOL,ITER, funk, varargin)
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
		[Y(i), eta, etai] = feval(funk,P(:,i), varargin{:});
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
	
		[YPR, eta, etai] = feval(funk,WORK(:,2), varargin{:});
		if (YPR <= Y(ILO))
			WORK(1:NDIM,3) = GAMMA*WORK(1:NDIM,2) + (1-GAMMA)*WORK(1:NDIM,1);
			[YPRR, eta, etai] = feval(funk, WORK(:,3), varargin{:});
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
			WORK(1:NDIM,3) = BETA*P(1:NDIM,IHI)+(1-BETA)*WORK(1:NDIM,1);
          	[YPRR, eta, etai] = feval(funk,WORK(:,3), varargin{:});
			if (YPRR < Y(IHI))
				P(1:NDIM,IHI)=WORK(1:NDIM,3);
				Y(IHI)=YPRR;
			else
				for (i = 1:MPTS)
					if (i ~= ILO)
						WORK(1:NDIM,2)=0.5*(P(1:NDIM,i)+P(1:NDIM,ILO));
						P(1:NDIM,i)=WORK(1:NDIM,2);
						[Y(i), eta, etai] = feval(funk,WORK(:,2), varargin{:});
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
function [ZBR, NER] = zbrent(FUNC,FVAL,DUMMY,X1,X2,TOL)
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
			if (2 * P < min(3.*XM*Q-abs(TOL1*Q),abs(E*Q)))
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
			XN = zbrent('xdn',QLEV,0,TABLE(I-1),TABLE(I),TOL);
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
		XCHI=0;
		IER=111111;
		return
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
					XCHI = zbrent('xdch',PLEV,D,A,B,TOL);
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
	XCHI = zbrent('xdch',PLEV,D,A,B,TOL);

% -------------------------------------------------------------------------------------
function [XF,IER] = xidf(PLEV,D1,D2)
	TOL = 0.0001;
	D = [D1 D2];
	if ((PLEV <= 0) || (PLEV >= 1) || (D1 <= 0) || (D2 <= 0))
		XF = 1;
		IER=111111;
		return
	end
	[XF,IER] = xidch(PLEV,D(1));
	A = XF;
	PA = xdf(D,A,IER);
	while (PA >= PLEV)
		A = A / 2;
		PA = xdf(D,A,IER);
	end
	B = XF;
	PB = xdf(D,B,IER);
	while (PB <= PLEV)
		B = 2 * B;
		PB = xdf(D,B,IER);
	end
	XF = zbrent('xdf',PLEV,D,A,B,TOL);

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
function [g, IER] = gammp(A,X,IER)
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
		GAMMCF = gcf__(A,X,GLN,IER);
		g = 1 - GAMMCF;
	end

% -------------------------------------------------------------------------------------
function [GAMSER, IER] = gser(A,X,GLN,IER)
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
function GAMMCF = gcf__(A,X,GLN,IER)
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
function [B, IER] = betai(A,B,X,IER)
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
function [B, IER] = betacf(A,B,X,IER)
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

% -------------------------------------------------------------------------------------
function mixed = classical_hellinger_order(flags_w, isow, flags_e, isoe)
% Join the two isochrons ISOW & ISOE and their flags indicating to which Hellinger
% segment belongs each point in the isoc into one single array where points are
% ranged by growing order os sgment number. That is, the classical hellinger file organization

	max_seg = max(flags_w(end), flags_e(end));		% Get the greater  segment number
	
	mixed = zeros(size(flags_w,1) + size(flags_e,1), 2);
	off = 0;
	for (k = 1:max_seg)
		ind = find(flags_w == k);
		mixed(off+(1:numel(ind))', :) = isow(ind, :);
		if (~isempty(ind)),		off = off + numel(ind);	end
		ind = find(flags_e == k);
		mixed(off+(1:numel(ind))', :) = isoe(ind, :);
		if (~isempty(ind)),		off = off + numel(ind);	end
	end

% -----------------------------------------------------------------
function flags = segmentate(isoc, segmented)
% Take the Mx2 (lon lat) ISOC array and compare it to the segmemted multi-segment
% line SEGMENTED. This line has actually all the segments delimited by NaNs, so we must
% first extract them into a cell array.
% Next mapproject does the magic of finding to which segment belongs each point of ISOC

	% First, extracts the segments from the NaN-delimited array
	% Make sure last row is a NaNs one
	if (~isnan(segmented(end,1))),	segmented = [segmented; NaN NaN];		end
	ind = find(isnan(segmented(:,1)));
	cells = cell(numel(ind), 1);		% Assumes first and last rows are NOT NaNs
	cells{1} = segmented(1:ind(1)-1,:);
	for (i = 2:numel(ind))
		cells{i} = segmented(ind(i-1)+1:ind(i)-1, :);
	end

	out = gmtmex('mapproject -L+uk+p', isoc, gmt('wrapseg', cells));

	% Now the flags are just the fourth column. Just beautifull
	flags = out.data(:,4) + 1;		% + 1 because mapproject is C and so 0 based

% ------------------------------------------------------------------------------
function set_extra_uicb_options(hLine, data)
% Reuse two entries (in this context) of the polygon's UIContextMenu to allow
% saving a .pick file or showing patches with the conjugated segments.

	h1 = get(get(hLine,'UIContextMenu'),'Children');
	h2 = findobj(h1,'-depth',0, 'Label','Join lines');	% Resuse this entry that has no sense in this context
	h3 = findobj(h1,'-depth',0, 'Label','Copy');
	set(h2,'Label','Save as .pick file')
	set(h3,'Label','Show conjugated picks')
	hMirFig = get(get(hLine, 'Parent'), 'Parent');
	set(h2,'Call', {@save_pick, data, hMirFig})
	set(h3,'Call', {@show_conjugates, hLine, data})
	h4 = findobj(h1,'-depth',0, 'Label','Spline Smooth');
	h5 = findobj(h1,'-depth',0, 'Label','Line length(s)');
	h6 = findobj(h1,'-depth',0, 'Label','Line azimuth(s)');
	h7 = findobj(h1,'-depth',0, 'Label','Extract profile');
	delete([h4 h5 h6 h7])		% Need to delay deletion because h1 has them all and can't have invalid handles

function save_pick(obj, evt, data, hMirFig)
% Save the contents of the DATA array in a .pick file (Hellinger stuff)
	handles = guidata(hMirFig);
	txt1 = 'Hellinger pick file (*.pick)';	txt2 = 'Select output Hellinger pick file format';
	[FileName,PathName] = put_or_get_file(handles,{'*.pick',txt1; '*.*', 'All Files (*.*)'},txt2,'put','.pick');
	if isequal(FileName,0),		return,		end
	f_name = [PathName FileName];
	fid = fopen(f_name, 'wt');
	fprintf(fid, '%d  %d  %.4f  %.4f  %.1f\n', data');
	fclose(fid);

function show_conjugates(obj, evt, hLine, data)
% Show colored patches illustrating the result of the automatic Hellinger segmentation.
	handles = guidata(hLine);
	cmenuHand = uicontextmenu('Parent',handles.figure1);
	left  = data(data(:,1) == 1, 2:4);
	right = data(data(:,1) == 2, 2:4);
	n_segs = min(left(end,1), right(end,1));
	for (k = 1:n_segs)
		ind_l = find(left(:,1)  == k);
		ind_r = find(right(:,1) == k);
		if (isempty(ind_l) || isempty(ind_r)),	continue,	end
		x = [left(ind_l,3); right(ind_r(end:-1:1),3)];		% data is in lat long
		y = [left(ind_l,2); right(ind_r(end:-1:1),2)];
		hP = patch('Parent',handles.axes1,'XData',x, 'YData',y, 'FaceColor',rand(1,3), 'Tag', 'HellPair');
		set(hP, 'UIContextMenu', cmenuHand);
	end
	uimenu(cmenuHand, 'Label', 'Delete all', 'Call', 'delete(findobj(''Tag'',''HellPair''))')

% --------------------------------------------------------------------------------------------
function update_patches(hLine, data, show_segs)
% Create patches showing the segmentation or update existent ones.
% In this later case just delete old ones and recreate.
	hP = findobj('Tag','HellPair');
	if (~isempty(hP))
		delete(hP)
		if (~show_segs),	return,		end		% Ok, do not reconstruct them then.
		show_conjugates([], [], hLine, data)
	elseif (show_segs)		% Here we don't have them yet but they were required.
		show_conjugates([], [], hLine, data)
	end

% last modified september 16, 1988

% subroutine to generate points on the confidence region boundary, for
% use by surface2, or other contouring program.
% input:
%   q    - the matrix q (in symmetric storage mode)
%   lhat - the smallest eigenvalue of q
%   fk   - constant used in determining the confidence region
%   uxx  - optimal axis of rotation
% output:
%   jer  - error indicator
%     0  - all is well
%     1  - an error has occurred
% output files:
%   30 - bounding curve
%   31 - upper surface
%   32 - lower surface

% The original fortran code had no License so I'm releasing this partially finished translation 
% that is part of Mirone, to the Public Domain
%
%	Copyright (c) 2004-2017 by J. Luis
%
%	Contact info: w3.ualg.pt/~jluis/mirone
% --------------------------------------------------------------------

% $Id$

function [blong,blat,jer] = bingham(q,lhat,fk,uxx)

	pmind=[];pmaxd=[];tmind=[];tmaxd=[];
	%fprintf('\nDetermining the confidence region in the form of min. and max. angles\n');
	%fprintf('of rotation expressed as functions of longitude and latitude axes.\n');

	%ind=0;
	ind=1;		% Do graphics	
	qf = zeros(4,4);
	nu = zeros(1,3);
	w  = zeros(3,3);
	mf = zeros(3,3);
	jer=0;
	% modify the matrix q by subtracting lhat times the identity from q, and dividing the result by fk.
	qt = q;
	qt(1) = qt(1) - lhat;
	qt(3) = qt(3) - lhat;
	qt(6) = qt(6) - lhat;
	qt(10)= qt(10) - lhat;
	qt = qt / fk;
	k = 0;
	for i=1:4
		for j=1:i
			k=k+1;
			qf(i,j)=qt(k);
			qf(j,i)=qf(i,j);
		end
	end

	%fprintf('The matrix Qtilde:\n')
	%for (i=1:4)
		%fprintf('%f\t%f\t%f\t%f\n', qf(i, :))
	%end
	% test whether the identity is in the confidence region.  if not, calculate the matrix
	%     m = (azero - 1)*b - dvect*transpose(dvect)
	% and find its eigenvalues and eigenvectors.
	azero=qf(1,1);
	if (azero <= 1.0)
		icase=7;
	else
		[nu,w,mf,icase,ier] = meig(qf,uxx,nu,w,mf);
		if (ier ~= 0)
			jer=1;
			return
		end
	end
	% print a description of the set a of admissible axes, depending on the value of icase.
	if (icase == 1)		
		fprintf('\nThe set a of admissible axes is a cap not containing either pole.\n');
		fprintf('In the longitude-latitude plane, this set is bounded by a closed curve.\n');
	elseif (icase == 2)
		fprintf('\nthe set a of admissible axes is a cap containing\n')
		fprintf('the north pole.  in the axis longitude-axis latitude\n')
		fprintf('plane this set is bounded by the lines: axis longi\n')
		fprintf('tude = -180 degrees, axis longitude = 180 degrees\n')
		fprintf('axis latitude = 90 degrees, and a curve which forms\n')
		fprintf('the southern border.\n')
	elseif (icase == 3)
		fprintf('\nthe set a of admissible axes is a cap containing');
		fprintf('the south pole.  in the axis longitude-axis latitude');
		fprintf('plane this set is bounded by the lines: axis longi-');
		fprintf('tude = -180 degrees, axis longitude = 180 degrees,');
		fprintf('axis latitude = -90 degrees, and a curve which forms');
		fprintf('the northern border.');
	elseif (icase == 4)
		fprintf('\nthe set a of admissible axes is the complement of');
		fprintf('two anti-podal caps which contain the poles.');
		fprintf('hence, this set is an equatorial belt.  in the axis');
		fprintf('longitude-axis latitude plane this set is bounded by');
		fprintf('the lines: axis longitude = -180 degrees, axis longi-');
		fprintf('tude = 180 degrees, and two curves which form the');
		fprintf('northern and southern borders of the belt.');
	elseif (icase == 5)
		fprintf('\nthe set a of admissible axes is the complement of');
		fprintf('two anti-podal caps which do not contain the');
		fprintf('poles.  in the axis longitude-axis latitude plane');
		fprintf('this set is bounded by the four lines: axis longi-');
		fprintf('tude = -180 degrees, axis longitude = 180 degrees,');
		fprintf('axis latitude = -90 degrees, and axis latitude =');
		fprintf('90 degrees; however, there are two holes in the set.');
	elseif (icase == 6)
		fprintf('\nany axis is admissible; however, the identity is not');
		fprintf('in the confidence region.  the set a of admissible');
		fprintf('axes is the entire axis longitude-axis latitude plane.');
	elseif (icase == 7)
		fprintf('\nthe identity is in the confidence region.  hence, each');
		fprintf('axis is admissible.  that is, the set a of admissible');
		fprintf('axes is the entire axis longitude-axis latitude plane.');
	end

	[w,icase,ind,pmind,pmaxd,tmind,tmaxd,blong,blat,ier] = boundc(azero,nu,w,icase,ind,pmind,pmaxd,tmind,tmaxd);
	if (ier ~= 0)
		jer=1;
		return
	end
	alat1=(largin(tmind));
	alat2=(-largin(-tmaxd));
	along1=(largin(pmind));
	along2=(-largin(-pmaxd));
% 	[alat1,alat2,along1,along2,qf,mf,w,drho,crho,icase,ind,ier] = grid(alat1,alat2,along1,along2,qf,mf,w,icase,ind);
% 	if (ier ~= 0)
% 		jer=1;
% 	end

% ------------------------------   boundc   -----------------------------------------------
% subroutine to calculate points on the bounding curve.
% inputs:
%   azero   - qf(1,1)
%   nu      - array containing the eigenvalues of m
%   w       - array containing the eigenvectors of m
%   icase   - indicator of the type of set a
%   ind     - indicator
%     ind=0:  box out the confidence region
%     ind=1:  write the boundary contour to file 30
% outputs:
%   pmind, pmaxd - min. and max. longitudes over the confidence region
%   tmind, tmaxd - min. and max. latitudes over the confidence region
%   jer          - error indicator
%       0        - all is well
%       1        - an error has occurred
% defined by parameter statements:
%   npart   - no. of parts into which to divide the bounding curve,
%             to get the step-size
%   kmax    - max. no. of points (must be at least 13)
% the subroutine boundc writes the longitude and latitude values for the
% points on the bounding curve to file 30, provided ind=1 and
% icase<=5.

function [w,icase,ind,pmind,pmaxd,tmind,tmaxd,blong,blat,jer] = ...
          boundc(azero,nu,w,icase,ind,pmind,pmaxd,tmind,tmaxd, jer,varargin)

format_2002=[ '\n ' ,'   Type name of boundary file (max.=10 characters.)', '\n ' ,'%7x','i.e., file to receive points on the boundary of the set',' of admissible axes.'];
format_703=['   A sequence of ','%3d',' points around the closed',' curve bounding the set of admissible', '\n ' ,'   axes has',' been written to the file ','%10c','.'];
format_712=[' ','  a sequence of ','%3d',' points along the southern',' border', '\n ' ,' ','  of the set of admissible axes has been',' written to', '\n ' ,' ','  the file ','%10c','.'];
format_713=[' ','  a sequence of ','%3d',' points along the northern',' border', '\n ' ,' ','  of the set of admissible axes has been',' written to', '\n ' ,' ','  the file ','%10c','.'];
format_714=[' ','  a sequence of ','%3d',' points along the curve',' which forms the', '\n ' ,' ','  northern border of the belt have been',' written to the file ','%10c','.', '\n ' ,' ','  these are followed by the',' same number of points along the', '\n ' ,' ','  curve which forms the',' southern border.'];
format_715=[' ','  a sequence of ','%3d',' points around each of',' the closed', '\n ' ,' ','  curves has been written to the file ','%10c','.'];
format_716=[' ','  a sequence of ','%3d',' points around the closed',' curve', '\n ' ,' ','  have been written to the file ','%10c','; these are',' followed', '\n ' ,' ','  by ','%3d',' points one one open curve and ','%3d',' points on the', '\n ' ,' ','  other.'];
format_801=['%4d','   1'];
format_802=['%1x',repmat('%16.8e',1,2)];
format_803=['%4d','   0'];
format_805=['%4d','   2'];
format_815=[' ','  neither of the holes straddles the line:',' axis', '\n ' ,' ','  longitude = 180 degrees.  hence, each hole is',' bounded by a', '\n ' ,' ','  closed curve.'];
format_816=[' ','  one of the holes straddles the line: axis',' longitude =', '\n ' ,' ','  180 degrees.  hence, this hole is bounded',' by two open curves', '\n ' ,' ','  (one near axis longitude = 180',' degrees, the other near', '\n ' ,' ','  axis longitude = -180',' degrees).  the other hole is bounded', '\n ' ,' ','  by a closed',' curve.'];

	D2R = pi / 180;		R2D = 180 / pi;		TWOPI = 2 * pi;
	jer=0;
	npart=300;
	kmax=400;
	pmind = 0;
	pmaxd = 0;
	tmind = 0;
	tmaxd = 0;

	phi=zeros(1,kmax);
	theta=zeros(1,kmax);
	phit=zeros(1,kmax);
	phim=zeros(1,kmax);
	thetm=zeros(1,kmax);
	oldu=zeros(1,3);

	if (icase >= 6)
		pmind= -180.0;
		pmaxd=  180.0;
		tmind=  -90.0;
		tmaxd=   90.0;
		fprintf('Min., Max. longitude: %f10.5\t%f10.5\nMin., Max. latitude: %f10.5\t%f10.5\n',pmind,pmaxd,tmind,tmaxd)
		return
	end

	if (ind == 1)
		% Printagens em fiche
% 		fprintf(format_2002);
% 		bfile=strAssign(bfile,[],[],input('','s'));
% 		thismlfid=fopen(strtrim(bfile),'r+');
% 		unit2fid=[unit2fid; 30,thismlfid];
% 		frewind(unit2fid(find(unit2fid(:,1)== 30,1,'last'),2));
	end

	% calculation of estimated length of the bounding curve.  this is
	% done by calculating thirteen points around the curve, at equal increments in phit.
	delpt = pi/6;
	ophi = 0;
	for (k=1:13)
		phit(k)=(k-1).*delpt;
		[fval,phi(k),theta(k),u,ier] = evalf(phit(k), nu, w, azero, k-1, oldu, ophi);		
		if (ier ~= 0)
			jer=1;
			return
		end
		oldu = u;
		ophi=phi(k);
	end

	eleng=0;
	for  k=1:12;
		dist=sqrt((phi(k+1) - phi(k)).*(phi(k+1) - phi(k))+(theta(k+1) - theta(k)).*(theta(k+1) - theta(k)));
		eleng=eleng + dist;
	end
	elengd = eleng * R2D;
	dels   = eleng./(npart);
	phit(1)=0;
	if (icase == 3)
		dels= -dels;
	end
	% main loop.  calculate the desired phit values by applying euler's
	% method to the differential equation:
	%     dphit/ds = f(phit) ,
	% where f(phit) = 1./sqrt((dphi/dphit)**2 + (dtheta/dphit)**2).
	% here s is the arc length along the bounding curve, measured in radians.
	did_break = false;
	for  k=2:kmax
		[fval,phi(k-1),theta(k-1),u,ier] = evalf(phit(k-1), nu, w, azero, k-2, oldu, ophi);
		if (ier ~= 0)
			jer=1;
			return
		end
		oldu = u;
		ophi = phi(k-1);
		phit(k) = phit(k-1) + fval * dels;
		
		if ((abs(phit(k)) >= TWOPI))
			did_break = true;			% Best GOTO simulation
			break
		end
	end

	if (~did_break)
		jer = 1;
		fprintf('error in boundc:\n');
		fprintf('unable to close the bounding curve with the prescribed\n');
		fprintf('npart and kmax.  either npart needs to be decreased or\n');
		fprintf('kmax needs to be increased.\n');
		return
	end

	if (icase == 1)		% The set a is a cap not containing either pole. Center the phi values, if necessary.
		npt=k;
		blong = zeros(npt, 1);
		blat = zeros(npt, 1);
		phit(npt) = TWOPI;
		phi(npt)  = phi(1);
		theta(npt)=theta(1);
		pmin=phi(1);
		pmax=phi(1);
		for  (k_=1:npt)
			if (phi(k_) < pmin)
				pmin=phi(k_);
			elseif (phi(k_) > pmax)
				pmax=phi(k_);
			end
		end
		pmean = (pmin + pmax)/2;
		if (pmean > pi)
			phi(1:npt) = phi(1:npt) -TWOPI;
			pmin=pmin - TWOPI;
			pmax=pmax - TWOPI;
		elseif (pmean <= -pi)
			phi(1:npt) = phi(1:npt) + TWOPI;
			pmin=pmin + TWOPI;
			pmax=pmax + TWOPI;
		end
		tmin=theta(1);
		tmax=theta(1);
		for (k_=1:npt)
			if (theta(k_) < tmin)
				tmin=theta(k_);
			elseif (theta(k_) > tmax)
				tmax=theta(k_);
			end
		end
		pmind = pmin * R2D;
		pmaxd = pmax * R2D;
		tmind = tmin * R2D;
		tmaxd = tmax * R2D;
		%fprintf('Min., Max. longitude: %10.5f\t%10.5f\nMin., Max. latitude: %10.5f\t%10.5f\n',pmind,pmaxd,tmind,tmaxd)
		if (ind == 1)
			% Printagens
% 				writeFmt(30,format_801,'npt');
				for  k_=1:npt
					blong(k_) = phi(k_) * R2D;
					blat(k_)  = theta(k_) * R2D;
% 					writeFmt(30,format_802,'blong','blat');
				end
% 				try
% 					fclose(unit2fid(find(unit2fid(:,1)==30,1,'last'),2));
% 					unit2fid=unit2fid(unit2fid(:,1)~=30,:);
% 				end
% 				fprintf(format_703,npt,bfile);
		end
		return
	elseif ((icase >= 2) && (icase <= 4))
		%   icase = 2:  the set a is a cap containing the north pole.
		%   icase = 3:  the set a is a cap containing the south pole.
		%   icase = 4:  the set a is the complement of two anti-podal caps
		%               which contain the poles.
		% calculate a sequence of points (phim(k),thetm(k),k=1,2,...,npt,
		% which are arranged in order on phi, with phi going from -180 to 180.

		npt=k - 1;
		for (k_=2:npt)
			if (phi(k_) <= pi)
				for  (k__ = 1:npt)
					phim(k__)=phi(k__);
					thetm(k__)=theta(k__);
				end
				go to 225;
			else
				% 215 continue;
				kz=k-1;
				for (l=1:npt-kz)
					phim(l)=phi(kz + l) - TWOPI;
					thetm(l)=theta(kz + l);
				end
				for (l=1:kz)
					k2 = l + npt - kz;
					phim(k2)=phi(l);
					thetm(k2)=theta(l);
				end
				% from now on we treat the cases icase=2,3 and icase=4 differently.					
			end
		end

		% 225 continue;
		if (icase ~= 4)
			tmin=thetm(1);
			tmax=thetm(1);
			for  k_=1:npt
				if (thetm(k_) < tmin)
					tmin=thetm(k_);
				elseif (thetm(k_) > tmax)
					tmax=thetm(k_);
				end
			end
			pmind= -180;
			pmaxd=  180;
			if (icase == 2)
				tmind=tmin * 180/pi;
				tmaxd=90;
			elseif (icase == 3)
				tmind= -90;
				tmaxd=tmax * 180/pi;
			end

			fprintf('Min., Max. longitude: %10.5f\t%10.5f\nMin., Max. latitude: %10.5f\t%10.5f\n',pmind,pmaxd,tmind,tmaxd)
			if (ind == 1)
				% Printagens
% 					writeFmt(30,format_803,'npt');
% 					for  k=1:npt
% 						blong=phim(k).*180/pi;
% 						blat=thetm(k).*180/pi;
% 						writeFmt(30,format_802,'blong','blat');
% 					end
% 					try
% 						fclose(unit2fid(find(unit2fid(:,1)==30,1,'last'),2));
% 						unit2fid=unit2fid(unit2fid(:,1)~=30,:);
% 					end
% 					if (icase == 2)
% 						[writeErrFlag]=writeFmt(1,[format_712],'npt','bfile');
% 					elseif (icase == 3) ;
% 						[writeErrFlag]=writeFmt(1,[format_713],'npt','bfile');
% 					end
			end
			return
		end

		% icase = 4.  the set a is the complement of two anti-podal caps
		% which contain the poles.  kzz is the number of points (phim(k), thetm(k)) with phi negative.

		% 400 continue;
		for (k_=1:npt)
			if (phim(k_) <= 0)
				jer=1;
				fprinf('error in boundc (icase=4):\n');
				fprinf('there are no points on the bounding curve with positve\n');
				fprinf('longitude.  the points on the bounding curve are too sparse\n');
				return
			end
			kzz = k-1;
			if (kzz <= 0)
				jer=1;
				fprinf('error in boundc (icase=4):\n');
				fprinf('there are no points on the bounding curve with\n');
				fprinf('negative longitude.  the points on the bounding\n');
				fprinf('curve are too sparse.\n');
				return
			end
		end

		tmax=thetm(1);
		for  k_=1:npt
			if (thetm(k_) > tmax)
				tmax=thetm(k_);
			end
		end
		pmind= -180.;
		pmaxd=  180.;
		tmaxd=tmax * 180/pi;
		tmind= -tmaxd;
		fprintf('Min., Max. longitude: %10.5f\t%10.5f\nMin., Max. latitude: %10.5f\t%10.5f\n',pmind,pmaxd,tmind,tmaxd)
		if (ind == 1)
			% Printagens
% 				writeFmt(30,format_803,'npt');
% 				for  k=1:npt;
% 					blong=phim(k) * 180/pi;
% 					blat=thetm(k) * 180/pi;
% 					writeFmt(30,format_802,'blong','blat');
% 				end
% 				writeFmt(30,format_803,'npt');
% 				for  l=kzz+1:npt;
% 					blong=phim(l).*180/pi - 180.;
% 					blat= -thetm(l).*180/pi;
% 					writeFmt(30,format_802,'blong','blat');
% 				end
% 				for  l=1:kzz
% 					blong=phim(l).*180/pi + 180.;
% 					blat= -thetm(l).*180/pi;
% 					writeFmt(30,format_802,'blong','blat');
% 				end
% 				try
% 					fclose(unit2fid(find(unit2fid(:,1)==30,1,'last'),2));
% 					unit2fid=unit2fid(find(unit2fid(:,1)~=30),:);
% 				end
% 				writeFmt(1,format_714,'npt','bfile');
		end
		return

	elseif (icase == 5)
		% The set a is the complement of two anti-podal caps
		% which do not contain the poles.  if the cap which we have crosses
		% either of the lines phi=pi or phi= -pi, we replace it by the
		% anti-podal cap.  we then have subcase a or subcase b, depending
		% on whether the cap which we have obtained does not cross the
		% line phi=0, or does cross it.

		npt=k;
		phi(npt)   = phi(1);
		theta(npt) = theta(1);
		pmind = -180;
		pmaxd =  180;
		tmind = -90;
		tmaxd =  90;
		fprintf('Min., Max. longitude: %10.5f\t%10.5f\nMin., Max. latitude: %10.5f\t%10.5f\n',pmind,pmaxd,tmind,tmaxd)
		for (k_=1:npt)
			if (phi(k_) > pi)
				for  (k__ =1:npt)
					phi(k__)   = phi(k__) - pi;
					theta(k__) = - theta(k__);
				end
			else
				for (k__=1:npt)
					if (phi(k__) <= -pi)
						for (k3 = 1:npt)
							phi(k3)   = phi(k3) + pi;
							theta(k3) = -theta(k3);
						end
					end
				end
				ipos=0;
				ineg=0;
				for  (k__=1:npt)
					if (phi(k__) > 0),	ipos=1;
					else				ineg=1;
					end
				end
			end

			if ~((ipos == 1) && (ineg == 1))
				% subcase a:  the cap which we have does not cross the line phi=0.
				% therefore, the bounding curve consists of two closed curves: the
				% one which we have and its anti-pode.  isign is a switch used in
				% determining the anti-podal curve.
				if ((ipos == 1) && (ineg == 0)),	isign= -1;	end
				if ((ipos == 0) && (ineg == 1)),	isign=  1;	end
				if ((ipos ~= 1) && (ineg ~= 1))
					jer=1;
					fprintf('error in boundc (icase=5, subcase a):\n');
					fprintf('the curve we have lies on the line phi=0, which\n');
					fprintf('is very unlikely.\n');
					return
				end

				if (ind == 1)
					% Printagens
% 						writeFmt(30,format_805,'npt');
% 						for  k=1:npt;
% 							blong=phi(k).*180/pi;
% 							blat=theta(k).*180/pi;
% 							writeFmt(30,format_802,'blong','blat');
% 						end
% 						writeFmt(30,format_805,'npt');
% 						for  k=1:npt;
% 							blong=phi(k).*180./pi + isign.*180.;
% 							blat= -theta(k).*180./pi;
% 							writeFmt(30,[format_802],'blong','blat');
% 						end
% 						try
% 							fclose(unit2fid(find(unit2fid(:,1)==30,1,'last'),2));
% 							unit2fid=unit2fid(unit2fid(:,1)~=30,:);
% 						end
% 						writeFmt(1,format_815);
% 						writeFmt(1,format_715,'npt','bfile');
				end
				return
			else
				% subcase b:  the cap which we have crosses the line phi=0.  therefore,
				% the anti-podal cap is broken into two pieces by the line phi= +/-pi.
				% the bounding curve consists of three pieces: one closed curve and two open curves.

				if (phi(1) >= 0),	isw = 1;
				else				isw = -1;
				end

				for  (kk=2:npt)
					if (phi(kk) >= 0),	ksw =  1;
					else				ksw = -1;
					end
					if (ksw * isw >= 0)
						jer=1;
						fprintf('error in boundc (icase=5, subcase b):\n');
						fprintf('impossible exit from loop, statement 575.\n');
						return
					else
						k1 = k-1;
						for (kk_ = k1+1:npt)
							if (phi(kk_) >= 0),	ksw=  1;
							else				ksw= -1;
							end
							if (ksw * isw <= 0)
								jer=1;
								fprintf('error in boundc (icase=5, subcase b):\n');
								fprintf('impossible exit from loop, statement 585.\n');
								return
							else
								k2=kk_-1;
								n1=k2 - k1;
								n2=npt - 1 -n1;
								n3=npt - 1 - k2;
								if (ind == 1)
									% Printagens
% 										writeFmt(30,format_805,'npt');
% 										for  k=1:npt
% 											blong=phi(k).*180/pi;
% 											blat=theta(k).*180/pi;
% 											writeFmt(30,format_802,'blong','blat');
% 										end
% 										writeFmt(30,format_803,'n1');
% 										for  l=1:n1
% 											k=k1 + l;
% 											blong=phi(k).*180/pi + isw.*180.;
% 											blat= -theta(k).*180/pi;
% 											writeFmt(30,format_802,'blong','blat');
% 										end
% 										writeFmt(30,format_803,'n2');
% 										if (n3 > 0)
% 											for  l=1:n3
% 												k=k2 + l;
% 												blong=phi(k).*180/pi - isw.*180.;
% 												blat= -theta(k).*180/pi;
% 												writeFmt(30,format_802,'blong','blat');
% 											end
% 										end
% 										for  l=n3+1:n2
% 											k=l - n3;
% 											blong=phi(k).*180../pi - isw.*180.;
% 											blat= -theta(k).*180../pi;
% 											writeFmt(30,format_802,'blong','blat');
% 										end
% 										try
% 											fclose(unit2fid(find(unit2fid(:,1)==30,1,'last'),2));
% 											unit2fid=unit2fid(find(unit2fid(:,1)~=30),:);
% 										end
% 										fprintf(format_816);
% 										writeFmt(1,format_716,'npt','bfile','n1','n2');
								end	
							end
						end
					end
				end	
			end
		end
	else
		error('Deu merda')
	end		% icase(s) == i

% -------------------------------------   grid   ----------------------------------------------
% subroutine to calculate the min. and max. values of rho (angle of
% rotation) on a rectangular grid in the axis longitude-axis latitude
% plane.
% inputs:
%   alat1   - latitude value for first row
%   alat2   - latitude value for last row
%   along1  - longitude value for first row
%   along2  - longitude value for last row
%     (these values are in degrees.)
%   qf      - the matrix qtilde in full storage mode
%   mf      - the matrix m in full storage mode
%   w       - array containing the eigenvectors of m
%   icase   - indicator for the type of set a
%   ind     - indicator
%     ind=0:  box out the confidence region
%     ind=1:  write upper and lower surfaces to files 31 and 32, resp.
% outputs:
%   jer     - error indicator
%       0   - all is well
%       1   - an error has occurred
% defined by parameter statements:
%   nlat    - no. of rows in the grid
%   nlong   - no. of cols. in the grid
% if ind=1, the subroutine grid writes the rho values for the upper
% surface to file 31; rho values for the lower surface to file 32.
% (an exception occurs when icase=7: in this case only file 31 is
% used; the min. value of rho is zero for each axis.)

function [alat1,alat2,along1,along2,qf,mf,w,drho,crho,icase,ind,jer] = ...
		grid(alat1,alat2,along1,along2,qf,mf,w,icase,ind,jer,varargin)

	nlat=101;
	nlong=201;

	u=zeros(1,3);
	mtu=zeros(1,3);
	af=zeros(2,2);
	bu=zeros(1,3);
	drho=zeros(1,nlong);
	crho=zeros(1,nlong);
	lfile=repmat(' ',1,10);
	ufile=repmat(' ',1,10);

	jer=0;
	if (ind == 1)
		% Isto deve ser pros nomes dos fiches. Mas nos so quwremos os dados, nao as escritas.
% 	fprintf(format_2003);
% 	ufile=strAssign(ufile,[],[],input('','s'));
% 	thismlfid=fopen(strtrim(ufile),'r+');
% 	unit2fid=[unit2fid;31,thismlfid];
% 	frewind(unit2fid(find(unit2fid(:,1)== 31,1,'last'),2));
% 	if (icase < 7)
% 		fprintf(format_2004);
% 		lfile=strAssign(lfile,[],[],input('','s'));
% 		thismlfid=fopen(strtrim(lfile),'r+');
% 		unit2fid=[unit2fid;32,thismlfid];
% 		frewind(unit2fid(find(unit2fid(:,1)== 32,1,'last'),2));
% 	end
	end

	dlong=(along2 - along1)./(nlong - 1.0);
	dlat=(alat2 - alat1)./(nlat - 1.0);

	% irho is a switch for calculation of min. and max. values of rho.
	%   irho= 0:  no values of rho have been calculated
	%   irho= 1:  at least one value of rho has been calculated
	% when the first value of rho is calculated, irho is switched from 0 to 1.
	irho=0;
	% outer loop: on the latitude values.  the arrays crho and drho
	% contain the values of rho for the lower and upper surfaces, resp.,
	% along the current row.  these arrays are initially set to the
	% value -1000000.; then the correct values for grid points which
	% coorespond to admissible axes are filled in.  the values
	% -1000000. which are left signal non-admissible axes.

for  i=1:nlat
	clat=alat1 +(i-1).*dlat;
	for  j=1:nlong
		drho(j)= -1000000;
		crho(j)= -1000000;
	end
	% inner loop: on the longitude values.  if icase<=5, an axis u is
	% skipped if transpose(u)*m*u>(azero - 1).  in addition, if icase
	% <=3, u is skipped if it is not in the desired cap.  for a partic-
	% ular axis u which is accepted, the permissible range of angles of
	% rotation is calculated based on the matrix a(u).
	for (j=1:nlong)
		clong=along1 +(j-1).*dlong;
		cphi=clong.*pi./180.;
		ctheta=clat.*pi./180.;
		cost=cos(ctheta);
		u(1)=cost.*cos(cphi);
		u(2)=cost.*sin(cphi);
		u(3)=sin(ctheta);
		if (icase <= 5)
			azero=qf(1,1);
			mlim=azero - 1.0;
			for  k=1:3
				mtu(k)=0.0;
				for  l=1:3;
					mtu(k)=mtu(k) + mf(k,l).*u(l);
				end
			end
			umu=0.0;
			for  k=1:3
				umu=umu + mtu(k).*u(k);
			end
			if (umu > mlim)
				continue
			end
			if (icase <= 3)
				dotp=0;
				for  k=1:3
					dotp=dotp + u(k).*w(k,3);
				end
				if (dotp <= 0)
					continue
				end
			end
		end
		af(1,1)=qf(1,1);
		sum=0;
		for  k=1:3
			sum=sum + qf(k+1,1).*u(k);
		end
		af(1,2)=sum;
		af(2,1)=af(1,2);
		for  k=1:3
			bu(k)=0.0;
			for  l=1:3;
				bu(k)=bu(k) + qf(k+1,l+1).*u(l);
			end
		end
		sum=0;
		for  k=1:3
			sum=sum + bu(k).*u(k);
		end
		af(2,2)=sum;
		%[af,dumvar2,dumvar3,mu,v,dumvar6,work,nrot] = jacob2(af,2,2,mu,v,2,work,nrot);
		% if icase<7, arrange that v(2,1)>=0., so that rhostar comes out in
		% the range from 0 to 360 degrees.
		[v,mu] = eig(af);
		[mu, ind] = sort(diag(mu));
		v = v(:,ind);
		if (icase < 7)
			if (v(2,1) < 0)
				for  k=1:2
					v(k,1)= -v(k,1);
				end
			end
		end
		% if icase==7, arrange that v(1,1)>=0., so that rhostar is in the
		% range from -180 to 180 degrees.
		if (icase == 7)
			if (v(1,1) < 0)
				for  k=1:2
					v(k,1)= -v(k,1);
				end
			end
		end
		angle=atan2(v(2,1),v(1,1));
		rhostar=(2*angle) * 180/pi;
		if (mu(1) >= 1);
			rhoinc=0;
		elseif (mu(2) <= 1)
			rhoinc=180;
		else
			td=sqrt((1.0 - mu(1))./(mu(2) - 1.0));
			delta=atan(td);
			rhoinc=(2*delta).*180/pi;
		end
		drho(j)=rhostar + rhoinc;
		crho(j)=rhostar - rhoinc;
		if (irho ~= 1)
			irho=1;
			rmind=crho(j);
			rmaxd=drho(j);
		else
			if (crho(j) < rmind);
				rmind=crho(j);
			end
			if (drho(j) > rmaxd)
				rmaxd=drho(j);
			end
		end
	end				% end of inner loop.

	if (ind == 1)
		% Aqui e que se dava a escrita
% 		writeFmt(31,format_905,{'drho(j)','j','1','1','nlong'});
% 		if (icase < 7);
% 			writeFmt(32,format_905,{'crho(j)','j','1','1','nlong'});
% 		end
	end
end					% end of outer loop.

	if (irho == 0)
		fprintf('unable to calculate min. and max. values of the angle\n');
		fprintf('of rotation, since none of the grid points correspond\n');
		fprintf('to admissible axes of rotation.\n');
		fprintf('Min., Max. axis longitude over grid: %g\t%g\tdegrees', along1, along2)
		fprintf('Min., Max. axis latitude over grid: %g\t%g\tdegrees', alat1, alat2)
		return
	end

	fprintf('Min., Max. angle of rotation over the confidence\n')
	fprintf('region: ", rmind, rmaxd ," degrees. Note: these\n')
	fprintf('values are the min. and max. over a rectangular\n')
	fprintf('grid superimposed on the set of admissible axes.\n')
	fprintf('Min., Max. axis longitude over grid: %g\t%g\tdegrees', along1, along2)
	fprintf('Min., Max. axis latitude over grid: %g\t%g\tdegrees', alat1, alat2)
	fprintf('Grid of %g longitude values\n', nlong)
	fprintf(' (cols.) and %g latitude values (rows)\n', nlat)

	if (ind == 1)
		% Mais escritas????
% 		if (icase < 7)
% 			try
% 				fclose(unit2fid(find(unit2fid(:,1)==31,1,'last'),2));
% 				unit2fid=unit2fid(unit2fid(:,1)~=31,:);
% 			end
% 			try
% 				fclose(unit2fid(find(unit2fid(:,1)==32,1,'last'),2));
% 				unit2fid=unit2fid(unit2fid(:,1)~=32,:);
% 			end
% 			writeFmt(1,format_916,'ufile','nlong','lfile');
% 		elseif (icase == 7)
% 			try
% 				fclose(unit2fid(find(unit2fid(:,1)==31,1,'last'),2));
% 				unit2fid=unit2fid(unit2fid(:,1)~=31,:);
% 			end
% 			writeFmt(1,[format_917],'ufile','nlong');
% 		end
	end

% --------------------------------------   meig   ---------------------------------------
% subroutine to calculate the matrix:
%     m = (azero - 1)*b - dvect*transpose(dvect)
% and to find its eigenvalues and eigenvectors.
% inputs:
%   qf - the matrix qtilde in full storage mode
%   uxx - optimal axis of rotation
% outputs:
%   nu - a 3-vector containing the eigenvalues of m
%   w  - a 3x3 matrix containing (as its columns) the eigenvectors of m
%   mf - the matrix m (in full storage mode)
%   icase - indicator for the type of region
%   jer - an error indicator
%      0 - all is well
%      1 - an error has occurred
function [nu,w,mf,icase,jer] = meig(qf,uxx,nu,w,mf,varargin)

	nrot=[];
	dvect=zeros(1,3);
	b=zeros(3,3);
	d=zeros(1,3);
	z=zeros(3,3);
	work=zeros(1,6);
	kpos=zeros(1,3);

	jer=0;
	azero = qf(1,1);
	for  i=1:3;
		dvect(i) = qf(i+1,1);
	end
	for  i=1:3
		for  j=1:3
			b(i,j)=qf(i+1,j+1);
		end
	end
	for  i=1:3
		for  j=1:3
			mf(i,j) = (azero - 1).*b(i,j) - dvect(i).*dvect(j);
		end
	end
	%[mf,dumvar2,dumvar3,d,z,dumvar6,work,nrot]=jacob2(mf,3,3,d,z,3,work,nrot);
	[z,d] = eig(mf);
	[d, ind] = sort(diag(d));
	z = z(:,ind);

	% check that m has one negative and two positive eigenvalues.
	nneg=0;
	npos=0;
	for  i=1:3
		if (d(i) < 0)
			nneg=nneg + 1;
			kneg=i;
		elseif (d(i) > 0)
			npos=npos + 1;
			kpos(npos)=fix(i);
		end
	end
	if ((nneg ~= 1) || (npos ~= 2))
		jer=1;
		fprintf('error in meig:\n');
		fprintf('it is not the case that m has one negative eigenvalue\n');
		fprintf('and two positive ones.  there is some error, as m\n');
		fprintf('should have eigenvalues of this type.\n');
		return
	end

	% count the no. of positive eigenvalues of m which are greater than azero - 1.
	nlarg=0;
	for  j=1:2
		k=kpos(j);
		if (d(k) >(azero - 1.0))
			nlarg=nlarg + 1;
		end
	end

	if (nlarg == 2)
		% The set of admissible axes consists of two (anti-podal) caps.
		% Set nu(3)=the negative eigenvalue and arrange that the eigen-
		% vector corresp. to nu(3) makes an acute angle with the optimal axis.
		% if the preferred cap does not contain one of the poles, icase=1.
		% otherwise, icase=2 or 3, depending on whether the pole involved is
		% the north or the south.

		nu(3)=d(kneg);
		w(1:3, 3)=z(1:3, kneg);

		for (j=1:2)
			k=kpos(j);
			nu(j)=d(k);
			w(1:3, j) = z(1:3, k);
		end
		dotp=0;
		for (i=1:3)
			dotp = dotp + w(i,3) * uxx(i);
		end
		if (dotp < 0)
			w(1:3, 3) = -w(1:3, 3);
		end

		if (mf(3,3) > (azero - 1))
			icase=1;
			% make sure that the columns of w form a right-handed system.
			det = w(1,1) * w(2,2) * w(3,3) + w(1,2) * w(2,3) * w(3,1) + w(1,3) * w(3,2) * w(2,1) - ...
				  w(1,3) * w(2,2) * w(3,1) - w(1,2) * w(2,1) * w(3,3) - w(1,1) * w(2,3) * w(3,2);
			if (det <= 0)
				w(1:3, 1) = -w(1:3, 1);
			end
		else
			if (w(3,3) >= 0)
				icase=2;
			end
			if (w(3,3) < 0)
				icase=3;
			end
			% make sure that the columns of w form a right-handed system.
			det = w(1,1) * w(2,2) * w(3,3) + w(1,2) * w(2,3) * w(3,1) + w(1,3) * w(3,2) * w(2,1) - ...
				  w(1,3) * w(2,2) * w(3,1) - w(1,2) * w(2,1) * w(3,3) - w(1,1) * w(2,3) * w(3,2);
			if (det <= 0)
				w(1:3, 1) = -w(1:3, 1);
			end
		end
	elseif (nlarg == 1)
		% The set of admissible axes is the complement of two
		% excluded (anti-podal) caps.  set nu(1)=the negative eigenvalue.
		% note: nu(2)<=nu(3), so nu(3) is the positive eigenvalue which
		% is greater than azero - 1.  if the excluded caps contain the
		% poles, icase=4.  otherwise, icase=5.  if icase=4, arrange that
		% the eigenvector corresp. to nu(3) is in the same cap as the north pole.

		nu(1)=d(kneg);
		w(1:3,1) = z(1:3,kneg);
		for (j=2:3)
			k=kpos(j-1);
			nu(j)=d(k);
			w(1:3, j) = z(1:3, k);
		end

		if (mf(3,3) <= (azero - 1))
			icase=5;
			det = w(1,1) * w(2,2) * w(3,3) + w(1,2) * w(2,3) * w(3,1) + w(1,3) * w(3,2) * w(2,1) - ...
				  w(1,3) * w(2,2) * w(3,1) - w(1,2) * w(2,1) * w(3,3) - w(1,1) * w(2,3) * w(3,2);
			if (det <= 0)
				w(1:3, 1) = -w(1:3, 1);
			end
		else
			icase = 4;
			if (w(3,3) < 0)
				for (i=1:3)
					w(i,3)= -w(i,3);
				end
			end
		end
	else
		% nlarg=0.  any axis is admissible.
		icase=6;
	end

% -----------------------------------------   evalf   -----------------------------------------
% subroutine to evaluate the function
%     f(phit) = 1./sqrt((dphi/dphit)**2 + (dtheta/dphit)**2)
% note:  f(phit) = dphit/ds, where s is the arc length along the
% bounding curve (measured in radians).
% additional note:  in the process of evaluating f(phit), evalf
% also calculates phi, theta, and the vector u.
% inputs:
%   phit   - parameter along the bounding curve
%   nu     - array containing the eigenvalues of the matrix m
%   w      - array containing the eigenvectors of the matrix m
%   azero  - qf(1,1)
%   icode  - indicator for mode of usage
%     0    - first point on bounding curve
%     pos. - subsequent points
%   oldu   - the u -vector for the immediately preceding point, in the
%            case icode is positive
%   ophi   - the phi value for the immediately preceding point, in the
%            case icode is positive
% outputs:
%   fval   - value of the function f(phit)
%   phi    - longitude corresp. to phit
%   theta  - latitude corresp. to phit
%   u      - axis of rotation corresp. to phit
%   jer    - error indicator
%     0    - all is well
%     1    - an error has occurred
% phit and thetat are the longitude and latitude, respectively, measured
% in the w-coordinate system (i.e., the coordinate system determined by
% the eigenvectors of m).

function [fval,phi,theta,u,jer] = evalf(phit, nu, w, azero, icode, oldu, ophi)

	eta=zeros(1,3);
	dudpt=zeros(1,3);
	dedpt=zeros(1,3);

	jer=0;
	cospt=cos(phit);
	sinpt=sin(phit);
	denom=nu(1).*cospt.*cospt + nu(2).*sinpt.*sinpt - nu(3);
	num=azero - 1.0 - nu(3);
	costt=sqrt(num./denom);
	thetat=acos(costt);
	eta(1)=costt.*cospt;
	eta(2)=costt.*sinpt;
	eta(3)=sin(thetat);
	for  k=1:3
		u(k)=0;
		for  l=1:3
			u(k)=u(k) + w(k,l) * eta(l);
		end
	end
	lambda=(nu(2) - nu(1))./num;
	dedpt(1)= -eta(2).*(lambda.*eta(1).*eta(1) + 1.);
	dedpt(2)= -eta(1).*(lambda.*eta(2).*eta(2) - 1.);
	ss=eta(1).*eta(1) + eta(2).*eta(2);
	dedpt(3)=lambda.*eta(1).*eta(2).*ss./eta(3);
	for  k=1:3
		dudpt(k)=0;
		for  l=1:3
			dudpt(k)=dudpt(k) + w(k,l).*dedpt(l);
		end
	end
	cost2=u(1).*u(1) + u(2).*u(2);
	dpdpt=( -u(2).*dudpt(1) + u(1).*dudpt(2))./cost2;
	dtdpt=dudpt(3)./sqrt(cost2);
	fval = 1./sqrt(dpdpt.*dpdpt + dtdpt.*dtdpt);
	phi=atan2(u(2),u(1));
	theta=asin(u(3));
	% modify phi if necessary.
	if (icode == 0)
		return
	end
	alpha=u(1).*oldu(1) + u(2).*oldu(2);
	beta= -u(1).*oldu(2) + u(2).*oldu(1);
	delphi=atan2(beta,alpha);
	ephi=ophi + delphi;
	if (abs(ephi - phi) < pi/180)
		return
	elseif (abs(ephi - phi - 2*pi) < pi/180)
		phi=phi + 2*pi;
	elseif (abs(ephi - phi + 2*pi) < pi/180)
		phi=phi - 2*pi;
	else
		jer=1;
		fprintf('error in evalf: unable to determine phi, for a\n');
		fprintf('point on the bounding curve.\n');
	end


% largin(x) = the largest integer which is less than or equal to x.
function [larginresult,x]=largin(x,varargin)
	larginresult=fix(x);
	if (x < 0)
		r = x - (larginresult);
		if (r < 0)
			larginresult=larginresult -1;
		end
	end

function in1=strAssign(in1,r1,r2,rhs)
	if isempty(r1),  r1=1;  end
	if isempty(r2),  r2=length(in1);  end
	rhsL=length(rhs);   lhsL=r2-r1+1;
	% rhs is too long for the left
	if lhsL<rhsL,  in1(r1:r2)=rhs(1:lhsL);  end
	% sizes match
	if lhsL==rhsL,  in1(r1:r2)=rhs;  end
	% rhs is too short for the left
	if lhsL>rhsL,  in1(r1:r1+rhsL-1)=rhs;  in1(r1+rhsL:lhsL)=' ';  end

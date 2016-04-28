function wts=bpass3d(nnx,nny,dx,dy,wlong,wshort)
% BPASS3D set up bandpass filter weights in 2 dimensions
%
% Usage:  wts3d=bpass3d(nnx,nny,dx,dy,wlong,wshort);
%
% Example:
%           nnx=64;nny=64;dx=0.5;dy=0.5;wlong=12;wshort=1;
%           wts=bpass3d(nnx,nny,dx,dy,wlong,wshort);
%           surf(wts);view(-30,65);
%           title('Bandpass filter in fourier domain');
%
%  Maurice A. Tivey     March 1996
%  Joaquim Luis         May   2004
%       Cleaned a bit of the code for saving RAM and removed the ploting calls
%---------------------------------------------------------------------------

% $Id: bpass3d.m 4247 2013-12-27 03:50:16Z j $

	twopi = pi*2;
	dk1 = twopi/((nnx-1)*dx);
	dk2 = twopi/((nny-1)*dy);

	% calculate wavenumber array
	% nx2 = nnx/2;            ny2 = nny/2;          % And if nnx|nny are odd?
	% dkx = pi/(nnx*dx);      dky = pi/(nny*dy);
	% kx = (-nx2:nx2-1).*dkx; ky = (-ny2:ny2-1).*dky;
	% X = ones(size(ky))'*kx; Y = ky'*ones(size(kx));
	% k = fftshift(2*sqrt(X.^2+Y.^2));  % wavenumber array
	nx2 = fix(nnx/2);       ny2 = fix(nny/2);
	if (rem(nnx,2) == 0),	sft_x = 1;
	else					sft_x = 0;
	end
	if (rem(nny,2) == 0),	sft_y = 1;
	else					sft_y = 0;
	end
	dkx = 2*pi / (nnx*dx);  dky = 2*pi / (nny*dy);
	kx = (-nx2:nx2-sft_x).*dkx;     ky = (-ny2:ny2-sft_y).*dky;
	X = repmat(kx,length(ky),1);    Y = repmat(ky',1,length(kx));
	k = ifftshift(sqrt(X.^2+Y.^2));      % wavenumber array
	clear X Y;

	if (wshort == 0),	wshort = max(dx*2,dy*2); end
	if (wlong == 0),	wlong = min(nnx*dx,nny*dy); end

	klo = twopi/wlong;      khi = twopi/wshort;
	khif = 0.5*khi;         klof = 2*klo;
	dkl = klof-klo;         dkh = khi-khif;

	% See if a window for verbose exists
	old_show = get(0,'ShowHiddenHandles');
	set(0,'ShowHiddenHandles','on')
	h = get(0,'Children'); 
	figdmsg = findobj(h,'flat','tag','Wdmsgfig');
	set(0,'ShowHiddenHandles',old_show)

	if (isempty(figdmsg))
		try,    message_win('create',' BPASS3D SET UP BANDPASS WEIGHTS ARRAY:','width',400,'height',300,'edit','y');
		catch,  fprintf(' BPASS3D\n SET UP BANDPASS WEIGHTS ARRAY :\n')
		end
	else
		try,    message_win('add',' BPASS3D SET UP BANDPASS WEIGHTS ARRAY:')
		catch,  fprintf(' BPASS3D\n SET UP BANDPASS WEIGHTS ARRAY :\n')
		end
	end
	try,    message_win('add',sprintf(' HIPASS COSINE TAPER FROM K= %10.6f TO K= %10.6f',klo,klof))
	catch,  fprintf(' HIPASS COSINE TAPER FROM K= %10.6f TO K= %10.6f\n',klo,klof)
	end
	try,    message_win('add',sprintf(' LOPASS COSINE TAPER FROM K= %10.6f TO K= %10.6f',khif,khi))
	catch,  fprintf(' LOPASS COSINE TAPER FROM K= %10.6f TO K= %10.6f\n',khif,khi)
	end
	try,    message_win('add',sprintf(' DK1,DK2= %10.4f  %10.4f',dk1,dk2))
	catch,  fprintf(' DK1,DK2= %10.4f  %10.4f\n',dk1,dk2)
	end

	wl1 = 1000;     wl2 = 1000;
	if (klo > 0)    wl1 = twopi/klo; end
	if (klof > 0)   wl2 = twopi/klof; end
	wl3 = twopi/khif;
	wl4 = twopi/khi;
	wnx = twopi/(dk1*(nnx-1)/2);
	wny = twopi/(dk2*(nny-1)/2);

	try,    message_win('add','IE BANDPASS OVER WAVELENGTHS')
	catch,  fprintf('IE BANDPASS OVER WAVELENGTHS\n')
	end
	try,    message_win('add',sprintf('   INF CUT-- %8.3f --TAPER-- %8.3f (PASS) %8.3f --TAPER--%8.3f',wl1,wl2,wl3,wl4))
	catch,  fprintf('   INF CUT-- %8.3f --TAPER-- %8.3f (PASS) %8.3f --TAPER--%8.3f\n',wl1,wl2,wl3,wl4)
	end
	try,    message_win('add',sprintf('   --  CUT TO NYQUIST X,Y= %8.3f  %8.3f',wnx,wny))
	catch,  fprintf('   --  CUT TO NYQUIST X,Y= %8.3f  %8.3f\n',wnx,wny)
	end

	%nnx2 = nnx/2+1;     nny2 = nny/2+1;
	wts = zeros(size(k));  % initialise to zero
	for i=1:nny,
		for j=1:nnx,
			if (k(i,j) > klo) 
				if (k(i,j) < khi)    wts(i,j) = 1;   end
			end
		end
	end
	for (i = 1:nny)
		for (j = 1:nnx)
			if (k(i,j) > klo) 
				if (k(i,j) < klof),		wts(i,j) = wts(i,j)*(1-cos(pi*(k(i,j)-klo)/dkl))/2; end
			end
			if (k(i,j) > khif) 
				if (k(i,j) < khi),		wts(i,j) = wts(i,j)*(1-cos(pi*(khi-k(i,j))/dkh))/2; end
			end
		end
	end

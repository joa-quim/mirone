function [fractfield, hdrStruct] = gen_UMF2d(alpha, C1, H, dim)
% Generation of a stocastic 2D multi-fractal random field
% Can be used to simulate DEM, turbulent fields (clouds etc)
%
% Based on:
% Lovejoy - Schertzer "Nonlinear variability in geophysics: multifractal simulations and analysis"
% 
% Inputs:
% alpha = [0..2] Levy parameter, governs the mutli-fractal behviour 
%    (0 = monofractals)
% C1 = governs the sparseness of the field
% H = fractional integration, parameters, governs the ruggedness or smoothness of the field
% dim = size of the fractal field

%   Author(s): Andrea Monti Guarnieri, P. Biancardi, , D. D'Aria et. al

% Example - syntehsis of a likely Digital Elevation Model
% H=1.9; dem = gen_UMF2d(1.8, 0.05, H, 150); 
% colormap('bone'); surfl(dem); shading('interp')
% H governs the ruggedness

% Copyright (c) 2010, Andrea monti guarnieri
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%       
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.

% Cleaned a lot of the unfortunatelly usual memory absurdity use

% $Id: gen_UMF2d.m 4728 2015-07-22 15:14:51Z j $

	if nargin == 0
		alpha = 1.8;	C1 = 0.05;		H = 0.9;	dim = 512;
	end

	% Pseudo corrdinates of this fractal field
	w = -12;	e = -2;		s = 35;		n = 45;

	d = 2;

	alphap = alpha / (alpha-1);

	% Generate Levy distrbuted noise
	noise = Salpha(alpha,dim,dim);
	noise = 2.3 * (C1/(alpha-1))^(1/alpha) * noise;

	% Noise transformation
	NOISE   = fft2(noise);
	asse    = assefr(1,dim)*2*pi;
	[kx ky] = meshgrid(asse);

	% Compute the "ultiscaling behaviour" filter
	FILTER = fftshift(abs(kx+1i*ky).^(-d/alphap));
	FILTER(1,1) = 0;

	% Lowpass filter
	FILTER2 = lambdafilt2d(dim,dim/2,.1);

	% Generate the field
	noise     = real(ifft2(NOISE.*FILTER.*FILTER2));
	clear FILTER2 NOISE
	noisetemp = noise-(C1/(alpha-1)*log(dim*dim));
	noisetemp = exp(noisetemp);
	FIELD     = fft2(noisetemp);
	clear noisetemp

	% Filter
	FILTER  = fftshift(abs(kx+1i*ky).^(-H));
	clear kx ky noise
	FILTER(1,1) = dim^(H);	%/8;
	fract  = single(real(ifft2(FIELD.*FILTER)));
	
	tmp.X = linspace(w,e,dim);
	tmp.Y = linspace(s,n,dim);
	tmp.head = [w e s n double(min(fract(:))) double(max(fract(:))) 0 (e-w)/(dim-1) (n-s)/(dim-1)];
	tmp.name = 'Fractal Surface';

	fract = c_grdtrend(fract,tmp.head,'-D','-N1');		% Remove mean
	tmp.head(5) = min(fract(:));	tmp.head(6) = max(fract(:));

	if (nargout == 0)
		mirone(fract, tmp)
	else
		hdrStruct = tmp;
		fractfield = fract;
	end

% ------------------------------------------------------------
function x=Salpha(alpha,dim1,dim2)
% Generation of a white noise with Levy pdf (parameter alpha)  
% rand('state',3);

	phi0 = -pi/2 * (1-abs(1-alpha)) / alpha;
	phi = rand(dim1,dim2) * pi - pi / 2;			% phi uniform in -pi/2 pi/2

	% W as an exponential, averge is 1
	a = rand(dim1, dim2);
	W = -log(1-a);
	x = (sin(alpha*(phi-phi0))./(cos(phi)).^(1/alpha) .* ((cos(phi-alpha*(phi-phi0))./W).^((1-alpha)/alpha)));

% ------------------------------------------------------------
function  fi = assefr(Fs,n)
	if rem(n,2)==0   
		fi=((0:n-1)-n/2)*Fs/n;      %n even   
	else
		fi=((0:n-1)-(n-1)/2)*Fs/n;  %n odd
	end

% ------------------------------------------------------------
function filtro=lambdafilt2d(dim,raggio,tau)
% Compute a low-pass filter unitary on a circle 'raggio'
% and exponential decay outside with time constant tau

	% crea un asse delle frequenze(dim pari o dispari)
	f = assefr(1,dim)*dim;
	[x,y] = meshgrid(f);

	filtro = sqrt(x.^2+y.^2);

	bool   = (filtro > raggio);
	filtro = -((filtro-raggio).*bool);

	filtro = exp(tau*filtro);
	filtro = ifftshift(filtro);

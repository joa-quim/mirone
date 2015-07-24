function [w, to_restore] = mboard(w,nx,ny,nnx,nny,mode)
% Pad a matrix in a civilized way (not just by adding zeros like Matlab funcs do) before FFTit
%
% W = MBOARD(W,NX,NY) mirror the matrix about last row and column
% [W,TO_RESTORE] = MBOARD(W,NX,NY,NNX,NNY) taper the matrix with a NNX, NNY hanning window.
% That means the matrix will be added NNY/2 rows on top; NNY/2 rows at bottom;
% NNX/2 at left and NNY/2 at right. The pading skirt will fall down to zero on each of these bands
% TO_RESTORE is a vector with the width of the pading bands [top, bot, left, right] that
% is used to restore the original matrix dimensions after the FFT
%
% If NNX or NNY == 0, than those are estimated as being the closest number to NX * 1.2 (or NY)
% IF MBOARD([],NX,NY,0,0) compute only the good NNX = NX * 1.2 & NNY = NY * 1.2 and return
% them in W. TO_RESTORE will contain "nlist". Note that there are no error testing. 

%	Copyright (c) 2004-2013 by J. Luis
%
% 	This program is part of Mirone and is free software; you can redistribute
% 	it and/or modify it under the terms of the GNU Lesser General Public
% 	License as published by the Free Software Foundation; either
% 	version 2.1 of the License, or any later version.
% 
% 	This program is distributed in the hope that it will be useful,
% 	but WITHOUT ANY WARRANTY; without even the implied warranty of
% 	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% 	Lesser General Public License for more details.
%
%	Contact info: w3.ualg.pt/~jluis/mirone
% --------------------------------------------------------------------

% $Id$

if (nargin == 5),    mode = 'taper';     end
if (nargin == 3),    mode = 'mirror';     end

% List of numbers with prime factors <= 5
nlist = {64,72,75,80,81,90,96,100,108,120,125,128,135,144,150,160,162,180,192,200,...
    216,225,240,243,250,256,270,288,300,320,324,360,375,384,400,405,432,450,480,...
    486,500,512,540,576,600,625,640,648,675,720,729,750,768,800,810,864,900,960,...
    972,1000,1024,1080,1125,1152,1200,1215,1250,1280,1296,1350,1440,1458,1500,...
    1536,1600,1620,1728,1800,1875,1920,1944,2000,2025,2048,2160,2187,2250,2304,...
    2400,2430,2500,2560,2592,2700,2880,2916,3000,3072,3125,3200,3240,3375,3456,...
    3600,3645,3750,3840,3888,4000,4096,4320,4374,4500,4608,4800,4860,5000}';

if (nargin >= 4 && nnx == 0)               % Find the good number that is about 20% larger than nx
    nnx = round(nx*1.2);
    if (nnx > nlist{end})   % Ohps, the list is too short
        nnx = nx;
    else
        nnx = cat(2,nlist{:}) - nnx;
        id = find(nnx > 0);
        nnx = nlist{id(1)};
    end
end
if (nargin >= 5 && nny == 0)               % Find the good number that is about 20% larger than ny
    nny = round(ny*1.2);
    if (nny > nlist{end})
        nny = ny;
    else
        nny = cat(2,nlist{:}) - nny;
        id = find(nny > 0);
        nny = nlist{id(1)};
    end
end

if (isempty(w))
    w = [nnx nny];    to_restore = nlist;    return
end

if (strcmp(mode,'mirror'))
    w(ny+1:2*ny,1:nx) = w(ny:-1:1,:);
    w(:,nx+1:2*nx) = w(:,nx:-1:1);
    to_restore = [];        % The compiler hates non initialized outputs
else    % pad with a hanning window
	dnx = nnx - nx;         dny = nny - ny;
    dnx_w = fix(dnx/2);     dnx_e = dnx - dnx_w;    % Need to do this if dnx is odd
    dny_n = fix(dny/2);     dny_s = dny - dny_n;    % Need to do this if dny is odd  
    to_restore = [dny_n dny_s dnx_w dnx_e];
    % Extend to South
	vhan = hanning(2*dny_s);    vhan = vhan(end/2+1:end);
	tmp1 = repmat(vhan,1,nx);       % Replicate the hanning vector by array's n_column
	tmp2 = repmat(w(end,:),dny_s,1);% Replicate the array's last row by the y_s pading size
    if (isa(w,'single')),     tmp12 = single(double(tmp1) .* double(tmp2));
    else                        tmp12 = tmp1 .* tmp2;   end
	w(ny+1:ny+dny_s,:) = tmp12;
    % Extend to East
	vhan = hanning(2*dnx_e);     vhan = vhan(end/2+1:end)';
	tmp1 = repmat(vhan,ny+dny_s,1); % Replicate the hanning vector by array's extended n_rows
	tmp2 = repmat(w(:,nx),1,dnx_e); % Replicate the array last column by the x_e pading size
    if (isa(w,'single')),     tmp12 = single(double(tmp1) .* double(tmp2));
    else                        tmp12 = tmp1 .* tmp2;   end
	w(:,nx+1:nx+dnx_e) = tmp12;
    % Extend to North
	vhan = hanning(2*dny_n);    vhan = vhan(1:end/2);
	tmp1 = repmat(vhan,1,nx+dnx_e); % Replicate the hanning vector by array's new n_column
	tmp2 = repmat(w(1,:),dny_n,1);  % Replicate the array's first row by the y_n pading size
    if (isa(w,'single')),     band_n = single(double(tmp1) .* double(tmp2));
    else                        band_n = tmp1 .* tmp2;   end
    w = [band_n; w];                % Remember that we are now adding before the original first row
    % Extend to West
	vhan = hanning(2*dnx_w);    vhan = vhan(1:end/2)';
	tmp1 = repmat(vhan,nny,1);      % Replicate the hanning vector by array's new n_row (nny)
	tmp2 = repmat(w(:,1),1,dnx_w);  % Replicate the array's first column by the x_w pading size
    if (isa(w,'single')),      band_w = single(double(tmp1) .* double(tmp2));
    else                        band_w = tmp1 .* tmp2;   end
    w = [band_w w];                 % Remember that we are now adding before the original first column
end

%---------------------------------------------------------------------
function w = hanning(n)
% HANNING Returns a symmetric N point hanning window 
a0 = 0.5;   a1 = 0.5;   a2 = 0;     a3 = 0;     a4 = 0;
if (~rem(n,2))  % Even length window
    half = n/2;
    x = (0:half-1)'/(n-1);
    w = a0 - a1*cos(2*pi*x) + a2*cos(4*pi*x) - a3*cos(6*pi*x) + a4*cos(8*pi*x);
    w = [w; w(end:-1:1)];
else            % Odd length window
    half = (n+1)/2;
    x = (0:half-1)'/(n-1);
    w = a0 - a1*cos(2*pi*x) + a2*cos(4*pi*x) - a3*cos(6*pi*x) + a4*cos(8*pi*x);
    w = [w; w(end-1:-1:1)];
end

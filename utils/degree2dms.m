function out = degree2dms(deg,format,n,dec_str)
% DEGREE2DMS Converts angles from degrees to deg, min, sec
% FORMAT may be one the following strings:
%   DDMM, DDMM.x, DDMMSS, DDMMSS.x
% where DD stands for degrees, MM for minutes, SS for seconds and .x means the precedent field
% (e.g. MM or SS) will be writen with the number of decimal figures as determined by N (default is 2)
% Note: this only applyies to char string outputs (see below). In case of numeric, figures are in double precision
% DEC_STR chooses between output as numeric (dec_str = 'numeric'), or as char strings (dec_str = 'str')
% Default is 'numeric'.
% The output comes into a structure OUT with the number of fields determined by FORMAT
% (e.g) out.dd out.mm out.ss if FORMAT = 'DDMMSS'
%
% Example:
%   out = degree2dms(10.3533,'DDMMSS.x',1,'str')
%    dd: '10'
%    mm: '21'
%    ss: '11.9'

%	Copyright (c) 2004-2012 by J. Luis
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

	if (nargin == 0)
		error('ERROR in DEGREE2DMS: no angle provided. You must think I''m bruxo')
	elseif (nargin == 1)
		format = 'DDMM.x';      n = '%.2f';      dec_str = 'numeric';
	elseif (nargin == 2) && ischar(format)		% not tested if it's a valid format string
		n = '%.2f';             dec_str = 'numeric';
	elseif (nargin == 3)						% not tested if n is numeric
		dec_str = 'numeric';
	end
	if isempty(deg)
		out = [];   return;
	end

	if isempty(format),     format = 'DDMM.x';      end
	if isempty(n)
		n = '%.2f';     n_decimal = 4;
	else
		n_decimal = n;
		n = ['%.' num2str(fix(n)) 'f'];			% Construct a C format string because the num2str(x,precision) is shit.
	end											% Half of the times it doesn't obey to the precision instruction.

	% Construct angle sign vector and also make sure it's = 1 when deg == 0
	ang_sign = sign(deg);   ang_sign = ang_sign + (ang_sign == 0);

	% Decompose the angle 'deg' into degrees, minutes and seconds
	deg = abs(deg);     dd = fix(deg);      ms = 60*(deg - dd);     mm = fix(ms);   ss = 60*(ms - mm);
	ss = round2n(ss,n_decimal);
	% To avoid the shaming case of (it happened to me) -9:29:60
	ind = find(ss >= 60);
	if ~isempty(ind);   mm(ind) = mm(ind) + 1;   ss(ind) = 0;   end
	% And now the minutes case
	ind = find(mm >= 60);
	if ~isempty(ind);   dd(ind) = dd(ind) + 1;   mm(ind) = 0;   end
	if (strcmp(format,'DDMM.x'))
		mm = round2n(mm+ss/60,n_decimal);
		ind = find(mm >= 60);
		if ~isempty(ind);   dd(ind) = dd(ind) + 1;   mm(ind) = 0;   end
	elseif (strcmp(format,'DDMM'))
		mm = round(mm+ss/60);
		ind = find(mm >= 60);
		if ~isempty(ind);   dd(ind) = dd(ind) + 1;   mm(ind) = 0;   end
	end

	% Store the angle sign into the largest of dd, mm or ss. Needed when (for example) dd = 0.
	dd_sign = ang_sign .* (dd ~= 0);        mm_sign = ang_sign .* (dd == 0 & mm ~= 0);
	ss_sign = ang_sign .* (dd == 0 & mm == 0 & ss ~= 0);

	% Now apply the signs to dd, mm & ss
	dd = ((dd_sign ==0 ) + dd_sign).*dd;    mm = ((mm_sign == 0) + mm_sign).*mm;
	ss = ((ss_sign == 0) + ss_sign).*ss;

	% And finally compute the output
	switch format
		case 'DDMM'
			if strcmp(dec_str,'numeric')
				out.dd = dd;            out.mm = mm;
			else    % string output
				out.dd = num2str(dd,'%02d');   out.mm = num2str(mm,'%02d');
			end
		case 'DDMM.x'
			if strcmp(dec_str,'numeric')
				out.dd = dd;            out.mm = mm;
			else    % string output
				out.dd = num2str(dd);   out.mm = num2str(mm,n);
			end
		case 'DDMMSS'
			if strcmp(dec_str,'numeric')
				out.dd = dd;            out.mm = mm;            out.ss = round(ss);
			else    % string output
				out.dd = num2str(dd);   out.mm = num2str(mm,'%02d');   out.ss = num2str(round(ss),'%02d');
			end
		case 'DDMMSS.x'
			if strcmp(dec_str,'numeric')
				out.dd = dd;            out.mm = mm;            out.ss = ss;
			else    % string output
				out.dd = num2str(dd);   out.mm = num2str(mm,'%02d');   out.ss = num2str(ss,n);
			end
	end

function x = round2n(x,n)
% Rounds input data at specified power of 10
	fact = 10 ^ (fix(n));
	x = round(x * fact) / fact;

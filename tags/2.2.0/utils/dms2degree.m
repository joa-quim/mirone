function deg = dms2degree(dd, mm, ss)
% DMS2DEGREE2 Converts angles from deg, min, sec to decimal degrees
% If only two inputs are provided, it is interpreted as deg, min

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

	% Do some error tests
	if (nargin < 2 || nargin > 3)
		errordlg('DMS2DEGREE: Incorrect number of arguments','Error')
		deg = [];   return;
	end

	msg = [];
	if (nargin == 2)
		if ~isequal(size(dd),size(mm));    msg = 'Vectors dd, mm, ss are not of same size';     end
		if any(mm < 0 & dd < 0)
			msg = 'Minus sign present in the mm (minutes) vector';
		end
	else
		if ~isequal(size(dd),size(mm),size(ss));    msg = 'Vectors dd & mm are not of same size';     end
		if any((ss < 0 & mm < 0) | (ss < 0 & dd < 0) | (mm < 0 & dd < 0))
			msg = 'Minus sign present in the mm (minutes) and/or ss (seconds) vector(s)';
		end
	end
	if (~isempty(msg))
		errordlg(['DMS2DEGREE: ' msg],'Error')
		deg = [];   return;
	end

	% Store the angle sign into the largest of dd, mm or ss. Needed when (for example) dd = 0.
	if (nargin == 2)
		neg = ((dd < 0) | (mm < 0));    sinal = ~neg - neg;
		deg  = sinal.*(abs(dd) + abs(mm)/60);
	else
		neg = ((dd < 0) | (mm < 0) | (ss < 0));    sinal = ~neg - neg;
		deg  = sinal.*(abs(dd) + abs(mm)/60 + abs(ss)/3600);
	end

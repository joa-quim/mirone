function [mm, dd] = doy(year, ddd, opt)
% Convert Day Of Year into month and day or vice-versa depending on the arguments
%
% DY = DOY('datevec')		% returns Day Of Year corresponding to date in DATEVEC
% DY = DOY(YEAR, MM, DD])	% returns Day Of Year corresponding to the date YYYY MM DD
% DY = DOY([YEAR MM DD])	% Same as above but with date in a vector
%
% Examples:
%		dy = doy('05-Aug-2014')     (DOY = 217)
%       dy = doy([2014 8 5])        (DOY = 217)
%       dy = doy(2014, 8, 5)        (DOY = 217)
%  [mm,dd] = doy(2014, 217)         (mm = 8; dd = 5)

%	Copyright (c) 2004-2014 by J. Luis
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

% $Id: doy.m 4614 2014-11-03 01:02:31Z j $

	if ((nargin == 1 && (isa(year, 'char')) || numel(year) == 3) || (nargin == 3))
		to_doy = true;		% Get DOY from a DATEVEC date string or [year month day]
	elseif (nargin > 1)
		to_doy = false;		% From DOY get [monthm, day]
	else
		error('DOY:Bad usage. Wrong number of arguments or their type')
	end

	if (to_doy)			% Get DOY
		if (isa(year, 'char'))
			d = datevec(year);
			mm = datenum([d(1:3), 0, 0, 0]) - datenum([d(1), 1, zeros(1, 4)]);
		elseif (numel(year) == 3)
			mm = datenum([year(1:3), 0, 0, 0]) - datenum([year(1), 1, zeros(1, 4)]);
		else
			mm = datenum([year, ddd, opt, 0, 0, 0]) - datenum([year, 1, zeros(1, 4)]);			
		end

	else				% Get month and day from DOY
		v = datevec(datenum(year, 1, ddd));
		mm = v(2);		dd = v(3);
	end

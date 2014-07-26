function [mm,dd] = doy(year, ddd, opt)
% Convert Day Of Year into month and day or vice-versa depending on the value of OPT
%
% OPT if not provided or equal to -1 ...

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

% $Id$

	if (nargin ~= 2 && nargin ~= 3)
		error('DOY:Wrong number of input args')
	elseif (nargin == 2)
		opt = -1;
	end

	if (opt == 1)		% Get DOY
		error('DOY:option not yet programmed')
	else				% Get month and day from DOY
		v = datevec(datenum(year, 1, ddd));
		mm = v(2);		dd = v(3);
	end
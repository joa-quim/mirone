function val = test_dms(s)
% Test if the input text string S is in the form dd:mm or dd:mm:ss
% val = test_dms(s) returns a cell array with the individual tokens
% if an error ocurrs or an empty imput is given, it returns an empty array

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

	if isempty(s),    val = [];    return;   end
	error = 0;
	if (strcmp(s(1),':') || strcmp(s(end),':') || strfind(s,'::'))
		errordlg('One or more of the dd, mm or ss fields is empty','Error')
		error = error + 1;
	end
	[t,r]=strtok(s,':');
	unit{1} = t;
	if isnan(str2double(unit{1}))
		str = [unit{1} '   is not a valid number'];
		error = error + 1;
		errordlg(str,'Error')
	end
	i = 2;
	val = cell(1);
	unit = cell(1);
	while (~isempty(r))
		[t,r] = strtok(r,':');
		unit{i} = t;
		if isnan(str2double(unit{i}))
			str = [unit{i} '   is not a valid number'];
			error = error + 1;
			errordlg(str,'Error')
		end
		i = i + 1;
	end
	if (~error), val = unit;	end

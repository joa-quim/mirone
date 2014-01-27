function out = check_cmop(fname)
% Check if the FNAME file is a CMOP data in csv format.
%
% These files come with time in ascii and convering with datenum is rather slow

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

	fid = fopen(fname);
	str = fread(fid,'*char');
	fclose(fid);
	str = str';
	str(str == '''') = '';			% Remove the annoying "'"

	t = str(1:300);					% Just a very small chunk to carrie on some guess work
	ind = strfind(t, char(10));		% Find the line breaks

	% Example after the "'" having been removed
	%time PST, saturn04.860.R.Oxygen [%],time PST, saturn04.30.F.Oxygen [%]
	%2012-09-01 00:01:00,87.91153935,2012-09-01 00:11:00,89.47921692

	if ~(strcmp(t(1:4),'time') && t(ind(1)+5) == '-' && t(ind(1)+20) == ',')
		out = true;
		return
	end

	year  = t(ind(1)+1:ind(1)+4);
	n_var = numel(strfind(t(ind(1)+1:ind(2)), year));
	frmt  = repmat('%s%f',1,n_var);

	t = cell(1, n_var * 2);
	[t{1:n_var*2}] = strread(str(ind(1)+1:end),frmt,'delimiter',',');

	for (k = 1:2:n_var*2)
		t1 = t{k};		t2 = t{k+1};
		n = numel(t1);
		if (k > 1)					% Second and on variables are not guarantied to have the same number of points
			while (t1{n} == '0')	% Loop from end up to search for 0s, that mean no data
				n = n - 1;
			end
			t1(n+1:end) = [];		t2(n+1:end) = [];
		end
		t{k} = DateStr2Num(t1, 31);
		t{k+1} = t2;
	end

	hf = ecran(t{1}, t{2}, 'CMOP file');
	if (n_var > 1)
		handEcran = guidata(hf);
		hLine = findobj(handEcran.figure1, 'type', 'line');
		set(hLine,'tag','CMOP')			% This tag will be used by ecran to set function handles to do cross scatter plots
		for (k = 3:2:n_var*2)
			ecran('add',t{3}, t{k+1});
			if (k == 3)
				t = findobj(handEcran.figure1, 'type', 'line');		% Fish the second line handle
				if (t(1) == hLine),		hLine(2) = t(2);
				else					hLine(2) = t(1);
				end
				set(hLine(2),'tag','CMOP')
			end
		end
	end

	h = findobj(hf,'Tag','add_uictx');
	cb = get(h, 'Call');
	feval(cb, h, guidata(hf))		% Call the ecran's add_uictx_CB function

	out = false;		% Means, no error

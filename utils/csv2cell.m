function data = csv2cell(fname)
% CSV2CELL - parses a CSV file into an NxM cell array, where N is the number
% of lines in file and M the number of fields in the longest.
%
% All fields from the CSV file are returned as strings. If the file contains lines
% with different numbers of fields, the empty fields appear as empties, [], in the output data.

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

% $Id: csv2cell.m 4437 2014-05-01 01:29:10Z j $

	fid = fopen(fname);
	text = fread(fid,inf,'*char')';
	fclose(fid);

	ind_nl = [0 regexp(text,'(\r\n|[\r\n])')];
	data = cell(numel(ind_nl)-1,1);
	for (n = 1:numel(ind_nl)-1)						% Loop over all lines in this file
		line = text(ind_nl(n)+1:ind_nl(n+1)-1);		% Extract the current line to simplify ops
		ind_c = [0 strfind(line,',') numel(line)+1];% First & last are to easy up the algo
		ind_q = strfind(line,'"');
		ind = false(1, numel(ind_c));
		if (rem(numel(ind_q), 2)),	ind_q(end) = [];	end		% If '"' are not paired, remove last (to not error down)
		for (k = 1:2:numel(ind_q))					% Find all commas inside pairs of double quotes
			ind = ind | ((ind_c > ind_q(k)) & (ind_c < ind_q(k+1)));
		end
		ind_c(ind) = [];							% and remove their reference
		for (k = 1:numel(ind_c)-1)					% Loop over all fields in this line
			str = line(ind_c(k)+1:ind_c(k+1)-1);
			if (~isempty(str) && str(1) == '"' && str(end) == '"')
				str = str(2:end-1);
			end
			data{n,k} = str;
		end
	end

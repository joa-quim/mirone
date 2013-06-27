	function [bin,n_column,multi_seg,n_headers,isGSHHS] = guess_file(fiche, opt1, opt2)
	% Guess characteristics of file "fiche"
	%
	% Input:
	%	OPT1, if given, will be MAXCHARS
	%	OPT2, if given, will be nl_max
	%
	% Output
	%	BIN        0 or 1, if file is ascii or binary, or empty if error loading file.
	%	N_COLUMN   Number of columns of the data section
	%	MULTI_SEG  0 or 1 depending if file is of the GMT multi-segment type or not
	%	N_HEADERS  Number of headers in file
	%	ISGSHHS    TRUE or FALSE depending if file is a GMT "GSHHS Master File" or not
	%
	% If it detects that "fiche" is ascii this function tries to find out wether the multisegment
	% symbol (">") is present, the number of columns in the file and if it has header lines.
	% NOTE that this last tests may not always give reliable results.

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

	% Error testing
	bin = 0;    multi_seg = 0;  n_headers = 0;  n_column = 0;
	n_args = nargin;
	if (~n_args)
		errordlg('function guess_file: must give an input file name','File Error')
		return
	elseif (n_args == 1)
		MAXCHARS = 4096;		% Maximum characters to load from file
		nl_max = 100;			% Maximum number of file lines to use in tests
	elseif (n_args == 2)
		MAXCHARS = opt1;
	elseif (n_args == 3)
		MAXCHARS = opt1;
		nl_max = opt2;
	end

	% Now remove any leading white space in file name "fiche"
	fiche = ddewhite(fiche);

	if (exist(fiche,'file') ~= 2)
		errordlg(['function guess_file: file "' fiche '" does not exist'],'File Error')
		return
	end

	% ----------------------------------------------------------------------------------------
	fid = fopen(fiche);
	[str,count] = fread(fid,MAXCHARS,'*char');
	fclose(fid);
	str = str';

	A = double(str);
	if (any(A > 126 & A < 192))		% Binary files have bytes with values greater than 126 (but so is the ç char)
		bin = guess_in_bin(fiche);
		return
	end

	str = strread(str,'%s','delimiter','\n');
	if (isempty(str)),		bin = [];	return,		end		% Will trigger a 'Don't know' message
	if (count == MAXCHARS)			% In these cases last line is normally incomplete
		str(end) = [];
	end

	nl_max = min(nl_max,numel(str));
	if (nl_max == 0),	bin = [];	return,		end

	% Make a crude test to find number of columns and number of headers
	delimiters = [9:13 32 44];	% White space characters plus comma
	n_col(1:nl_max) = 0;    n_multi = 0;
	idM = false(1,nl_max);
	for (i = 1:nl_max)
		if isempty(str{i});  continue,   end		% Jump blank lines
		if (str{i}(1) == '#' || str{i}(1) == '%'),		continue,   end		% Jump known comment lines
		if strfind(str{i}(1:min(2,length(str{i}))),'>')   % multisegmet line, so the rest of it is of no interest (but count them)
			n_multi = n_multi + 1;
			idM(i) = true;      % Tag it to deletion
			continue
		end
		str{i} = deblank(str{i});          % Blanks make a mess to the guessing code
		[tok,rem] = strtok(str{i}, delimiters);
		if ~isempty(rem);   n_col(i) = n_col(i) + 1;    end     % count first column
		while ~isempty(rem)
			[tok,rem]=strtok(rem, delimiters);
			n_col(i) = n_col(i) + 1;
		end
	end
	multi_seg = n_multi;
	str(idM) = [];			% idM is the index of the multisegment lines
	n_col(idM) = [];
	nl_max = min(nl_max,numel(str));	% Update because str may have become empty

	% Now decide how many columns have the data lines. The easeast is to assume that the info is in the last line
	% However this may fail if last line contains, for example, the multisegment symbol (">").
	% So, do another test.
	n_col(~n_col) = [];				% Remove zeros which correspond to comment lines (the ones starting by '#') 
	m = min(nl_max,numel(n_col));
	n_headers_candidate = max(0, nl_max - numel(n_col));	% This value will be confirmed later (I'm messing this up)

	if (m > 1)
		n_c1 = n_col(m);	n_c2 = n_col(m-1);		% With bad luck n_c2 can be zero
		if (n_c2 && n_c1 ~= n_c2 && isempty(find(str{n_col(m-1)} > 57 & str{n_col(m-1)} < 127)))
			n_column = max(n_c1,n_c2);
		else
			n_column = n_c1;
		end
	else
		n_column = n_col;
	end

	% Well, this is a stupid patch for the case where the file has only one column.
	% In that case the above test failed
	if (isempty(n_column) || n_column == 0),     n_column = 1;   end     % (unless the file is empty!!!)

	% If no header number request was made, return right away
	if (nargout <= 3),		return,		end

	% Now test if header lines are present (ascii 65:122 contain upper and lower case letters)
	for (j = 1:n_headers_candidate)
		head = find((str{j} > 32 & str{j} < 43) | str{j} > 58);
		tmp = find(str{j} == 78);    % I'm searching for a NaN string (ascii 78,97,78)
		if ~isempty(tmp) && length(tmp) >= 2
			tmp = find(str{j}(tmp(1)+1) == 97);
		else
			tmp = [];
		end
		if ~isempty(tmp)
			% NaN string found. But now we have a problem. If we reach here it's because a text string
			% was found. But that may equaly happen on a header line, or a data line with NaNs.
			% So how to decide which was the case? The way that occurs to me is to decide on basis
			% of the number of columns, but this is far from being 100% full proof.
			if n_col(j) == n_column
				head = [];              % NaNs on a data line, or in a header line with same n columns as data line
			else
				head = 1;               % NaNs on a header line
			end
			% However, exclude eventualy NaNs in a header line with the same n columns as a data line
			if isempty(head);
				w = strrep(str{j},'NaN','');
				head = find((w > 31 & w < 43) | w > 58);
			end
		end
		if ~isempty(head),	n_headers = n_headers + 1;	end
	end

	% Another trouble we may face is when there is a text column in file (e.g. a date string)
	% In that case the above test failed miserably and I don't really know what clever thing to do.
	% So we do a more dumb test. If first element of last line is a number we re-test again for
	% headers, but this time will accept only those lines whose first character is '#'
	if ((n_headers == m) && (m > 0))
		t = strtok(str{m});
		if (~isnan(t))
			n_headers = 0;
			j = 1;
			while (str{j}(1) == '#' && j <= m)
				n_headers = n_headers + 1;
				j = j + 1;
			end
		end
	end

	if (nargout == 5)
		isGSHHS = false;
		for (k = 1:n_headers)
			if (strfind(str{k}, 'GSHHS Master'))
				isGSHHS = true;		break
			end
		end
		if (~isGSHHS)		% Try if WDBII
			for (k = 1:n_headers)
				if (strfind(str{k}, 'WDBII Borders'))
					isGSHHS = true;		break
				end
			end
		end
	end

	% ----------------------------------------------------------------------------------------
	function guessed = guess_in_bin(fiche)
	% ...

	fid = fopen(fiche);
	out = fread(fid,24,'single');	out = reshape(out,2,12)';
	if ( (abs(out(1,1)) < 1e10) && (abs(out(1,2)) < 1e10) )
		rel = abs(std(out) ./ out(1,:));
		if (all(rel < 0.1 ))
			guessed.nCols = 2;		guessed.type = 'single';
			fclose(fid);			return
		end
	elseif ( isnan(out(1,1)) && isnan(out(1,2)) && ~isnan(out(2,1)) )
		guessed.nCols = 2;			guessed.type = 'single';
		fclose(fid);				return
	end

	frewind(fid);
	out = fread(fid,24,'double');	out = reshape(out,2,12)';
	if ( (abs(out(1,1)) < 1e10) && (abs(out(1,2)) < 1e10) )
		rel = abs(std(out) ./ out(1,:));
		if (all(rel < 0.1 ))
			guessed.nCols = 2;		guessed.type = 'double';
			fclose(fid);			return
		end
	elseif ( isnan(out(1,1)) && isnan(out(1,2)) && ~isnan(out(2,1)) )
		guessed.nCols = 2;			guessed.type = 'double';
		fclose(fid);				return
	end

	frewind(fid);
	out = fread(fid,36,'single');	out = reshape(out,3,12)';
	if ( (abs(out(1,1)) < 1e10) && (abs(out(1,2)) < 1e10) && (abs(out(1,3)) < 1e10) )
		rel = abs(std(out) ./ out(1,:));
		if (all(rel(1:2) < 0.1 & rel(3) < 0.5))	% Relax condition on 3rth column because data can be more disperse
			guessed.nCols = 3;		guessed.type = 'single';
			fclose(fid);			return
		end
	elseif ( isnan(out(1,1)) && isnan(out(1,2)) && isnan(out(1,3)) )
		guessed.nCols = 3;			guessed.type = 'single';
		fclose(fid);				return
	end

	frewind(fid);
	out = fread(fid,36,'double');	out = reshape(out,3,12)';
	if ( (abs(out(1,1)) < 1e10) && (abs(out(1,2)) < 1e10) && (abs(out(1,3)) < 1e10) )
		rel = abs(std(out) ./ out(1,:));
		if (all(rel(1:2) < 0.1 & rel(3) < 0.5))	% Relax condition on 3rth column
			guessed.nCols = 3;		guessed.type = 'double';
			fclose(fid);			return
		end
	elseif ( isnan(out(1,1)) && isnan(out(1,2)) && isnan(out(1,3)) )
		guessed.nCols = 3;			guessed.type = 'double';
		fclose(fid);				return
	end

	guessed = true;				% TRUE to mean that it's binary anyway though we didn't guess its organization
	fclose(fid);

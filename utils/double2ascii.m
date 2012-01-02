function double2ascii(filename, X, formatStr, multiseg)
% DOUBLE2ASCII - Writes double array X to an output ASCII file.
%
% If a single format specifier is specified in the input format
% string, that format will be used for all columns of X.
% The user may also specify different formats for each column of the double array X.
%
% Inputs:  filename  - Name of the output ASCII file
%		X  - double array can be a vector (1-D), a matrix (2-D) or a cell array.
%		In later case each cell must contain a Mx2 array and output file will be multisegment
%		formatStr  - OPTIONAL format string  (Default = '%f')
%					Provide one string with formats if you want to use different ones
%		multiseg   - OPTIONAL
%					If X is a cell array and MULTISEG is a char string (whatever)
%					write a multisegment file separated with the GMT '>' multisegment flag.
%					The same result is acomplished if X is an array with NaNs as separators.
%					If X and MULTISEG are both cells with the same number of elements,
%					use MULTISEG{k} contents as multiseg flag separator. 
%
% Example 1: Export array X to ASCII file with the same format for all columns
%		X = rand(300,10);
%		double2ascii('foo1.txt', X, '%4.2f ')
%
% Example 2: Export array X to ASCII file with different formats for each column
%		year = (1991:2000)';
%		x = (1:10)';
%		column2 = x / 100;         column3 = x * 1e27;
%		X = [ year column2 column3];
%		double2ascii('foo2.txt', X, '%d  %5.2f  %10.3e');

% Distantly rooted on file with the same name by Denis Gilbert
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

	% Set default format string
	n_arg = nargin;
	if (n_arg < 3),		formatStr = '%f';	end
	if (~(isnumeric(X) || isa(X,'cell')))
		error('Input variable ''X'' must be numeric or cell array')
	end
	if (n_arg < 2)
		error('At least two input arguments are required: ''filename'' and ''X''')
	end
	if (~isa(X,'cell') && ndims(X) > 2)
		error('Input variable ''X'' cannot have more than 2 dimensions')	% Cell case is not tested
	end
	if (isa(X,'cell'))				% May or not be printed as multisegment files
		do_multiseg = true;			% Force multisegment, unless unset
		if (n_arg <= 3)
			n_arg = 4;		do_multiseg = false;	multiseg = [];
		elseif (~isa(multiseg,'cell'))
			multiseg = [];
		end
	end

	%Find all occurrences of the percent (%) format specifier within the input format string
	kpercent = numel(strfind(formatStr,'%'));		% IF 'kpercent' == 1, ==> Same format for ALL columns

	%Open and write to ASCII file
	if (ispc),			fid = fopen(filename,'wt');
	elseif (isunix)		fid = fopen(filename,'w');
	else				error('DOUBLE2ASCII: Unknown platform.');
	end

	if (~isa(X,'cell'))
		ncols  = size(X,2);					% Determine the number of rows and columns in array X
		fmt = make_print_format(formatStr, ncols, kpercent);
	end

	if (n_arg < 4)							% Original form. No eventual NaN cleaning
		fprintf(fid, fmt, X');
	else									% We might have NaNs (that is multi-segments files)
		if (isa(X,'cell'))
			for (k = 1:length(X))			% Currently deals only with Mx2 arrays case
				if (do_multiseg && isempty(multiseg))
					fprintf(fid,'%s\n','>');
				elseif (do_multiseg)
					fprintf(fid,'%s\n', multiseg{k});	% Write out the multisegment info we got ininput					
				end
				fmt = make_print_format(formatStr, size(X{k},2), kpercent);
				fprintf(fid, fmt, X{k}');
			end
		elseif ( ~any(isnan(X)) )		% NO, we haven't
			fprintf(fid, fmt, X');
		else							% YES, we have them (then multisegs)
			if (ncols == 2)
				[y_cell,x_cell] = aux_funs('localPolysplit',X(:,2),X(:,1));
				for (k = 1:numel(x_cell))
					fprintf(fid,'%s\n','>');
					fprintf(fid, fmt, [x_cell{k}(:)'; y_cell{k}(:)']);
				end
			elseif (ncols == 3)
				[y_cell,x_cell,z_cell] = aux_funs('localPolysplit',X(:,2), X(:,1), X(:,3));
				for (k = 1:numel(x_cell))
					fprintf(fid,'%s\n','>');
					fprintf(fid, fmt, [x_cell{k}(:)'; y_cell{k}(:)'; z_cell{k}(:)']);
				end
			else
				warndlg('Saving multiseg by NaN works only with 2 or 3 columns. Nothing saved.','warning')
			end
		end
	end

	fclose(fid);

%---------------------------------------------------------------------------------
function fmt = make_print_format(formatStr, ncols, kpercent)
% Create the format string for fprintf by repetition of one base format type
	if (kpercent == 1)
		if ( strcmp(formatStr(end-1:end),'\n') ),	formatStr = formatStr(1:end-2);		end
		if ( ~strcmp(formatStr(end-1:end),'\t') )
			formatStr = [formatStr '\t'];
		end
		fmt = repmat(formatStr, [1,ncols-1]);
		fmt = [fmt formatStr(1:end-2) '\n'];
	else
		if ( strcmp(formatStr(end-1:end),'\n') )
			fmt = formatStr;
		else
			fmt = [formatStr '\n'];
		end
	end

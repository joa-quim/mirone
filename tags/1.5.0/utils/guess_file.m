function [bin,n_column,multi_seg,n_headers] = guess_file(fiche, opt1, opt2)
% [bin,n_column,multi_seg,n_headers] = guess_file(fiche, opt1, opt2) trys to guess if file "fiche"
% is ascii or binary. 
% If it detects that "fiche" is ascii this function tries to find out wether the multisegment
% symbol (">") is present, the number of columns in the file and if it has header lines.
% NOTE that this last tests may not always give reliable results.
% OPT1, if given, will be MAXCHARS
% OPT2, if given, will be nl_max

% Error testing
bin = 0;    multi_seg = 0;  n_headers = 0;  n_column = 0;
	n_args = nargin;
	if (~n_args)
		errordlg('function guess_file: must give an input file name','File Error')
		return
	elseif (n_args == 1)
		MAXCHARS = 1024;        % Maximum characters to load from file
		nl_max = 30;            % Maximum number of file lines to use in tests
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
	if (any(A > 126 & A < 192))    % Binary files have bytes with values greater than 126 (but so is the ç char)
        bin = 1;
        return
	end
	clear A;

	str = strread(str,'%s','delimiter','\n');
	if (count == MAXCHARS)      % In these cases last line is normally incomplete
        str(end) = [];
	end
	
	nl_max = min(nl_max,numel(str));
	if (nl_max == 0),	bin = [];	return,		end

    % Make a crude test to find number of columns and number of headers
    n_col(1:nl_max) = 0;    n_multi = 0;
    idM = false(1,nl_max);
    for (i =1:nl_max)
        if isempty(str{i});  continue,   end		% Jump blank lines
        if (str{i}(1) == '#' || str{i}(1) == '%'),		continue,   end		% Jump known comment lines
        if strfind(str{i}(1:min(2,length(str{i}))),'>')   % multisegmet line, so the rest of it is of no interest (but count them)
            n_multi = n_multi + 1;
            idM(i) = true;      % Tag it to deletion
            continue
        end
        str{i} = deblank(str{i});          % Blanks make a mess to the guessing code
        [tok,rem]=strtok(str{i});
        if ~isempty(rem);   n_col(i) = n_col(i) + 1;    end     % count first column
        while ~isempty(rem)
            [tok,rem]=strtok(rem);
            n_col(i) = n_col(i) + 1;
        end
    end
    multi_seg = n_multi;
    str(idM) = [];
    n_col(idM) = [];

    % Now decide how many columns have the data lines. The easeast is to assume that the info is in the last line
    % However this may fail if last line contains, for example, the multisegment symbol (">").
    % So, do another test.
	n_col(~n_col) = [];				% Remove zeros which correspond to comment lines (the ones starting by '#') 
    m = min(nl_max,numel(n_col));

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
    for j=1:m
        head = find((str{j} > 32 & str{j} < 43) | str{j} > 58);
        tmp = find(str{j} == 78);    % I'm searching for a NaN string (ascii 78,97,78)
        if ~isempty(tmp) && length(tmp) == 2
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
        if ~isempty(head);   n_headers = n_headers + 1;  end
    end

function sout = ddewhite(s,opt)
%DDEWHITE Double dewhite. Strip both leading and trailing whitespace.
%
%   DDEWHITE(S) removes leading and trailing white space and any null
%   characters from the string S.  A null character is one that has an
%   absolute value of 0.
%   OPT = '0' will remove trailing zeros from the string S (and not the blanks)
%
%   See also DEWHITE, DEBLANK, DDEBLANK.

%   Author:      Peter J. Acklam
%   Time-stamp:  2002-03-03 13:45:06 +0100
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

%   Added an option to remove trailing zeros from a string
%   J. Luis

%error(nargchk(1, 1, nargin));
if ~ischar(s)
    error('Input must be a string (char array).');
end

if (isempty(s)),    sout = s;   return;      end
if (nargin == 1),   opt = [];   end

if (isempty(opt))
	[r, c] = find(~isspace(s));
	if (size(s, 1) == 1)
        sout = s(min(c):max(c));
	else
        sout = s(:,min(c):max(c));
	end
else
	if (size(s, 1) == 1)
        ii = length(s);     jj = ii;  n = 0;
        while (strcmp(s(ii),'0'))
            ii = ii - 1;    n = n + 1;
        end
        if (n > 0),     sout = s(1:jj-n);
        else            sout = s;       end
    else
        error('Function ddewhite is not programed to deal with multiline strings.');
    end    
end
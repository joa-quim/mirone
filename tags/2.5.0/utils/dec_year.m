function yd = dec_year(varargin)
%   DEC_YEAR(YEAR, MONTH, DAY, HOUR, MINUTE, SECOND) returns the ordinal year
%   number plus a fractional part depending on the month, day, and time of day
%
%   Any missing MONTH or DAY will be replaced by 1.  HOUR, MINUTE or SECOND
%   will be replaced by zeros.
%
%   If no date is specified, the current date and time is used.

%   Adapted from timeutil functions of Peter J. Acklam by Joaquim Luis

nargsin = nargin;
if nargsin
    argv = { 1 1 1 0 0 0 };
    argv(1:nargsin) = varargin;
else
    argv = num2cell(clock);
end
[year, month, day, hour, minute, second] = deal(argv{:});

days_in_prev_months = [0 31 59 90 120 151 181 212 243 273 304 334]';

% Day in given month.
try
    yd = days_in_prev_months(month) ...               % days in prev. months
         + ( isleapyear(year) & ( month > 2 ) ) ...   % leap day
         + day ...                                    % day in month
         + ( second + 60*minute + 3600*hour )/86400;  % part of day
catch
    yd = [];    return
end

yd = year + (yd - 1) ./ (365 + isleapyear(year));

%--------------------------------------------------------------------------
function t = isleapyear(year)
t = ( ~rem(year, 4) & rem(year, 100) ) | ~rem(year, 400);

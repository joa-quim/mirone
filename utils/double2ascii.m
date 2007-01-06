function double2ascii(filename, X, formatString, multiseg)
% DOUBLE2ASCII - Writes double array X to an output ASCII file.
%
% If a single format specifier is specified in the input format
% string, that format will be used for all columns of X.
% The user may also specify different formats for each column of the double array X.
%
% Syntax:  double2ascii(filename, X, formatString, multiseg)
%
% Inputs:  filename  - Name of the output ASCII file
%          X  - double array can be a vector (1-D) or matrix (2-D)
%          formatString  - OPTIONAL format string  (Default = '%f')
%          multiseg      - OPTIONAL Replace NaNs lines with the GMT '>' multisegment flag
%
% Example 1: Export array X to ASCII file with the same format for all columns
%            X = rand(300,10);
%            double2ascii('foo1.txt', X, '%4.2f ')
%
% Example 2: Export array X to ASCII file with different formats for each column
%            year = (1991:2000)';
%            x = (1:10)';
%            column2 = x / 100;         column3 = x * 1e27;
%            X = [ year column2 column3];
%            double2ascii('foo2.txt', X, '%d  %5.2f  %10.3e');

% Author: Denis Gilbert, Ph.D., physical oceanography
% Maurice Lamontagne Institute, Dept. of Fisheries and Oceans Canada
% email: gilbertd@dfo-mpo.gc.ca  Web: http://www.qc.dfo-mpo.gc.ca/iml/
% September 2001; Revision: 25-Apr-2002
%
% J. Luis   06-01-2007  Added multisegment writing option

% Set default format string
if (nargin < 3),        formatString = '%f ';  end
if (~isnumeric(X)),     error('Input variable ''X'' must be numeric');   end
if (nargin < 2),        error('At least two input arguments are required: ''filename'' and ''X'''); end
if (ndims(X) > 2),      error('Input variable ''X'' cannot have more than 2 dimensions');     end

%Find all occurrences of the percent (%) format specifier within the input format string
kpercent = strfind(formatString,'%');

%Open and write to ASCII file
if (ispc),      fid = fopen(filename,'wt');
elseif (isunix) fid = fopen(filename,'w');
else            fclose(fid);    error('DOUBLE2ASCII: Unknown platform.');
end
ncols  = size(X,2);     % Determine the number of rows and columns in array X

if (nargin < 4)                 % Original form. No eventual NaN cleaning
	if (kpercent == 1)          % Same format for ALL columns
        fprintf(fid,[repmat(formatString,[1,ncols]) '\n'], X');
	else                        % Different format for each column
       fprintf(fid,[formatString '\n'], X');
	end
else                            % We might have NaNs (that is multi-segments files)
    if ( ~any(isnan(X)) )       % NO, we haven't
		if (kpercent == 1),     fprintf(fid,[repmat(formatString,[1,ncols]) '\n'], X');
		else                    fprintf(fid,[formatString '\n'], X');
		end
    else                        % YES, we have them (then multisegs)
        [y_cell,x_cell] = localPolysplit(X(:,2),X(:,1));
        for (k=1:numel(x_cell))
            fprintf(fid,'%s\n','>');
    		if (kpercent == 1),     fprintf(fid,[repmat(formatString,[1,ncols]) '\n'], [x_cell{k}(:)'; y_cell{k}(:)']);
		    else                    fprintf(fid,[formatString '\n'], [x_cell{k}(:)'; y_cell{k}(:)']);
		    end
            %fprintf(fid,'%.5f\t%.5f\n',[x_cell{k}(:)'; y_cell{k}(:)']);
        end
    end
end

fclose(fid);

% --------------------------------------------------------------------------------
function [latcells,loncells] = localPolysplit(lat,lon)
%POLYSPLIT Extract segments of NaN-delimited polygon vectors to cell arrays
%
%   [LATCELLS,LONCELLS] = POLYSPLIT(LAT,LON) returns the NaN-delimited
%   segments of the vectors LAT and LON as N-by-1 cell arrays with one
%   polygon segment per cell.  LAT and LON must be the same size and have
%   identically-placed NaNs.  The polygon segments are column vectors if
%   LAT and LON are column vectors, and row vectors otherwise.

% Copyright 1996-2006 The MathWorks, Inc.
% $Revision: 1.4.4.5 $    $Date: 2006/05/24 03:35:26 $

	[lat, lon] = localRemoveExtraNanSeps(lat, lon);
	indx = find(isnan(lat(:)));         % Find NaN locations.
	
	% Simulate the trailing NaN if it's missing.
	if ~isempty(lat) && ~isnan(lat(end))
        indx(end+1,1) = numel(lat) + 1;
	end
	
	%  Extract each segment into pre-allocated N-by-1 cell arrays, where N is
	%  the number of polygon segments.  (Add a leading zero to the indx array
	%  to make indexing work for the first segment.)
	N = numel(indx);
	latcells = cell(N,1);       loncells = cell(N,1);
	indx = [0; indx];
	for k = 1:N
        iStart = indx(k)   + 1;
        iEnd   = indx(k+1) - 1;
        latcells{k} = lat(iStart:iEnd);
        loncells{k} = lon(iStart:iEnd);
	end

% --------------------------------------------------------------------------------
function [xdata, ydata, zdata] = localRemoveExtraNanSeps(xdata, ydata, zdata)
    %removeExtraNanSeps  Clean up NaN separators in polygons and lines

	p = find(isnan(xdata(:)'));     % Determing the positions of each NaN.
	
	% Determine the position of each NaN that is not the final element in a sequence of contiguous NaNs.
	q = p(diff(p) == 1);
	
	% If there's a leading sequence of NaNs (a sequence starting with a NaN in
	% position 1), determine the position of each NaN in this sequence.
	if isempty(p),      r = [];
	else                r = find((p - (1:numel(p))) == 0);
	end
	
	% Determine the position of each excess NaN.
	if isempty(r),      s = q;
	else                s = [r q(q > r(end))];
	end
	
	% Remove the excess NaNs.
	xdata(s) = [];      ydata(s) = [];
	if (nargin >= 3),   zdata(s) = [];  end

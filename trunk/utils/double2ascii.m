function double2ascii(filename, X, formatString)
%DOUBLE2ASCII - Writes double array X to an output ASCII file.
% If a single format specifier is specified in the input format
% string, that format will be used for all columns of X.
%
% The user may also specify different formats for 
% each column of the double array X.
%
% Syntax:  double2ascii(filename, X, formatString);
%
% Inputs:  filename  - Name of the output ASCII file
%          X  - double array can be a vector (1-D) or matrix (2-D)
%          formatString  - OPTIONAL format string  (Default = '%f')
%
% Output:  none
%
% Example 1: Export array X to ASCII file with the same format for all columns
%
%            X = rand(300,10);
%            double2ascii('foo1.txt', X, '%4.2f ')
%
% Example 2: Export array X to ASCII file with different formats for each column
%
%            year = (1991:2000)';
%            x = (1:10)';
%            column2 = x / 100;
%            column3 = x * 1e27;
%            X = [ year column2 column3];
%            double2ascii('foo2.txt', X, '%d  %5.2f  %10.3e');
% 
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: SAVE  (several  -ascii options are available)

% Author: Denis Gilbert, Ph.D., physical oceanography
% Maurice Lamontagne Institute, Dept. of Fisheries and Oceans Canada
% email: gilbertd@dfo-mpo.gc.ca  Web: http://www.qc.dfo-mpo.gc.ca/iml/
% September 2001; Last revision: 25-Apr-2002

% Set default format string
if nargin < 3 , formatString = '%f ';  end

if nargin < 2
    error('At least two input arguments are required:  ''filename'' and ''X''')
end

if ~isnumeric(X)
   error('Input variable ''X'' must be numeric');
end

if ndims(X) > 2
   disp('Double array ''X'' must be a vector or a matrix (2-D array)')
   error('Input variable ''X'' cannot have more than 2 dimensions')
end

%Determine the number of rows and columns in array X
ncols  = size(X,2);

%Find all occurrences of the percent (%) format specifier
%within the input format string
kpercent = strfind(formatString,'%');

%Open and write to ASCII file
if ispc;        fid = fopen(filename,'wt');
elseif isunix;  fid = fopen(filename,'w');
else    errordlg('Unknown platform.','Error');
end
if kpercent == 1
   %Same format for ALL columns
   fprintf(fid,[repmat(formatString,[1,ncols]) '\n'], X');
else
   %Different format for each column
   fprintf(fid,[formatString '\n'], X');
end

fclose(fid);


function matstr = mat2clip(matrix,n)
% MAT2CLIP Formats a 2-D matrix to clipboard formated text array
%
%   STR = MAT2CLIP(MAT) converts the 2-D matrix MAT to a MATLAB 
%   string so that the Windows clipboard can paste it as a matrix
%	in any text editor, or better, in EXCEL.
%	MAT must be a numeric array (doubles or logicals).
%
%	Default reproduces the original matrix to within 15 digits of precision.
%   STR = MAT2STR(MAT,N) uses N digits of precision. 
%   MAT2STR(MAT) without output puts the result directly in the clipboard
%				 using the CLIPBD_MEX mex file.
%
%   Example
%       mat2clip(rand(10,2)) produces the string simulating a 10x2 matrix
%							that can be directly pasted in excel.
%
%	QUALTY-CONTROL: THIS PROGRAM IS JAVA-FREE
%
%	Joaquim Luis (jluis@ualg.pt), Sept 2008
%	Based in MAT2STR

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

	if (nargout && ~strncmp(computer,'PC',2))
		error('mat2clip:badusage','Without output, this function will works on Windows');
	end
	if (nargin < 1 || nargin > 2)
		error('mat2clip:badinput','Wrong number of arguments.');
	end
	if ischar(matrix),
		matstr = matrix;		return
	end
	if isempty(matrix)
		matstr = '';	return
	end
	if (~isnumeric(matrix))
		error('mat2clip:badinput','Input mmatrix must be of numeric type.');
	end

	if (ndims(matrix) > 2)
		error('mat2clip:badinput','Input matrix must be 2-D.');
	end

	if ( ~isa(matrix,'double') || ~islogical(matrix) )
		% Detect which matlab version is in use. R13 doesn't accept other than doubles on sprintf - beautiful stupidity
		ver = version;
		if (double(ver(1)) < 55)		% Older than 7
			matrix = double(matrix);
		end
	end

	[nrows, ncols] = size(matrix);

	if (nargin < 2),	n = 15;		end
	form = sprintf('%%.%df%%c',n);

	% now make an estimate on how big string will need to be 
	% covers (tab) between columns, the decimal point and the unknown number of integer chars
	n_int_chars = numel(sprintf('%d',fix(matrix(1,:))));	% Estimate number of integer characters
	spaceRequired = nrows * ((n+n_int_chars+2) * ncols + 2);

	string = '';
	string(1,spaceRequired) = char(0);

	CR_LF = [char(13) char(10)];	% Windows clipboard wants this to separate lines
	TAB   = char(9);

	pos = 1;
	for i = 1:nrows
		for j = 1:ncols
			if (matrix(i,j) == Inf)
				string(pos:pos+2) = 'Inf';
				pos = pos + 3;
			elseif (matrix(i,j) == -Inf)
				string(pos:pos+3) = '-Inf';
				pos = pos + 4;
			elseif islogical(matrix(i,j))
				if matrix(i,j)		% == true
					string(pos:pos+3) = 'true';
					pos = pos + 4;
				else
					string(pos:pos+4) = 'false';
					pos = pos + 5;
				end
			else
				tempStr = sprintf(form,matrix(i,j),TAB);
				len = numel(tempStr)-1;
				string(pos:pos+len) = tempStr;
				pos = pos+len+1;
			end
		end
		string(pos-1:pos) = CR_LF;
		pos = pos + 1;
	end
	
	% clean up the end of the string
	string = string(1:pos-3);		% 3 = 1 (extra pos + 1) + 2 (= length CR_LF)

	% Without output, send it directly to clipboard
	if (nargout == 0),		clipbd_mex(string)
	else					matstr = string;
	end

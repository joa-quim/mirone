function val = string_token(s,st)
% STRING_TOKEN(S,ST) Decomposes the string S in tokens delimited by ST
% STRING_TOKEN(S) Decomposes the string S in tokens delimited by blanks
% val = STRING_TOKEN(S,ST) (or STRING_TOKEN(S)) returns a cell matrix with the individual
% tokens. If an error ocurrs or with an empty imput, it returns an empty cell matrix

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

if isempty(s)
    val = [];
    return
end
s = deblank(s);
if nargin == 1         % Called like in the MATLAB strtok with delimeter == blanks
    error = 0;
    [t,r]=strtok(s);
    unit{1} = t;
    if isnan(str2double(unit{1}))
        str = [unit{1} '   is not a valid number'];
        error = error + 1;
        errordlg(str,'Error in string_token')
    end
    i = 2;
    while ~isempty(r)
        [t,r]=strtok(r);
        unit{i} = t;
        if isnan(str2double(unit{i}))
            str = [unit{i} '   is not a valid number'];
            error = error + 1;
            errordlg(str,'Error in string_token')
        end
        i = i + 1;
    end
    if error
        val = [];
    else
        val = unit;
    end
elseif nargin == 2
    error = 0;
    [t,r]=strtok(s,st);
    unit{1} = t;
    if isnan(str2double(unit{1}))
        str = [unit{1} '   is not a valid number'];
        error = error + 1;
        errordlg(str,'Error')
    end
    i = 2;
    while ~isempty(r)
        [t,r]=strtok(r,st);
        unit{i} = t;
        if isnan(str2double(unit{i}))
            str = [unit{i} '   is not a valid number'];
            error = error + 1;
            errordlg(str,'Error')
        end
        i = i + 1;
    end
    if error
        val = [];
    else
        val = unit;
    end
else
    errordlg('Error in string_token, called with a wrong number of arguments','Error')
    val = [];
    return
end

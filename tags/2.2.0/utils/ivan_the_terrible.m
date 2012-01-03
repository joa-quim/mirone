function inc = ivan_the_terrible(dx,n,mod,one_or_zero)
% This function was writen to deal with the TERRORIST dictatorship of GMT_grd_RI_verify.
% INC = IVAN_THE_TERRIBLE(DX,N,MODE)
% If MODE == 1, where DX is a real and N an integer (number of lines), returns a text string with INC
% that complies with GMT_grd_RI_verify.
% If MODE == 2, N now contains the candidate INC. The new INC is computed by a very simple algorithm
% (see below) that, on the contrary of a previous one (supposadly clever but very slow), works prety well.
% ONE_OR_ZERO is a scalar with either 1 (grid registration) or 0 (pixel registration)

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

if (nargin < 3)
    errordlg('function ivan_the_terrible called with a wrong number of arguments','Error')
    return
end
if  (nargin == 3)           % If not given, defaults to grid registration
    one_or_zero = 1;
end

if mod == 1
    inc = dx / (n - one_or_zero);
elseif mod == 2
    t_inc = n;
    if (dx / t_inc) < 2
        msgbox(['You must be joking. With this increment your grid would have only two lines ' ...
            'along this dimension. I''ll give you the least idiot increment, which is the one that ' ...
            'alows three lines. But you better reconsider your choice.'],'Chico Clever');
        inc = dx / 2;
        return
    end
    nl = round(dx / t_inc) + one_or_zero;
    inc = dx / (nl - one_or_zero);
else
    errordlg('function ivan_the_terrible called in a wrong mode','Error')
end

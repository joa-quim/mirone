function [chunk_s,chunk] = tile(m,height,skirt)
% TILE chunks the 1:m vector into a nx2 matrix. Each line contains the indexes
% of successif [(height*n_chunk - 2) (height*n_chunk + 2)]. SKIRT is the overlapping 
% zone between adjacent chunks.
% This function is used to split a mxn matrix in horizontal tiles of height 'height', but
% with an overlaping zone. The reason I wrote this is to deal with the Matlab voracity
% for RAM. So instead of calling extremelly memory consumption routines (like surfnorm
% used in illumination) with the full matrix, I can call it several times (one for
% each tile).
% CHUNK_S contains the indexes of the successive chunks including the overlaping zone
% CHUNK   contains the indexes of CHUNK_S that do not overlap when put side-by-side
% (or vertically) to reconstruc the hole image.
% Example:      [xx yy] = tile(402,100,2);
% produces the following result
% xx =   1     102      yy =    1   100
%        99    202              3   102  
%        199   302              3   102
%        299   402              3   104

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

if (nargin == 1),       height = 200;     skirt = 2;
elseif (nargin == 2),   skirt = 2;
end

chunk_s(1,1) = 1;                                  % First index is allways one
tmp = 1:height:m;  tmp(2:end) = tmp(2:end) - 1;
if length(tmp) > 1
    chunk_s(1:length(tmp)) = tmp + skirt;
    yy = chunk_s(2:end);
    xx = yy - skirt*2 + 1;
    chunk_s = [[1 xx(1:end-1)]' yy'];
else
    chunk_s(1,2) = m;
end
nl = size(chunk_s,1);
if chunk_s(nl,2) > m                               % If last index is supperior to m, set it to m
    chunk_s(nl,2) = m;
elseif chunk_s(nl,2) == m                          % Do nothing. We have the result
else                                               % Else build a last line with the remainder.
    chunk_s(nl+1,1) = chunk_s(nl,2) - skirt*2 + 1;
    chunk_s(nl+1,2) = m;
end

if size(chunk_s,1) > 1                              % Build a second nx2 matrix with the indexes of
    chunk = chunk_s;                                % Just make a copy of chunk_s to initialize chunk
    chunk(1,1) = 1;                chunk(1,2) = height;
    chunk(2:end,1) = skirt + 1;    chunk(2:end,2) = height + skirt;
    chunk(end,2) = m - (size(chunk_s,1)-1)*height + skirt;
else
    chunk(1,1) = 1;     chunk(1,2) = m;    
end

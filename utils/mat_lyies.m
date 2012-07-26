function [s,w] = mat_lyies(str,opt)
% [s,w] = mat_lyies(str) is a replacement to the "dos" or "unix" command when
% used in codes that will be compiled. The Matlab original doesn't accept two
% outputs when they are compiled. 
% OPT is an optional text string with the file name that will hold the output
% to stdout. This file will not be deleted upon exit.

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

if nargin == 2
    text_stdout = opt;
else
    text_stdout = 'stdout_txt';
end

w = '';   s = 0;
run_comm = [str ' > ' text_stdout ' 2> text_stderr'];

if (str(end) == '&')
	str(end) = [];
	run_comm = [str ' > ' text_stdout ' 2> text_stderr &'];
end

if isunix           % UNIX
    unix(run_comm);
elseif ispc         % Windows
    dos(run_comm);
else
    errordlg('Unknown platform.','Error');  return;
end

% Now find out if it worked and get the output message into w

fid = fopen('text_stderr');
while 1
    tline = fgetl(fid);
    if ~ischar(tline),	break,	end
    if strfind(tline,'is not recognized')			% case windows error
        s = 1;        w = [str '. Command not found'];
    elseif strfind(tline,'command not found')		% case unix (or cygwin) error
        s = 1;        w = [str '. Command not found'];
    elseif strfind(tline,'ERROR')
        s = 1;        w = tline;
    elseif ~isempty(tline)          % Gave up. I cannot guess all error messages, so if any report it
        s = 1;        w = tline;
    end
end
fclose(fid);

if (s == 1),	return,		end

fid = fopen(text_stdout);
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    w = [w ' ' tline];
end

fclose(fid);
if (nargin ~= 2),    delete stdout_txt;     end
delete text_stderr;

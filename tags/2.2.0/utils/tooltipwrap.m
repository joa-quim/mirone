function tooltipwrap(hObject,texto,nCols)
%TOOLTIPWRAP Uses the information in TEXTO as the 'TooltipString' for the HOBJECT UI Control.
%
% HOBJECT, handle of the object to which the tooltip string will be applyied.
% TEXTO can be a linestring, a cell array of strings or the name of a file containing
% the desired text info.
% NCOLS, is the width in chars of the wrapped TEXTO matrix. If absent it defaults to 20.
% When TEXTO is a filename and NCOLS = Inf, the file contents will be used with no change
% (that is, no text wrapping)
%
%   Example:
%           h=figure;
%           hui=uicontrol('Parent',h,'Units','normalized','Position',[.05 .08 .9 .9],'Style','edit');
%           txt= ['Xi, so much room here to write my thoughts. And this ',...
%               'little boring tooltip window that doesn''t leave me alone']
%           tooltipwrap(hui,txt,15)

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

% Error check
msg = [];
if (nargin < 2 || nargin > 3)
    msg = 'Wrong number of arguments';
elseif (~ishandle(hObject) || ~strcmp(get(hObject,'Type'),'uicontrol') )
    msg = 'First argument is not a valid uicontrol handle.';
elseif ( ~(ischar(texto) || iscell(texto)) )
    msg = 'Second argument must be a line string or a cell array of chars.';
end

if (~isempty(msg))
    h = errordlg(['Tooltipwrap: ' msg],'Error');
    hButt = findobj(h,'Style','pushbutton');
    tooltipwrap(hButt,{msg},10)
    return
end
if (nargin ~= 3),   nCols = 20;     end
    
% Try to open 'texto' as a file. If we can, use the file contents as 'TooltipString'
if (~iscell(texto))
    fid = fopen(texto);
    if (fid > 0)                % Yes, 'texto' is a file name
        str = fscanf(fid,'%c');  fclose(fid);
        if (~isinf(nCols))
            str = textwrap({str},nCols);
        end
    else                        % No, 'texto' is a text string
        str = textwrap({texto},nCols);
    end
else
    str = textwrap(texto,nCols);
end

if (iscell(str))
    str_tip = [];
    for (i = 1:length(str)-1)
        str_tip = sprintf([str_tip str{i} '\n']);
    end
    str_tip = sprintf([str_tip str{end}]);
else
    str_tip = str;
end
set(hObject,'TooltipString',str_tip)

function frame3D(parent,pos,color,bg_color,usr_dat)
% frame3D.  Define a frame with a 3D effect.

%	Copyright (c) 2004-2015 by J. Luis
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

% $Id$


% Build the rectangle's left vertical side
x1 = [pos(1) pos(2) 1 pos(4)];
uicontrol('Parent',parent, 'Style','frame', 'Position',x1, 'ForegroundColor',color,'UserData',usr_dat);
x2 = [pos(1)+1 pos(2)+1 1 pos(4)-2];
uicontrol('Parent',parent, 'Style','frame', 'Position',x2, 'ForegroundColor',[1 1 1],'UserData',usr_dat);

% Build the rectangle's right vertical side
x1 = [pos(1)+pos(3)-1 pos(2) 1 pos(4)];
uicontrol('Parent',parent, 'Style','frame', 'Position',x1, 'ForegroundColor',color,'UserData',usr_dat);
x2 = [pos(1)+pos(3) pos(2)-1 1 pos(4)+1];
uicontrol('Parent',parent, 'Style','frame', 'Position',x2, 'ForegroundColor',[1 1 1],'UserData',usr_dat);

% Build the rectangle's bottom side
x1 = [pos(1) pos(2) pos(3) 1];
uicontrol('Parent',parent, 'Style','frame', 'Position',x1, 'ForegroundColor',color,'UserData',usr_dat,'Tag','B');
x2 = [pos(1) pos(2)-1 pos(3)+1 1];
uicontrol('Parent',parent, 'Style','frame', 'Position',x2, 'ForegroundColor',[1 1 1],'UserData',usr_dat,'Tag','BB');

% Build the rectangle's top side
x1 = [pos(1) pos(2)+pos(4)-1 pos(3) 1];
uicontrol('Parent',parent, 'Style','frame', 'Position',x1, 'ForegroundColor',color,'UserData',usr_dat,'Tag','TT');
x2 = [pos(1)+1 pos(2)+pos(4)-2 pos(3)-2 1];
uicontrol('Parent',parent, 'Style','frame', 'Position',x2, 'ForegroundColor',[1 1 1],'UserData',usr_dat,'Tag','T');

%h1=uicontrol('Parent',parent, 'Style','frame', 'Position',position+[0 -1 1 1], 'ForegroundColor',[1 1 1]);
%h2=uicontrol('Parent',parent, 'Style','frame', 'Position',position, 'ForegroundColor',color);
%h3=uicontrol('Parent',parent, 'Style','frame', 'Position',position+[1 1 1-position(3) -2],'ForegroundColor',[1 1 1]);
%h4=uicontrol('Parent',parent, 'Style','frame', 'Position',position+[1 position(4)-2 -2 1-position(4)],...
%          'ForegroundColor',[1 1 1]);
%uistack(h4, 'bottom');  uistack(h3, 'bottom');  uistack(h2, 'bottom');  uistack(h1, 'bottom');

%if ~isempty(bg_color)
%    set(h2,'BackgroundColor',bg_color)    
%end

function out = new_frame3D(hFig, hText, hFrame)
% Give a 3D Pro Look to the miserable looking frame uis
%
% HTEXT eventualy contains handles to texts that need to be recreated
%	If HTEXT = [], fish all texts in Fig and recreate them above the 3D frames
%	If HTEXT = NaN, ignore texts and deal only with frames
%
%	new_frame3D(hFig, hText)	fish all frames in Fig
%
%	out = new_frame3D(...)	returns the new handles of the recreated text uis

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

	% Give a Pro look (3D) to the frame boxes 
	bgcolor = get(0,'DefaultUicontrolBackgroundColor');
	framecolor = max(min(0.65*bgcolor,[1 1 1]),[0 0 0]);
	if (nargin < 3)
		hFrame = findobj(hFig,'Style','Frame');
	end
	for (i = 1:numel(hFrame))
        frame_size = get(hFrame(i),'Position');
        f_bgc = get(hFrame(i),'BackgroundColor');
        usr_d = get(hFrame(i),'UserData');
		frame3D(hFig,frame_size,framecolor,f_bgc,usr_d)
		delete(hFrame(i))
	end

	% Recopy the text fields on top of previously created frames (uistack is too slow)
	if (isempty(hText))
		hText = findobj(hFig,'Style','Text');
	end
	if (isnan(hText)),		hText = [];		end		% So that next loop is not executed
	hText_new = zeros(1,numel(hText));
	for (i = 1:numel(hText))
        usr_d = get(hText(i),'UserData');
        t_size = get(hText(i),'Position');
		t_str = get(hText(i),'String');
		t_just = get(hText(i),'HorizontalAlignment');
		t_tag = get (hText(i),'Tag');
		fs = get(hText(i),'FontSize');
		fa = get(hText(i),'FontAngle');
		fw = get(hText(i),'FontWeight');
        fn = get(hText(i),'FontName');
		bgc = get (hText(i),'BackgroundColor');
		fgc = get (hText(i),'ForegroundColor');
        hText_new(i) = uicontrol('Parent',hFig, 'Style','text', 'Position',t_size,'String',t_str, ...
			'BackgroundColor',bgc, 'ForegroundColor',fgc, 'FontSize',fs, 'FontAngle',fa, 'FontWeight',fw, ...
			'FontName',fn, 'HorizontalAlignment',t_just, 'Tag',t_tag, 'UserData',usr_d);
	end
	delete(hText)
	if (nargout)	out = hText_new;	end

% ----------------------------------------------------------------------	
function frame3D(hFig,pos,color,bg_color,usr_dat)
%=======================================================================
% frame3D.  Define a frame with a 3D effect.
%=======================================================================
% 
% Build the rectangle's left vertical side
x1 = [pos(1) pos(2) 1 pos(4)];
uicontrol('Parent',hFig, 'Style','frame', 'Position',x1, 'ForegroundColor',color,'BackgroundColor',bg_color,'UserData',usr_dat);
x2 = [pos(1)+1 pos(2)+1 1 pos(4)-2];
uicontrol('Parent',hFig, 'Style','frame', 'Position',x2, 'ForegroundColor',[1 1 1],'BackgroundColor',bg_color,'UserData',usr_dat);

% Build the rectangle's right vertical side
x1 = [pos(1)+pos(3)-1 pos(2) 1 pos(4)];
uicontrol('Parent',hFig, 'Style','frame', 'Position',x1, 'ForegroundColor',color,'BackgroundColor',bg_color,'UserData',usr_dat);
x2 = [pos(1)+pos(3) pos(2)-1 1 pos(4)+1];
uicontrol('Parent',hFig, 'Style','frame', 'Position',x2, 'ForegroundColor',[1 1 1],'BackgroundColor',bg_color,'UserData',usr_dat);

% Build the rectangle's bottom side
x1 = [pos(1) pos(2) pos(3) 1];
uicontrol('Parent',hFig, 'Style','frame', 'Position',x1, 'ForegroundColor',color,'BackgroundColor',bg_color,'UserData',usr_dat,'Tag','B');
x2 = [pos(1) pos(2)-1 pos(3)+1 1];
uicontrol('Parent',hFig, 'Style','frame', 'Position',x2, 'ForegroundColor',[1 1 1],'BackgroundColor',bg_color,'UserData',usr_dat,'Tag','BB');

% Build the rectangle's top side
x1 = [pos(1) pos(2)+pos(4)-1 pos(3) 1];
uicontrol('Parent',hFig, 'Style','frame', 'Position',x1, 'ForegroundColor',color,'BackgroundColor',bg_color,'UserData',usr_dat,'Tag','TT');
x2 = [pos(1)+1 pos(2)+pos(4)-2 pos(3)-2 1];
uicontrol('Parent',hFig, 'Style','frame', 'Position',x2, 'ForegroundColor',[1 1 1],'BackgroundColor',bg_color,'UserData',usr_dat,'Tag','T');

%h1=uicontrol('Parent',hFig, 'Style','frame', 'Position',position+[0 -1 1 1], 'ForegroundColor',[1 1 1]);
%h2=uicontrol('Parent',hFig, 'Style','frame', 'Position',position, 'ForegroundColor',color);
%h3=uicontrol('Parent',hFig, 'Style','frame', 'Position',position+[1 1 1-position(3) -2],'ForegroundColor',[1 1 1]);
%h4=uicontrol('Parent',hFig, 'Style','frame', 'Position',position+[1 position(4)-2 -2 1-position(4)],...
%          'ForegroundColor',[1 1 1]);
%uistack(h4, 'bottom');  uistack(h3, 'bottom');  uistack(h2, 'bottom');  uistack(h1, 'bottom');

%if ~isempty(bg_color)
%    set(h2,'BackgroundColor',bg_color)    
%end
	

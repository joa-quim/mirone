function handFig = message_win(option, texto, varargin)
%MESSAGE_WIN  Message function.
%
%   Usage:  message_win('create','Some text')	For the first call
%			message_win('Some text')			For the first call (but VARARGIN mus be empty)
%           message_win('add','More text')      for subsequent calls. When the window message
%                                               is full of text a slider will be added.
%			message_win(...,varargin)			Where varargin contains PN/PV pairs. Valid PN are:
%					PN				Default value
% 				'figname',		'Message window'
% 				'width',		Guessed from text string
% 				'height',		Guessed from text string
% 				'bgcolor',		[.95 .95 .95]
% 				'fwcolor',		'k' (black)
% 				'fontweight',	'demi'
% 				'fontsize',		9
% 				'fontname',		'Helvetica'
% 				'position',		a position key as accepted by move2side (e.g, 'right')
%				'button'		'yes' creates a button that will plot the Grid's Min/Max
%								locations on the figure that called message_win.
%								But for that the necessary info must have been stored in this figure's
%								'UserData'. That is the responsability of the calling function.
%							Example:
%								ud.hAxe = handles.axes1;	ud.p1 = [xx_min yy_min];	
%								ud.p2 = [xx_max yy_max];	set(hMsg, 'UserData', ud)

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

tagFig = 'Wdmsgfig';
tagTxt = 'textTag';

if (nargin == 1)
	texto = option;		option = 'create';
end
hFig = findmsgwin(tagFig);
if ( isempty(hFig) && strcmp(option, 'add') )
	option = 'create';
end

if (~isa(texto,'cell')),		texto = cellstr(texto);	end

switch option
	case 'create'
		[figpos, figName, posTxt, movepos, bgcolor, fwcolor, Font, addButt, winMaxH, modal] = ...
			parse_inputs(texto, varargin{:});
		if (isempty(movepos))			% Reposition figure on screen
			figpos = getnicelocation(figpos, 'pixels');
		end
		hFig = figure('MenuBar','none', 'Name',figName, 'HandleVisibility', 'off', ...
				'Vis','off', 'Unit','pix', 'Pos',figpos,...
				'Color',bgcolor, 'NumberTitle','off', 'DoubleBuffer','on', 'Tag',tagFig);
		hTxt = uicontrol('Parent',hFig, 'Style','text', ...
				'Unit','pix', 'Pos', posTxt, ...
				'HorizontalAlignment','left',...
				'FontWeight',Font.Weight, 'FontSize',Font.Size, 'FontName',Font.Name,...
				'BackgroundColor',bgcolor, 'ForegroundColor',fwcolor, ...
				'Tag',tagTxt);
		posBut = [figpos(3)-44 4 40 21];
		hBut = uicontrol('Parent',hFig, 'Style','pushbutton', ...
				'Unit','pix', 'Position', posBut, ...
				'Str', 'OK', ...
				'Call', 'delete(gcbf)');
		if (addButt)
			uicontrol('Parent',hFig, 'Style','pushbutton', ...
				'Unit','pix', 'Position', [4 4 80 21], ...
				'Str', 'Plot min/max', ...
				'Call', @plotMiMa, ...
				'Tag','plotaMiMa');
		end

		if (modal)		% Make the Fig modal
			set(hFig,'WindowStyle','modal')
		end
		if (~isempty(movepos))			% Reposition figure on screen
			move2side(hFig, movepos)
		end
		set([hFig,hTxt],'units','norm');
		set(hTxt, 'String',texto)
		extent = get(hTxt,'Extent');
		if (isempty(extent)),	extent = zeros(1,4);	end		% Protection against empty texts

		% See if we should shrink the figure horizontally
		if ( extent(3) < 1 )					% Yes
			set([hFig,hTxt], 'unit', 'pix')
			extent_pix = get(hTxt,'Extent');
			figpos(3) = extent_pix(3)+20;
			figpos = getnicelocation(figpos, 'pixels');			% Reposition again
			set(hFig, 'Pos', figpos)
			posBut(1) = figpos(3)-44;
			set(hBut, 'Pos', posBut)			% We need to recenter the OK button too
			set([hFig,hTxt], 'unit', 'norm')
		end

		if ( extent(4) > 1 || addButt)		% Text too big to fit in?
			if (figpos(4) == winMaxH)		% Yes, add a slider to the figure
				posTxt = get(hTxt,'Position');
				set_slider(hFig, hTxt, posTxt, ceil(extent(4)));
				figpos(3) = figpos(3) + 15;	% Extend fig width by the slider width (15 pixels)
				set(hFig, 'Pos', figpos)
			else
				% Previous size estimates failed. Estimate again.
				set([hFig,hTxt], 'unit', 'pix')
				figpos = get(hFig, 'Pos');
				extent = get(hTxt,'Extent');
				extent(4) = extent(4) + 30;							% Need to take button size into account
				figpos(4) = round(extent(4)+10);					% New fig height
				figpos = getnicelocation(figpos, 'pixels');			% Reposition again
				set(hFig, 'Pos', figpos)
				set(hTxt, 'Pos', [posTxt(1:3) extent(4)])			% Update text position after fig resizing
				set([hFig,hTxt], 'unit', 'norm')
			end
		end

		set(hFig, 'Vis', 'on')
		if (nargout)	handFig = hFig;		end

    case 'add'
		hTxt = findobj(hFig,'Style','Text');
		hSlider = findobj(hFig,'Style','slider');
		Txt = get(hTxt,'String');
		set(hTxt, 'String',[Txt; texto])
		posTxt = get(hTxt,'Position');
		extent = get(hTxt,'Extent');
		
		if (extent(4) > 1 && isempty(hSlider))		% Text to big to fit. Add a slider to the figure
			old_u = get(hFig,'Units');		set(hFig,'Units','pixel')
			set_slider(hFig, hTxt, posTxt, (extent(4)))
			figpos = get(hFig, 'Pos');
			figpos(3) = figpos(3) + 15;	% Extend fig width by the slider width (15 pixels)
			set(hFig, 'Pos', figpos),		set(hFig, 'Units',old_u)
		elseif (extent(4) > 1)						% Slider already exists; scroll it down
			scal = get(hSlider, 'Max');
			if (scal < extent(4))
				scal = ceil(extent(4));			% Round up to next multiple of text display size
				set(hSlider, 'Max', scal);
			end
			slid_val = scal - (extent(4) - 1);
			set(hSlider,'Value',slid_val)
		end
		if (nargout)	handFig = hFig;		end

    case 'close'
		delete(hFig);
end

% ------------------------------------------------------------------------------------	
function [figpos, figName, posTxt, movepos, bgcolor, fwcolor, Font, addButt, winMaxH, modal] = parse_inputs(texto, varargin)
% Parse inputs and compute also fig and text sizes

	if ( rem(numel(varargin), 2) )
		error('msgwin:parse_input','Property value/number must come in pairs')
	end
	% Assign defaults
	win_width = 0;		win_height = 0;		movepos = [];	fwcolor = 'k';	bgcolor = [.95 .95 .95];
	addButt = false;	modal = false;
	figName = 'Message window';
	Font.Size = 9;
	Font.Weight = 'demi';
	Font.Name = 'Helvetica';
	
	for (k = 1:numel(varargin))
		if ( ischar(varargin{k}) )
			switch lower(varargin{k})
				case 'figname',		figName = varargin{k+1};
				case 'width',		win_width = varargin{k+1};
				case 'height',		win_height = varargin{k+1};
				case 'position',	movepos = varargin{k+1};
				case 'bgcolor',		bgcolor = varargin{k+1};
				case 'fwcolor',		bgcolor = varargin{k+1};
				case 'fontweight',	Font.Weight = varargin{k+1};
				case 'fontsize',	Font.Size = varargin{k+1};
				case 'fontname',	Font.Name = varargin{k+1};
				case 'button',		addButt = true;
				case 'modal',		modal = true;
			end
		end
	end

	screensize = get(0,'ScreenSize');
	if (~win_width)
		defFigPos = get(0,'DefaultfigurePosition');
		win_width = defFigPos(3);
	end
	if (~win_height)
		numLines = size(texto, 1);
		win_height = numLines * 15;
		winMaxH = round(screensize(4)*.9);		% Max win height is 90% of screen height
		win_height = min(win_height, winMaxH);
	else
		winMaxH = win_height;
	end
	if (~isempty(movepos) && ~ischar(movepos))
		error('msgwin:parse_input','POSITIOIN Property must be a char string')
	end
	
	% Calculate figure size
	figpos = [screensize(3)-5-win_width 40 win_width win_height];
	if (screensize(4) < 800),	figpos(2) = 20;		end

	% Calculate text position in figure
	bord  = 10;
 	posTxt = [bord bord/5 win_width-2*bord win_height-2*bord/5];
	
	if (addButt),	figpos(4) = figpos(4) + 30;		end

% ------------------------------------------------------------------------------------	
function figSize = getnicelocation(figSize, figUnits)
% adjust the figure position to fit nicely over GCBF or into the upper 3rd of the screen
%  Mofified from Matlab's getnicedialoglocation function

	hFigParent = gcbf;
	propName = 'Position';
	if (isempty(hFigParent))
		hFigParent = 0;
		propName = 'ScreenSize';
	end

	old_u = get(hFigParent,'Units');
	set(hFigParent,'Units',figUnits);
	parentSize = get(hFigParent,propName);
	set(hFigParent,'Units',old_u);

	figSize(1) = parentSize(1) + 1/2*(parentSize(3) - figSize(3));
	figSize(2) = parentSize(2) + 2/3*(parentSize(4) - figSize(4));

	% Now check if window does go out the ceiling
	SS = get(0,'ScreenSize');
	if ( (figSize(2) + figSize(4) + 30) > SS(4))
		figSize(2) = figSize(2) - (figSize(2) + figSize(4) + 30 - SS(4));	% Longer but clearer
	end

% ------------------------------------------------------------------------------------	
function set_slider(hFig, hTxt, posTxt, scal)
% 
	set(hFig, 'Units', 'pixels')
	pos = get(hFig, 'Pos');
	sliderW = 15 / pos(3);			% Make the slider 15 pixels wide
	pos = [1-sliderW 0 sliderW 1];
	if (scal > 1)	scal = scal - 1;	end
	cb_slide_step = {@slide_step, hTxt, posTxt};
	uicontrol(hFig,'style','slider','unit','normalized','pos',pos,...
        'call',cb_slide_step,'min',0,'max',scal,'Value',scal);

% ------------------------------------------------------------------------------------	
function slide_step(obj,event, h, pos)
	maxVal = get(obj, 'Max');
	new_pos = [pos(1:3) maxVal-get(obj,'value')+1];
	set(h,'Position',new_pos)

% ------------------------------------------------------------------------------------	
function hFig = findmsgwin(tagFig)
% H = findmsgwin('tagFig') returns the fig handle of the figure whose tag is 'tagFig'
	showBak = get(0,'ShowHiddenHandles');
	set(0,'ShowHiddenHandles','on');
	hFig = findobj(get(0,'Children'),'flat', 'tag',tagFig);
	set(0,'ShowHiddenHandles',showBak);

% --------------------------------------------------------------------
function plotMiMa(obj,evt)
% Plot the min/max locations on the figure that called message_win.
% But for that the necessary info must have been stored in this figure's
% 'UserData'. That is responsability of the calling function.

	ud = get(get(obj,'Parent'),'UserData');
	if (isempty(ud))
		errordlg('The UserData container is empty, so I don''t know where Min/Max are.','Error')
		return
	end
	if (~ishandle(ud.hAxe))
		errordlg('Invalid Parent axes handle. Did you kill that Window?','Error')
		return
	end
	h1 = line('XData',ud.p1(1),'YData',ud.p1(2),'parent',ud.hAxe, ...
		'Marker','p','MarkerFaceColor','r','MarkerEdgeColor','y','MarkerSize',10,'Tag','Symbol');
	h2 = line('XData',ud.p2(1),'YData',ud.p2(2),'parent',ud.hAxe, ...
		'Marker','p','MarkerFaceColor','b','MarkerEdgeColor','y','MarkerSize',10,'Tag','Symbol');

	draw_funs(h1,'DrawSymbol')
	draw_funs(h2,'DrawSymbol')

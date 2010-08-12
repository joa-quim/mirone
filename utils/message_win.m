function handFig = message_win(option, texto, varargin)
%MESSAGE_WIN  Message function.
%
%   Usage:  message_win('create','Some text')   for the first call
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

%	Copyright (c) 2004-2010 by J. Luis
%
%	This program is free software; you can redistribute it and/or modify
%	it under the terms of the GNU General Public License as published by
%	the Free Software Foundation; version 2 of the License.
%
%	This program is distributed in the hope that it will be useful,
%	but WITHOUT ANY WARRANTY; without even the implied warranty of
%	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%	GNU General Public License for more details.
%
%	Contact info: w3.ualg.pt/~jluis/mirone
% --------------------------------------------------------------------

tagFig = 'Wdmsgfig';
tagTxt = 'textTag';

hFig = findmsgwin(tagFig);
if ( isempty(hFig) && strcmp(option, 'add') )
	option = 'create';
end

if (~isa(texto,'cell')),		texto = cellstr(texto);	end

switch option
	case 'create'
		[figpos, figName, posTxt, movepos, bgcolor, fwcolor, Font, add_button] = parse_inputs(texto, varargin{:});
		if (isempty(movepos))			% Reposition figure on screen
			figpos = getnicelocation(figpos, 'pixels');
		end
		hFig = figure('MenuBar','none', 'Name',figName, 'HandleVisibility', 'off', ...
				'Vis','on', 'Unit','pix', 'Pos',figpos,...
				'Color',bgcolor, 'NumberTitle','off', 'DoubleBuffer','on', 'Tag',tagFig);
		hTxt = uicontrol('Parent',hFig, 'Style','text', ...
				'Unit','pix', 'Pos', posTxt, ...
				'HorizontalAlignment','left',...
				'FontWeight',Font.Weight, 'FontSize',Font.Size, 'FontName',Font.Name,...
				'BackgroundColor',bgcolor, 'ForegroundColor',fwcolor, ...
				'Tag',tagTxt);
		if (add_button)
			posBut = [(figpos(3)/2 - 40) 4 80 21];
			hBut = uicontrol('Parent',hFig, 'Style','pushbutton', ...
				'Unit','pix', 'Position', posBut, ...
				'Str', 'Plot min/max', ...
				'Call', @plotMiMa, ...
				'Tag','plotaMiMa');
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
			if (add_button)						% We need to recenter the button too
				posBut(1) = (figpos(3)/2 - 40);
				set(hBut, 'Pos', posBut)
			end
			set([hFig,hTxt], 'unit', 'norm')
		end

		if ( extent(4) > 1 || add_button)		% Text too big to fit in?
			if (figpos(4) > 500)				% Yes, add a slider to the figure
				posTxt = get(hTxt,'Position');
				set_slider(hFig, hTxt, posTxt, ceil(extent(4)))
			else
				% Previous size estimates failed. Estimate again.
				set([hFig,hTxt], 'unit', 'pix')
				figpos = get(hFig, 'Pos');
				extent = get(hTxt,'Extent');
				if (add_button)		extent(4) = extent(4) + 30;		end	% Need to take button size into account
				figpos(4) = round(extent(4)+10);						% New fig height
				figpos = getnicelocation(figpos, 'pixels');				% Reposition again
				set(hFig, 'Pos', figpos)
				set(hTxt, 'Pos', [posTxt(1:3) extent(4)])				% Update text position after fig resizing
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
			set_slider(hFig, hTxt, posTxt, ceil(extent(4)))
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
function [figpos, figName, posTxt, movepos, bgcolor, fwcolor, Font, add_button] = parse_inputs(texto, varargin)
% Parse inputs and compute also fig and text sizes

	if ( rem(numel(varargin), 2) )
		error('msgwin:parse_input','Property value/number must come in pairs')
	end
	% Assign defaults
	win_width = 0;		win_height = 0;		movepos = [];	fwcolor = 'k';	bgcolor = [.95 .95 .95];
	add_button = false;
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
				case 'button',		add_button = true;
			end
		end
	end

	if (~win_width)
		defFigPos = get(0,'DefaultfigurePosition');
		win_width = defFigPos(3);
	end
	if (~win_height)
		numLines = size(texto, 1);
		win_height = numLines * 15;
		win_height = min(win_height, 500);
	end
	if (~isempty(movepos) && ~ischar(movepos))
		error('msgwin:parse_input','POSITIOIN Property must be a char string')
	end
	
	% Calculate figure size
	screensize = get(0,'ScreenSize');
	figpos = [screensize(3)-5-win_width 40 win_width win_height];
	if (screensize(4) < 800),	figpos(2) = 20;		end

	% Calculate text position in figure
	bord  = 10;
 	posTxt = [bord bord/5 win_width-2*bord win_height-2*bord/5];
	
	if (add_button),	figpos(4) = figpos(4) + 30;		end

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
	container_size = get(hFigParent,propName);
	set(hFigParent,'Units',old_u);

	figSize(1) = container_size(1)  + 1/2*(container_size(3) - figSize(3));
	figSize(2) = container_size(2)  + 2/3*(container_size(4) - figSize(4));

% ------------------------------------------------------------------------------------	
function set_slider(hFig, hTxt, posTxt, scal)
	pos = [0.97 0 .03 1];
	cb_slide_step = {@slide_step, hTxt, posTxt};
	uicontrol(hFig,'style','slider','unit','normalized','position',pos,...
        'call',cb_slide_step,'min',0,'max',scal,'Value',scal);

% ------------------------------------------------------------------------------------	
function slide_step(obj,event, h, pos)
	maxVal = get(obj, 'Max');
	new_pos = [pos(1) pos(2) pos(3) maxVal-get(obj,'value')+1.05];
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

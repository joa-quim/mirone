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
		[figpos, figName, posTxt, movepos, bgcolor, fwcolor, fWeight, fSize, fName] = parse_inputs(texto, varargin{:});
		if (isempty(movepos))			% Reposition figure on screen
			figpos = getnicelocation(figpos, 'pixels');
		end
		hFig = figure('MenuBar','none', 'Name',figName, 'HandleVisibility', 'off', ...
				'Visible','off', 'Unit','pix', 'Position',figpos,...
				'Color',bgcolor, 'NumberTitle','off', 'DoubleBuffer','on', 'Tag',tagFig);
		hTxt = uicontrol('Parent',hFig, 'Style','text', ...
				'Unit','pix', 'Position', posTxt, ...
				'HorizontalAlignment','left',...
				'FontWeight',fWeight, 'FontSize',fSize, 'FontName',fName,...
				'BackgroundColor',bgcolor, 'ForegroundColor',fwcolor, ...
				'Tag',tagTxt);
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
			set(hFig, 'Pos', [figpos(1:2) extent_pix(3)+20 figpos(4)])
			set([hFig,hTxt], 'unit', 'norm')
		end

		if ( extent(4) > 1 )					% Text too big to fit in?
			if (figpos(4) > 500)				% Yes, add a slider to the figure
				posTxt = get(hTxt,'Position');
				set_slider(hFig, hTxt, posTxt, ceil(extent(4)))
			else
				% Previous size estimates failed. Estimate again.
				set([hFig,hTxt], 'unit', 'pix')
				figpos = get(hFig, 'Pos');
				extent = get(hTxt,'Extent');
				figpos(4) = round(extent(4)+10);	% New fig height
				figpos = getnicelocation(figpos, 'pixels');		% Reposition again
				set(hFig, 'Pos', figpos)
				set(hTxt, 'Pos', [posTxt(1:3) extent(4)])	% Update text position after figure resizing
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
function [figpos, figName, posTxt, movepos, bgcolor, fwcolor, fWeight, fSize, fName] = parse_inputs(texto, varargin)
% Parse inputs and compute also fig and text sizes

	if ( rem(numel(varargin), 2) )
		error('msgwin:parse_input','Property value/number must come in pairs')
	end
	% Assign defaults
	win_width = 0;		win_height = 0;		movepos = [];	fwcolor = 'k';	bgcolor = [.95 .95 .95];
	fWeight = 'demi';	fSize = 9;			fName = 'Helvetica';	figName = 'Message window';
	
	for (k = 1:numel(varargin))
		if ( ischar(varargin{k}) )
			switch lower(varargin{k})
				case 'figname',		figName = varargin{k+1};
				case 'width',		win_width = varargin{k+1};
				case 'height',		win_height = varargin{k+1};
				case 'position',	movepos = varargin{k+1};
				case 'bgcolor',		bgcolor = varargin{k+1};
				case 'fwcolor',		bgcolor = varargin{k+1};
				case 'fontweight',	fWeight = varargin{k+1};
				case 'fontsize',	fSize = varargin{k+1};
				case 'fontname',	fName = varargin{k+1};
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

% ------------------------------------------------------------------------------------	
function figure_size = getnicelocation(figure_size, figure_units)
% adjust the specified figure position to fig nicely over GCBF
% or into the upper 3rd of the screen
%  Copyright 1999-2006 The MathWorks, Inc.

	parentHandle = gcbf;
	propName = 'Position';
	if isempty(parentHandle)
		parentHandle = 0;
		propName = 'ScreenSize';
	end

	old_u = get(parentHandle,'Units');
	set(parentHandle,'Units',figure_units);
	container_size=get(parentHandle,propName);
	set(parentHandle,'Units',old_u);

	figure_size(1) = container_size(1)  + 1/2*(container_size(3) - figure_size(3));
	figure_size(2) = container_size(2)  + 2/3*(container_size(4) - figure_size(4));

% ------------------------------------------------------------------------------------	
function set_slider(hFig, hTxt, posTxt, scal)
	pos = [0.97 0 .03 1];
	cb_slide_step = {@slide_step, hTxt, posTxt};
	uicontrol(hFig,'style','slider','units','normalized','position',pos,...
        'callback',cb_slide_step,'min',0,'max',scal,'Value',scal);

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

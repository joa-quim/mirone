function hFig = message_win(option, texto, varargin)
%MESSAGE_WIN  Message function.
%
%   Usage:  message_win('create','Some text')   for the first call
%           message_win('add','More text')      for subsequent calls. When the window message
%                                               is full of text a slider will be added.

%   Joaquim Luis

%------------------------
tagFig = 'Wdmsgfig';
tagTxt = 'textTag';

hFig = findmsgwin(tagFig);
if ( isempty(hFig) && strcmp(option, 'add') )
	option = 'create';
end

if (~isa(texto,'cell')),		texto = cellstr(texto);	end

switch option
	case 'create'
		if isempty(hFig)
			[figpos, posTxt, movepos, bgcolor, fwcolor, fWeight, fSize] = parse_inputs(texto, varargin{:});
			hFig = figure('MenuBar','none', 'Name','Message window', 'HandleVisibility', 'off', ...
					'Visible','off', 'Unit','pixels', 'Position',figpos,...
					'Color',bgcolor, 'NumberTitle','off', 'Tag',tagFig);
			hTxt = uicontrol('Parent',hFig, 'Style','text', ...
					'Units','pixels', 'Position', posTxt, ...
					'HorizontalAlignment','left',...
					'FontWeight',fWeight, 'FontSize',fSize, ...
					'BackgroundColor',bgcolor, 'ForegroundColor',fwcolor, ...
					'Tag',tagTxt);
			if (~isempty(movepos)),		movegui(hFig, movepos),		end		% Reposition figure on screen
			set([hFig,hTxt],'units','normalized','Visible','on');
		else		% User said 'create' but figure already exists
			hTxt = findobj(hFig,'Tag',tagTxt);
			figure(hFig);
        end
		set(hTxt, 'String',texto)
		extent = get(hTxt,'Extent');
		if (isempty(extent)),	extent = zeros(1,4);	end		% Protection against empty texts
		if ( extent(4) > 1 )				% Text to big to fit. Add a slider to the figure
			posTxt = get(hTxt,'Position');
			set_slider(hFig, hTxt, posTxt, ceil(extent(4)))
		end

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

    case 'close'
		delete(hFig);
end

% ------------------------------------------------------------------------------------	
function [figpos, posTxt, movepos, bgcolor, fwcolor, fWeight, fSize] = parse_inputs(texto, varargin)
% Parse inputs and compute also fig and text sizes

	if ( rem(numel(varargin), 2) )
		error('msgwin','Property value/number must come in pairs')
	end
	% Assign defaults
	win_width = 0;		win_height = 0;		movepos = [];	bgcolor = 'w';	fwcolor = 'k';
	fWeight = 'bold';	fSize = 8;
	
	for (k = 1:numel(varargin))
		if ( ischar(varargin{k}) )
			switch lower(varargin{k})
				case 'width',		win_width = varargin{k+1};
				case 'height',		win_height = varargin{k+1};
				case 'position',	movepos = varargin{k+1};
				case 'bgcolor',		bgcolor = varargin{k+1};
				case 'fwcolor',		bgcolor = varargin{k+1};
				case 'fontweight',	fWeight = varargin{k+1};
				case 'fontsize',	fSize = varargin{k+1};
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
		%win_height = 500;
	end
	if (~isempty(movepos) && ~ischar(movepos))
		error('msgwin','POSITIOIN Property must be a char string')
	end
	
	% Calculate figure size
	screensize = get(0,'ScreenSize');
	figpos = [screensize(3)-5-win_width 40 win_width win_height];
	if (screensize(4) < 800),	figpos(2) = 20;		end

	% Calculate test position in figure
	bord  = 10;
 	posTxt = [bord bord/5 win_width-2*bord win_height-2*bord/5];
	
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

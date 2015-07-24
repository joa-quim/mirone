function varargout = resizetrue(handles, opt, axis_t, opt_pad)
%   RESIZETRUE Adjust display size of image.
%   This results in the display having one screen pixel for each image pixel (if not zoomed).
%
%	OPT      A empty var or: a 1x2 vector with the desired image size OR a scalar with dx/dy anisotropy
%
%	AXIS_T   A string with either 'xy' (sets 'YDir','normal') or 'off' (sets axes invisible)
%
%	OPT_PAD  Optional arg with the [width height] in pixels of a zone to be left clear on the
%            left side of the image. Ignored if == 0. Used only by TINTOL
%
%   This function is distantly rooted on TRUESIZE, but havily hacked to take into account
%   the left and bottom margins containing (if they exist) Xlabel & Ylabel, etc...

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

	if (nargin < 4),	opt_pad = 0;	end

	% If we already have a colorbar, remove it
	if (strcmp(get(handles.PalIn,'Check'),'on'))
		delete(get(handles.PalIn,'Userdata'))
		set(handles.PalIn,'Checked','off')
	end
	if (strcmp(get(handles.PalAt,'Check'),'on'))
		delete(get(handles.PalAt,'Userdata'))
		set(handles.PalAt,'Checked','off')
	end

	if (strcmp(axis_t,'xy')),			set(handles.axes1,'YDir','normal')
	elseif (strcmp(axis_t,'off')),		set(handles.axes1,'Visible','off')
	else    warndlg('Warning: Unknown axes setting in RESIZETRUE','Warning')
	end

	hFig = handles.figure1;
	[hAxes, hImg, msg] = ParseInputs(hFig);
	if (~isempty(msg)),		errordlg(msg,'Error'),	return,		end

	% When we open a new image to replace a previous existing one, the new image size may
	% very well be different and the status bar must be relocated. It is easier to just
	% remove the old status bar and rebuild it (later down) again.
	handsStBar = getappdata(hFig,'CoordsStBar');
	if (~isempty(handsStBar))
		delete(handsStBar);         rmappdata(hFig,'CoordsStBar');
		set(0,'CurrentFigure',hFig)     % This may be need if another figure (e.g. a warning figure) is the gcf
		pixval_stsbar('exit')
	end

	xfac = 1;
	if (handles.geog && handles.scale2meanLat)
		xfac = abs(cos(sum(handles.head(3:4)) / 2 * pi/180));
	end

	DAR = [1 1 1];      imSize = [];
	if (~ischar(opt) && numel(opt) == 2)                % imSize was transmited in input (e.g. histograms)
		imSize = opt;
		opt = 'fixed_size';
	elseif (~ischar(opt) && numel(opt) == 1)            % Case of anysotropic dx/dy
		oldDAR = get(hAxes, 'DataAspectRatio');
		% Attempt heuristic to try to discriminate between same units on both axes or different units
		if ((oldDAR(1) / oldDAR(2) < 8) && (oldDAR(1) / oldDAR(2) > 1/8))	% Assume same units
			DAR(2) = xfac;
			opt = sprintf('adjust_size_%.12f', opt * xfac);
		else
			DAR(2) = opt * xfac;
			if (opt > 1),	DAR(2) = 1/opt * xfac;	end
		end
	end
	if (isempty(opt) && xfac ~= 1)                      % Case of isotropic geog grid rescaled to mean lat
		DAR(2) = xfac;
		opt = sprintf('adjust_size_%.12f', xfac); 
	end

	if strcmp(opt,'sCapture')
		set(hAxes,'Visible','off');
		set(get(hAxes,'Title'),'Visible','on');
		delete(findobj(hFig,'Tag','sbAxes'))
	elseif strcmp(opt,'after_sCapture')
		set(hAxes,'Visible','on');
	end

	resize(hAxes, hImg, imSize, opt, handles.withSliders, handles.oldSize(2,:), opt_pad);
	set(hAxes, 'DataAspectRatio', DAR);		% Need to set it explicitly because of compiler bugginess

	if (nargout)		% Compute magnification ratio
		imgWidth  = size(get(hImg, 'CData'), 2);
		imgHeight = size(get(hImg, 'CData'), 1);
		axUnits = get(hAxes, 'Units');       set(hAxes, 'Units', 'pixels');
		axPos = get(hAxes,'Pos');            set(hAxes, 'Units', axUnits);
		varargout{1} =  round((axPos(3) / imgWidth + axPos(4) / imgHeight) * 0.5 * 100);
	end

%--------------------------------------------------------------------------------
function [hAxes,hImg,msg] = ParseInputs(varargin)

	msg = '';   hAxes = [];      hImg = [];

	if (~ishandle(varargin{1}) || ~strcmp(get(varargin{1},'type'),'figure'))
		msg = 'FIG must be a valid figure handle';		return
	else
		hFig = varargin{1};
	end
	hAxes = get(hFig, 'CurrentAxes');
	if (isempty(hAxes));
		msg = 'Current figure has no axes';		return
	end

	% Find all the images and texturemapped surfaces in the current figure.  These are the candidates.
	h = [findobj(hAxes, '-depth',1, 'Type', 'image');
		findobj(hAxes, '-depth',1, 'Type', 'surface', 'FaceColor', 'texturemap')];

	if (isempty(h))
		msg = 'No images or texturemapped surfaces in the figure';	return
	end

	% Start with the first object on the list as the initial candidate.
	% If it's not in the current axes, look for another one that is.
	hImg = h(1);
	if (get(hImg,'Parent') ~= hAxes)
		for k = 2:numel(h)
			if (get(h(k),'Parent') == hAxes),    hImg = h(k);    break;        end
		end
	end

%--------------------------------------------
function resize(hAxes, hImg, imSize, opt, withSliders, firstFigSize, pad_left)
% Resize figure containing a single axes object with a single image.
%
% PAD_LEFT   Optional arg with the width in pixels of a zone to be left clear on the
%            left side of the image. Ignored if == 0.

	if (nargin < 7),	pad_left = 0;		end

	hFig = get(hAxes, 'Parent');
	set(hAxes,'Units','normalized','Position',[0 0 1 1]) % Don't realy understand why, but I need this

	if (~pad_left(1))			% Check if this figure in padded already
		am_I_padded = getappdata(hFig,'SidePadded');
		if (~isempty(am_I_padded)),		pad_left = am_I_padded;		end
	end

	if (isempty(imSize))    % Get image's dimensions
		imgWidth  = size(get(hImg, 'CData'), 2);
		imgHeight = size(get(hImg, 'CData'), 1);
	else
		imgWidth  = imSize(2);
		imgHeight = imSize(1);
	end

	if (isa(opt,'char') && strncmp(opt,'adjust_size_',12))	% We have an anisotropic dx/dy. OPT is of the form adjust_size_[aniso]
		aniso = sscanf(opt(13:end),'%f');
		if (aniso == 0),	aniso = 1;		end				% Troublematic SeaWiFS files had made this happen
		if (aniso > 1)
			imgWidth  = imgWidth * aniso;
		else
			imgHeight = imgHeight / aniso;
		end
	end

	% Screen dimensions
	rootUnits = get(0, 'Units');			set(0, 'Units', 'pixels');
	screenSize = get(0, 'ScreenSize');      screenWidth = screenSize(3)-8;    screenHeight = screenSize(4);

	% Use a trick to simplify the case when a padding zone was required on the left side.
	% Just pretend the screen is narrower to not let the extend fig size 'overflow' the true
	% screen width and at the end we are safe to resize the fig to the account for the required extension
	if (pad_left),		screenWidth = screenWidth - pad_left(1);	end

	minFigWidth = max(firstFigSize(3), 581);	minFigHeight = 128;      % don't try to make the display image smaller than this.

	% For small images, compute the minimum side as 60% of largest of the screen dimensions
	% Except in the case of croped images, where 512 is enough for the pushbuttons (if the croped
	% image aspect ratio permits so)
	croped = getappdata(hFig,'Croped');
	if ~isempty(croped)
        LeastImageWidth = 512;    rmappdata(hFig,'Croped');
	end
	%LeastImageSide = fix(max([screenWidth screenHeight] * 0.6));

	% Mind change. For a while I'll try this way.
	LeastImageSide = 512;

	% Check that small images are grown till one of its dimensions >= LeastImageSide
	if (imgWidth < LeastImageSide && imgHeight < LeastImageSide && ~strcmp(opt,'fixed_size'))
        if ~isempty(croped)     % Croped image
            while (imgWidth < LeastImageWidth)  % Here is enough to have 1 side
                imgWidth = imgWidth*1.05;  imgHeight = imgHeight*1.05;
            end
        else                    % Full image
            while (imgWidth < LeastImageSide && imgHeight < LeastImageSide)
                imgWidth = imgWidth*1.05;  imgHeight = imgHeight*1.05;
            end
        end
        imgWidth = fix(imgWidth);   imgHeight = fix(imgHeight);
        % Large aspect ratio figures may still need to have their size adjusted
		if (imgWidth < 512)
			while (imgWidth < LeastImageSide && imgHeight < screenHeight-50)
				imgWidth = imgWidth*1.05;  imgHeight = imgHeight*1.05;
			end
		end
	end

	set(hAxes, 'Units', 'pixels');				axPos = get(hAxes, 'Position');
	figUnits = get(hFig, 'Units');				set(hFig, 'Units', 'pixels');

% ----------------------------------------------
    hXlabel = get(hAxes,'Xlabel');				hYlabel = get(hAxes,'Ylabel');
    units_save = get(hXlabel,'units');
    set(hXlabel,'units','pixels');				set(hYlabel,'units','pixels');
    Xlabel_pos = get(hXlabel,'position');		Ylabel_pos = get(hYlabel,'Extent');

    % One more atempt to make any sense out of this non-sense
    tenSizeX = 0;       tenSizeY = 0;   % When axes labels have 10^n this will hold its ~ text height
    XTickLabel = get(hAxes,'XTickLabel');    XTick = get(hAxes,'XTick');
	if (isa(XTickLabel,'cell')),		XTickLabel = XTickLabel{end};		end		% In Octave it is
    if (XTick(end) ~= 0)				% See that we do not devide by zero
		test_tick = XTick(end);			test_tick_str = sscanf(XTickLabel(end,:),'%f');
	else								% They cannot be both zero
		test_tick = XTick(end-1);		test_tick_str = sscanf(XTickLabel(end-1,:),'%f');
    end
    if ( test_tick_str && (test_tick_str / test_tick < 0.1) )	% Also test if first guy is != 0
		% We have a 10 power. That's the only way I found to detect
		% the presence of this otherwise completely ghost text.
		tenSizeX = 1;       % Take into account the 10 power text size when creating the pixval stsbar
    end

    % OK, here the problem is that YTickLabel still does not exist (imgHeight +- 2 or 3)
    set(hAxes, 'Position', axPos+[0 -500 0 500]);        % So, use this trick to set it up
    YTickLabel = get(hAxes,'YTickLabel');    YTick = get(hAxes,'YTick');
	if (isa(YTickLabel,'cell')),		YTickLabel = YTickLabel{end};		end		% In Octave it is
    if (YTick(end) ~= 0)				% See that we do not devide by zero
        test_tick = YTick(end);			test_tick_str = sscanf(YTickLabel(end,:),'%f');
	else								% They cannot be both zero
        test_tick = YTick(end-1);		test_tick_str = sscanf(YTickLabel(end-1,:),'%f');
    end
    if ( test_tick_str / test_tick < 0.1 )
        tenSizeY = 20;
    end

	% assume figure decorations are ?? pixels (!!)
	figBottomBorder = 30;       figTopBorder = 80;
	figTopBorder = figTopBorder + tenSizeY;

    sldT = 0;
    if (withSliders),       sldT = 9;      end      % Slider thickness

	% What are the gutter sizes?
	gutterLeft = max(axPos(1) - 1, 0);

	if (gutterLeft > 0)
        defAxesPos = get(0,'DefaultAxesPosition');
        gutterWidth  = round((1 - defAxesPos(3)) * imgWidth / defAxesPos(3));
        gutterHeight = round((1 - defAxesPos(4)) * imgHeight / defAxesPos(4));
        newFigWidth  = imgWidth + gutterWidth;
        newFigHeight = imgHeight + gutterHeight;
	else
        newFigWidth = imgWidth;				newFigHeight = imgHeight;
	end
	while ((newFigWidth > screenWidth) || ((newFigHeight + figBottomBorder + figTopBorder) > (screenHeight - 40)))
        imgWidth  = imgWidth * 0.98;		imgHeight  = imgHeight * 0.98;
        newFigWidth = newFigWidth * 0.98;	newFigHeight = newFigHeight * 0.98;
	end
	imgWidth  = round(imgWidth);			imgHeight  = round(imgHeight);
	newFigWidth = round(newFigWidth);		newFigHeight = round(newFigHeight);

    old_FU = get(hAxes,'FontUnits');		set(hAxes,'FontUnits','points')
    FontSize = get(hAxes,'FontSize');		set(hAxes,'FontUnits',old_FU)
    nYchars = size(YTickLabel,2);
	t = max(abs(YTick));
	if (t - fix(t) == 0 && ~tenSizeY),		nYchars = nYchars + 2;		end
	% This is kitchen sizing, but what else can it be done with such can of bugs?
	Ylabel_pos(1) = max(abs(Ylabel_pos(1)), nYchars * FontSize * 0.8 + 2);

	if strcmp(opt,'sCapture'),		stsbr_height = 0;    
	else							stsbr_height = 20;
	end

	y_margin = abs(Xlabel_pos(2))+get(hXlabel,'Margin') + tenSizeY + stsbr_height;    % To hold the Xlabel height
	x_margin = max( abs(Ylabel_pos(1)),15 );			% To hold the Ylabel width
	if (y_margin > 70)									% Play safe. LabelPos non-sense is always ready to strike 
		y_margin = 30 + tenSizeY + stsbr_height;
	end

	topMarg = 0;
	if (~tenSizeY),     topMarg = 7;    end					% To account for Ylabels exceeding image height
	axVisible = strcmp(get(hAxes,'Visible'),'on');
	if (~axVisible)				% No Labels, give only a 20 pixels margin to account for Status bar
		x_margin = 0;   y_margin = stsbr_height;
		topMarg  = 0;
	elseif (minFigWidth - x_margin > imgWidth + x_margin)	% Image + x_margin still fits inside minFigWidth
		x_margin = 0;
	end
	set(hXlabel,'units',units_save);     set(hYlabel,'units',units_save);

	newFigWidth  = max(newFigWidth + x_margin, minFigWidth);
	if (newFigWidth >= screenWidth)     % Larger than screen. The == isn't allowed either due to the 'elastic' thing
		newFigWidth = screenWidth - 0;
		imgWidth = screenWidth - x_margin;
	end
	newFigHeight = max(newFigHeight, minFigHeight) + y_margin + topMarg;

	figPos = get(hFig, 'Position');
	figPos(1) = max(1, figPos(1) - floor((newFigWidth  - figPos(3))/2));
	figPos(2) = max(1, figPos(2) - floor((newFigHeight - figPos(4))/2));
	figPos(3) = newFigWidth;
	figPos(4) = newFigHeight;

	% Figure out where to place the axes object in the resized figure
	gutterWidth  = newFigWidth  - imgWidth;
	gutterHeight = newFigHeight - imgHeight;
	gutterLeft   = floor(gutterWidth/2)  + x_margin/2;
	gutterBottom = floor(gutterHeight/2) + y_margin/2;
    
	axPos(1) = gutterLeft;      axPos(2) = gutterBottom - tenSizeY;
	axPos(3) = imgWidth;		axPos(4) = imgHeight;

	H = 22;		% The status bar (simulated box at bottom) height
	if (~axVisible && ((newFigHeight - H) > (1.25 * imgHeight)))	% A sign of a potentialy large aspect ratio
		xLim = get(hAxes, 'XLim');		yLim = get(hAxes, 'YLim');
		aspectWH = diff(xLim) / diff(yLim);
		if (aspectWH > 3)			% Yes, large aspect ratio. Try to make it a bit more usable (when zooming)
			ampFact = (1 + (newFigHeight / imgHeight - 1));
			newH = round(axPos(4) * ampFact) - H;
			axPos(2) = round(axPos(2) - (newH - axPos(4)) / 2);
			axPos(4) = round((axPos(4) - H) * ampFact);			
			yLim(2) = yLim(2) * ampFact * 1;
			set(hAxes, 'YLim', yLim)
		end
	end

	% --------- See if we need to apply further resizings bcause of left padding -----------------
	shift_frame_x = 0;	% Amount used to shift the status bar axes when pad_left != 0
	if (pad_left(1))
		dx = round(gutterWidth / 2 - x_margin - (sldT + 2));		% Slider thickness + 2
		if (dx > 0)		% We have unused space left on the right side. Use it
			axPos(1)  = axPos(1) + pad_left(1);
			figPos(3) = figPos(3) + pad_left(1) - dx;
			shift_frame_x = pad_left(1) + dx;
		else
			axPos(1) = axPos(1) + pad_left(1);
			figPos(3) = figPos(3) + pad_left(1);
		end
		if ((figPos(1) + figPos(3)) > screenSize(3))	% Fig is partially out on the right side
			figPos(1) = screenSize(3) - figPos(3);
		end
		setappdata(hFig,'SidePadded',pad_left)
	end
	% ----------------------------------------------------------------------------------------------

	% Force the window to be in the "north" position. 78 is the height of the blue Bar + ...
	figPos(2) = screenHeight - figPos(4) - 78;
	if (pad_left(1) && figPos(4) < pad_left(2))			% For the time being only TINTOL triggers this
		figPos(4) = pad_left(2);
		figPos(2) = screenHeight - figPos(4) - 78;
		offset_y = -axPos(4)/2 + pad_left(2)/2 - 20;
		axPos(2) = axPos(2) + offset_y;					% Aproximate but should be good enough
		hSliders = getappdata(hAxes,'SliderAxes');		% Need to move the horizontal slider as well
		pos = get(hSliders(1),'Pos');
		pos(2) = pos(2) + offset_y;
		set(hSliders(2), 'Pos', pos)	% THIS IS NOT WORKING
	end
	axPos = round(axPos);
	set(hFig, 'Position', figPos);
	set(hAxes, 'Position', axPos);

	hFanim = [];		% To play safe while I don't settle the 'sCapture' issue once and for all
	if ~strncmp(opt,'sCap',4)		% sCapture. I'm not sure this is used anymore
        %-------------- This section simulates a box at the bottom of the figure
        sbPos(1) = 1 + shift_frame_x;               sbPos(2) = 1;
        sbPos(3) = figPos(3) -2 - shift_frame_x;    sbPos(4) = H-1;
		hFanim = findobj(hFig,'-depth',1,'Style', 'frame', 'tag','cinanima');	% This guy exists only whith animations
		if (~isempty(hFanim))
			pos = get(hFanim, 'Pos');	sbPos(3) = pos(1)-1;		% Make the sbAxes shorter to not overlap with the frame
		end
        h = axes('Parent',hFig,'Box','off','Visible','off','Tag','sbAxes','Units','Pixels',...
			'Position',sbPos,'XLim',[0 sbPos(3)],'YLim',[0 H-1]);
        tenXMargin = 1;
        if (tenSizeX),     tenXMargin = 30;     end
        hFieldFrame = createframe(h, [1 (sbPos(3) + 2 - tenXMargin)], H);
        setappdata(hFig,'CoordsStBar',[h hFieldFrame]);  % Save it for use in ...
        set(hFieldFrame,'Visible','on')
        set(h,'HandleVisibility','off')
        if (withSliders),        setSliders(hFig, hAxes, figPos, axPos, sldT, H);    end   
	end
	%------------------------------------

	% Restore the units
	set(hFig, 'Units', figUnits);
	set(hAxes, 'Units', 'normalized');   % So that resizing the Fig also resizes the image
	set(0, 'Units', rootUnits);

	if ~strcmp(opt,'sCapture'),   pixval_stsbar(hFig);  end
	if(~isempty(hFanim))
		set(getappdata(hFanim,'AnimaFamily'), 'units', 'normalized')	% So that they remain in their corner uppon resizing
	end

%--------------------------------------------------------------------------
function hFrame = createframe(ah, fieldPos, H)
% Creates a virtual panel surrounding the field starting at fieldPos(1) and
% ending end fieldPos(2) pixels. AH is the sb's handle (axes).
% It returns a handle array designating the frame.
	
	from = fieldPos(1);     to = fieldPos(2);
	lightColor = [0.9 0.891509433962 0.874528301887];
	darkColor  = [0.4 0.392452830189 0.377358490566];
	
	hFrame(1) = line([from to],[H-2 H-2],'Color',darkColor,'Visible','off','Tag','Sts_T','parent',ah);    % Top line
	hFrame(2) = line([from from],[1 H-2],'Color',darkColor,'Visible','off','Tag','Sts_L','parent',ah);    % Left line
	hFrame(3) = line([from+1 to-1],[1 1],'Color',lightColor,'Visible','off','Tag','Sts_B','parent',ah);   % Bottom line
	hFrame(4) = line([to-1 to-1],[1 H-2],'Color',lightColor,'Visible','off','Tag','Sts_R','parent',ah);   % Right line

%--------------------------------------------------------------------------
function figPos = setSliders(hFig, hAxes, figPos, axPos, sldT, H)
% Create a pair of sliders and register them to use in 'imscroll_j'
	gutterRight = figPos(3) - (axPos(1) + axPos(3));
	if (gutterRight < sldT+1)       % Grow Figure's width to acomodate the slider
		figPos(3) = figPos(3) + (sldT-gutterRight) + 1;
		set(hFig, 'Position', figPos);
	end

    hSliders = getappdata(hAxes,'SliderAxes');
	if (isempty(hSliders))
		sliderVer = uicontrol('Units','pixels','Style','slider','Parent',hFig, 'Visible','off',...
            'Pos',[axPos(1)+axPos(3)+2 axPos(2) sldT axPos(4)+1],'Background',[.9 .9 .9]);
		sliderHor = uicontrol('Units','pixels','Style','slider','Parent',hFig, 'Visible','off',...
            'Pos',[axPos(1) H-1 axPos(3)+1 sldT],'Background',[.95 .95 .95]);
		set(sliderHor,'Min',0,'Max',1,'Value',0,'Tag','HOR','Callback',{@slider_Cb,hAxes,'SetSliderHor'})
		set(sliderVer,'Min',0,'Max',1,'Value',0,'Tag','VER','Callback',{@slider_Cb,hAxes,'SetSliderVer'})
		%set(sliderVer,'Min',0,'Max',1,'Value',0,'Tag','VER','Callback','imscroll_j(gca,''SetSliderVer'')')
		% Register the sliders in the axe's appdata
		drawnow			% Needed because otherwise - BUG BUG - next line would change the VER slider thickness !!!!!
		setappdata(hAxes,'SliderAxes',[sliderHor sliderVer])
		imscroll_j(hAxes,'ZoomSetSliders')              % ...
	else			% We have them already. They just need to be updated
        set(hSliders(1), 'Pos',[axPos(1)+axPos(3)+1 axPos(2) sldT axPos(4)+1],'Vis','off')
        set(hSliders(2), 'Pos',[axPos(1) H-1 axPos(3)+1 sldT],'Vis','off')
	end

% -----------------------------------------------------------------------------------------
function slider_Cb(obj,evt,ax,opt)
	imscroll_j(ax,opt)

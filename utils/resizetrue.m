function varargout = resizetrue(handles, opt, axis_t)
%   RESIZETRUE Adjust display size of image.
%   This results in the display having one screen pixel for each image pixel.
%
%   This function is based on Matlab's TRUESIZE, but havily hacked to take into account
%   the left and bottom margins containing (if they exist) Xlabel & Ylabel, etc...

%	Copyright (c) 2004-2011 by J. Luis
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
	else    warndlg('Warning: Unknown axes setting in show_image','Warning')
	end

	hFig = handles.figure1;
	[axHandle, imHandle, msg] = ParseInputs(hFig);
	if (~isempty(msg));    errordlg(msg,'Error');  return;      end

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
	if (nargin == 1)
		opt = [];
	elseif (~ischar(opt) && numel(opt) == 2)            % imSize was transmited in input (e.g. histograms)
		imSize = opt;
		opt = 'fixed_size';
	elseif (~ischar(opt) && numel(opt) == 1)            % Case of anysotropic dx/dy
		DAR(2) = xfac;
		opt = ['adjust_size_' sprintf('%.12f', opt * xfac)];
	end
	if (isempty(opt) && xfac ~= 1)                      % Case of isotropic geog grid rescaled to mean lat
		DAR(2) = xfac;
		opt = ['adjust_size_' sprintf('%.12f', xfac)]; 
	end

	if strcmp(opt,'sCapture')
		set(axHandle,'Visible','off');
		set(get(axHandle,'Title'),'Visible','on');
		delete(findobj(hFig,'Tag','sbAxes'))
	elseif strcmp(opt,'after_sCapture')
		set(axHandle,'Visible','on');
	end

	Resize1(axHandle, imHandle, imSize, opt, handles.withSliders, handles.oldSize(2,:));
	set(axHandle, 'DataAspectRatio', DAR);		% Need to set it explicitly because of compiler bugginess

	if (nargout)        % Compute magnification ratio
		imageWidth  = size(get(imHandle, 'CData'), 2);
		imageHeight = size(get(imHandle, 'CData'), 1);
		axUnits = get(axHandle, 'Units');       set(axHandle, 'Units', 'pixels');
		axPos = get(axHandle,'Pos');            set(axHandle, 'Units', axUnits);
		varargout{1} =  round((axPos(3) / imageWidth + axPos(4) / imageHeight) * 0.5 * 100);
	end

%--------------------------------------------
% Subfunction ParseInputs
%--------------------------------------------
function [axHandle,imHandle,msg] = ParseInputs(varargin)

msg = '';   axHandle = [];      imHandle = [];          figHandle = [];

if (nargin == 0);    axHandle = gca;    end
    
if (nargin >= 1)
    if (~ishandle(varargin{1}) || ~strcmp(get(varargin{1},'type'),'figure'))
        msg = 'FIG must be a valid figure handle';      return;
    else
        figHandle = varargin{1};
    end
    axHandle = get(figHandle, 'CurrentAxes');
    if (isempty(axHandle));     msg = 'Current figure has no axes';     return;    end
end

if (isempty(figHandle))
    figHandle = get(axHandle, 'Parent');
end

% Find all the images and texturemapped surfaces in the current figure.  These are the candidates.
h = [findobj(figHandle, 'Type', 'image') ;
    findobj(figHandle, 'Type', 'surface', 'FaceColor', 'texturemap')];

if (isempty(h)),    msg = 'No images or texturemapped surfaces in the figure';    return;   end

% Start with the first object on the list as the initial candidate.
% If it's not in the current axes, look for another one that is.
imHandle = h(1);
if (get(imHandle,'Parent') ~= axHandle)
    for k = 2:length(h)
        if (get(h(k),'Parent') == axHandle),    imHandle = h(k);    break;        end
    end
end

%--------------------------------------------
% Subfunction Resize1
%--------------------------------------------
function Resize1(axHandle, imHandle, imSize, opt, withSliders, firstFigSize)
% Resize figure containing a single axes object with a single image.

	figHandle = get(axHandle, 'Parent');
	set(axHandle,'Units','normalized','Position',[0 0 1 1]) % Don't realy understand why, but I need this

	if (isempty(imSize))    % How big is the image?
		imageWidth  = size(get(imHandle, 'CData'), 2);
		imageHeight = size(get(imHandle, 'CData'), 1);
	else
		imageWidth  = imSize(2);
		imageHeight = imSize(1);
	end

	if (strncmp(opt,'adjust_size_',12))			% We have an anisotropic dx/dy. OPT is of the form adjust_size_[aniso]
		aniso = sscanf(opt(13:end),'%f');
		if (aniso == 0),	aniso = 1;		end				% Troublematic SeaWiFS files had made this happen
		if (aniso > 1)
			imageWidth  = imageWidth * aniso;
		else
			imageHeight = imageHeight / aniso;
		end
	end

	% Don't try to handle the degenerate case.
	if (imageWidth * imageHeight == 0),    return;  end
	
	% What are the screen dimensions?
	rootUnits = get(0, 'Units');			set(0, 'Units', 'pixels');
	screenSize = get(0, 'ScreenSize');      screenWidth = screenSize(3)-4;    screenHeight = screenSize(4);
	
	minFigWidth = max(firstFigSize(3), 581);	minFigHeight = 128;      % don't try to display a figure smaller than this.
	
	% For small images, compute the minimum side as 60% of largest of the screen dimensions
	% Except in the case of croped images, where 512 is enough for the pushbuttons (if the croped
	% image aspect ratio permits so)
	croped = getappdata(figHandle,'Croped');
	if ~isempty(croped)
        LeastImageWidth = 512;    rmappdata(figHandle,'Croped');
	end
	%LeastImageSide = fix(max([screenWidth screenHeight] * 0.6));
	
	% Mind change. For a while I'll try this way.
	LeastImageSide = 512;

	if (imageWidth < LeastImageSide && imageHeight < LeastImageSide && ~strcmp(opt,'fixed_size'))   % Augment very small images
        if ~isempty(croped)     % Croped image
            while (imageWidth < LeastImageWidth)  % Here is enough to have 1 side
                imageWidth = imageWidth*1.05;  imageHeight = imageHeight*1.05;
            end
        else                    % Full image
            while (imageWidth < LeastImageSide && imageHeight < LeastImageSide)
                imageWidth = imageWidth*1.05;  imageHeight = imageHeight*1.05;
            end
        end
        imageWidth = fix(imageWidth);   imageHeight = fix(imageHeight);
        % Large aspect ratio figures may still need to have their size adjusted
        if (imageWidth < 512)
            while (imageWidth < LeastImageSide && imageHeight < screenHeight-50)
                imageWidth = imageWidth*1.05;  imageHeight = imageHeight*1.05;
            end
        end
	end

	set(axHandle, 'Units', 'pixels');	        axPos = get(axHandle, 'Position');
	figUnits = get(figHandle, 'Units');			set(figHandle, 'Units', 'pixels');

% ---------------------------------------------
    h_Xlabel = get(axHandle,'Xlabel');			h_Ylabel = get(axHandle,'Ylabel');
    units_save = get(h_Xlabel,'units');
    set(h_Xlabel,'units','pixels');				set(h_Ylabel,'units','pixels');
    Xlabel_pos = get(h_Xlabel,'position');		Ylabel_pos = get(h_Ylabel,'Extent');

    % One more atempt to make any sense out of this non-sense
    tenSizeX = 0;       tenSizeY = 0;   % When axes labels have 10^n this will hold its ~ text height
    XTickLabel = get(axHandle,'XTickLabel');    XTick = get(axHandle,'XTick');
	if (isa(XTickLabel,'cell')),		XTickLabel = XTickLabel{end};		end		% In Octave it is
    if (XTick(end) ~= 0)				% See that we do not devide by zero
		test_tick = XTick(end);			test_tick_str = sscanf(XTickLabel(end,:),'%f');
	else								% They cannot be both zero
		test_tick = XTick(end-1);		test_tick_str = sscanf(XTickLabel(end-1,:),'%f');
    end
    if ( test_tick_str / test_tick < 0.1 )
		% We have a 10 power. That's the only way I found to detect
		% the presence of this otherwise completely ghost text.
		tenSizeX = 1;       % Take into account the 10 power text size when creating the pixval stsbar
    end

    % OK, here the problem is that YTickLabel still does not exist (imageHeight +- 2 or 3)
    set(axHandle, 'Position', axPos+[0 -500 0 500]);        % So, use this trick to set it up
    YTickLabel = get(axHandle,'YTickLabel');    YTick = get(axHandle,'YTick');
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
    if (withSliders),       sldT = 7;      end      % Slider thickness

	% What are the gutter sizes?
	gutterLeft = max(axPos(1) - 1, 0);
	nonzeroGutters = (gutterLeft > 0);

	if (nonzeroGutters)
        defAxesPos = get(0,'DefaultAxesPosition');
        gutterWidth  = round((1 - defAxesPos(3)) * imageWidth / defAxesPos(3));
        gutterHeight = round((1 - defAxesPos(4)) * imageHeight / defAxesPos(4));
        newFigWidth  = imageWidth + gutterWidth;
        newFigHeight = imageHeight + gutterHeight;
	else
        newFigWidth = imageWidth;				newFigHeight = imageHeight;
	end
	while ((newFigWidth > screenWidth) || ((newFigHeight + figBottomBorder + figTopBorder) > (screenHeight - 40)))
        imageWidth  = imageWidth * 0.98;		imageHeight  = imageHeight * 0.98;
        newFigWidth = newFigWidth * 0.98;		newFigHeight = newFigHeight * 0.98;
	end
	imageWidth  = round(imageWidth);			imageHeight  = round(imageHeight);
	newFigWidth = round(newFigWidth);			newFigHeight = round(newFigHeight);

    old_FU = get(axHandle,'FontUnits');			set(axHandle,'FontUnits','points')
    FontSize = get(axHandle,'FontSize');		set(axHandle,'FontUnits',old_FU)
    nYchars = size(YTickLabel,2);
	t = max(abs(YTick));
	if (t - fix(t) == 0 && ~tenSizeY),			nYchars = nYchars + 2;		end
	% This is kitchen sizing, but what else can it be done with such can of bugs?
	Ylabel_pos(1) = max(abs(Ylabel_pos(1)), nYchars * FontSize * 0.8 + 2);

	if strcmp(opt,'sCapture'),		stsbr_height = 0;    
	else							stsbr_height = 20;
	end

	y_margin = abs(Xlabel_pos(2))+get(h_Xlabel,'Margin') + tenSizeY + stsbr_height;    % To hold the Xlabel height
	x_margin = max( abs(Ylabel_pos(1)),15 );			% To hold the Ylabel width
	if (y_margin > 70)									% Play safe. LabelPos non-sense is always ready to strike 
		y_margin = 30 + tenSizeY + stsbr_height;
	end

	topMarg = 0;
	if (~tenSizeY),     topMarg = 5;    end					% To account for Ylabels exceeding image height
	axVisible = strcmp(get(axHandle,'Visible'),'on');
	if (~axVisible)				% No Labels, give only a 20 pixels margin to account for Status bar
		x_margin = 0;   y_margin = stsbr_height;
		topMarg  = 0;
	elseif (minFigWidth - x_margin > imageWidth + x_margin)	% Image + x_margin still fits inside minFigWidth
		x_margin = 0;
	end
	set(h_Xlabel,'units',units_save);     set(h_Ylabel,'units',units_save);

	newFigWidth  = max(newFigWidth + x_margin, minFigWidth);
	if (newFigWidth >= screenWidth)     % Larger than screen. The == isn't allowed either due to the 'elastic' thing
		%x_margin = x_margin - (newFigWidth-screenWidth) - 2;    % 'Discount' the difference on x_margin
		newFigWidth = screenWidth - 0;
		imageWidth = screenWidth - x_margin;
	end
	newFigHeight = max(newFigHeight, minFigHeight) + y_margin + topMarg;

	figPos = get(figHandle, 'Position');
	figPos(1) = max(1, figPos(1) - floor((newFigWidth  - figPos(3))/2));
	figPos(2) = max(1, figPos(2) - floor((newFigHeight - figPos(4))/2));
	figPos(3) = newFigWidth;
	figPos(4) = newFigHeight;

	% Figure out where to place the axes object in the resized figure
	gutterWidth  = newFigWidth  - imageWidth;
	gutterHeight = newFigHeight - imageHeight;
	gutterLeft   = floor(gutterWidth/2)  + x_margin/2;
	gutterBottom = floor(gutterHeight/2) + y_margin/2;
    
	axPos(1) = gutterLeft;      axPos(2) = gutterBottom - tenSizeY;
	axPos(3) = imageWidth;      axPos(4) = imageHeight;

	H = 22;		% The status bar (simulated box at bottom) height
	if ( ~axVisible && ((newFigHeight - H) > (1.25 * imageHeight)) )	% A sign of a potentialy large aspect ratio
		xLim = get(axHandle, 'XLim');		yLim = get(axHandle, 'YLim');
		aspectWH = diff(xLim) / diff(yLim);
		if (aspectWH > 3)			% Yes, large aspect ratio. Try to make it a bit more usable (when zooming)
			ampFact = (1 + (newFigHeight / imageHeight - 1));
			newH = round(axPos(4) * ampFact) - H;
			axPos(2) = round(axPos(2) - (newH - axPos(4)) / 2);
			axPos(4) = round((axPos(4) - H) * ampFact);			
			yLim(2) = yLim(2) * ampFact * 1;
			set(axHandle, 'YLim', yLim)
		end
	end

	% Force the window to be in the "north" position. 73 is the height of the blue Bar + ...
	figPos(2) = screenHeight - figPos(4) - 73;
	set(figHandle, 'Position', figPos);
	set(axHandle, 'Position', axPos);

	if ~strncmp(opt,'sCap',4)		% sCapture. I'm not sure this is used anymore
        %-------------- This section simulates a box at the bottom of the figure
        sbPos(1) = 1;               sbPos(2) = 2;
        sbPos(3) = figPos(3)-2;     sbPos(4) = H-1;
        h = axes('Parent',figHandle,'Box','off','Visible','off','Tag','sbAxes','Units','Pixels',...
			'Position',sbPos,'XLim',[0 sbPos(3)],'YLim',[0 H-1]);
        tenXMargin = 1;
        if (tenSizeX),     tenXMargin = 30;     end
        hFieldFrame = createframe(h,[1 (figPos(3) - tenXMargin)],H);
        setappdata(figHandle,'CoordsStBar',[h hFieldFrame]);  % Save it for use in ...
        set(hFieldFrame,'Visible','on')
        set(h,'HandleVisibility','off')
        if (withSliders),        setSliders(figHandle, axHandle, figPos, axPos, sldT, H);    end   
	end
	%------------------------------------
	  
	% Restore the units
	set(figHandle, 'Units', figUnits);
	set(axHandle, 'Units', 'normalized');   % So that resizing the Fig also resizes the image
	set(0, 'Units', rootUnits);
	pause(0.01)				% Needed for example when creating bg on which slow things will be rendered
	
	if ~strcmp(opt,'sCapture'),   pixval_stsbar(figHandle);  end

%--------------------------------------------------------------------------
function hFrame = createframe(ah,fieldPos,H)
% Creates a virtual panel surrounding the field starting at fieldPos(1) and
% ending end fieldPos(2) pixels. ah is the sb's handle (axes).
% It returns a handle array designating the frame.
	
	from = fieldPos(1);     to = fieldPos(2);
	% col = rgb2hsv(get(fh,'Color'));       % fh was the figure's handle
	% lightColor = col;   lightColor(2) = 0.5*lightColor(2);  lightColor(3) = 0.9; lightColor = hsv2rgb(lightColor);
	% darkColor = col;    darkColor(3) = 0.4;  darkColor = hsv2rgb(darkColor);
	% This is the result of the above. I just don't want the extra burden of compiling those routines
	lightColor = [0.9 0.89150943396226 0.87452830188679];
	darkColor  = [0.4 0.39245283018868 0.37735849056604];
	
	hFrame(1) = line([from to],[H-2 H-2],'Color',darkColor,'Visible','off','Tag','Sts_T','parent',ah);    % Top line
	hFrame(2) = line([from from],[1 H-2],'Color',darkColor,'Visible','off','Tag','Sts_L','parent',ah);    % Left line
	hFrame(3) = line([from+1 to-1],[1 1],'Color',lightColor,'Visible','off','Tag','Sts_B','parent',ah);   % Bottom line
	hFrame(4) = line([to-1 to-1],[1 H-2],'Color',lightColor,'Visible','off','Tag','Sts_R','parent',ah);   % Right line

%--------------------------------------------------------------------------
function figPos = setSliders(figHandle, axHandle, figPos, axPos, sldT, H)
% Create a pair of sliders and register them to use in 'imscroll_j'
	gutterRight = figPos(3) - (axPos(1) + axPos(3));
	if (gutterRight < sldT+1)       % Grow Figure's width to acomodate the slider
		figPos(3) = figPos(3) + (sldT-gutterRight) + 1;
		set(figHandle, 'Position', figPos);
	end

    hSliders = getappdata(axHandle,'SliderAxes');
	if (isempty(hSliders))
		sliderVer = uicontrol('Units','pixels','Style','slider','Parent',figHandle,...
            'Pos',[axPos(1)+axPos(3)+1 axPos(2) sldT axPos(4)+1],'Background',[.9 .9 .9]);
		sliderHor = uicontrol('Units','pixels','Style','slider','Parent',figHandle,...
            'Pos',[axPos(1) H-1 axPos(3)+1 sldT],'Background',[.95 .95 .95]);
		set(sliderHor,'Min',0,'Max',1,'Value',0,'Tag','HOR','Callback',{@slider_Cb,axHandle,'SetSliderHor'})
		%set(sliderVer,'Min',0,'Max',1,'Value',0,'Tag','VER','Callback',{@slider_Cb,axHandle,'SetSliderVer'})
		set(sliderVer,'Min',0,'Max',1,'Value',0,'Tag','VER','Callback','imscroll_j(gca,''SetSliderVer'')')
		% Register the sliders in the axe's appdata
		setappdata(axHandle,'SliderAxes',[sliderHor sliderVer])
		imscroll_j(axHandle,'ZoomSetSliders')              % ...
	else			% We have them already. They just need to be updated
        set(hSliders(1), 'Pos',[axPos(1)+axPos(3)+1 axPos(2) sldT axPos(4)+1],'Vis','off')
        set(hSliders(2), 'Pos',[axPos(1) H-1 axPos(3)+1 sldT],'Vis','off')
	end

% -----------------------------------------------------------------------------------------
function slider_Cb(obj,evt,ax,opt)
	imscroll_j(ax,opt)

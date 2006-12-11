function resizetrue(hFig, opt)
%RESIZETRUE Adjust display size of image.
%   RESIZETRUE(FIG) uses the image height and width for [MROWS MCOLS].
%   This results in the display having one screen pixel for each image pixel.
%   If you omit the figure argument, RESIZETRUE works on the current figure.
%
%   This function is based on Matlab's TRUESIZE, but hacked to take into account
%   the left and bottom margins containing (if they exist) Xlabel & Ylabel

%   Coffeeright 2002-2003 J. Luis & mathworks

[axHandle, imHandle, colorbarHandle, imSize, resizeType, msg] = ParseInputs(hFig);
if (~isempty(msg));    errordlg(msg,'Error');  return;      end

% When we open a new image to replace a previous existing one, the new image size may
% very well be different and the status bar must be relocated. It is easier to just
% remove the old status bar and rebuild it (later down) again.
handsStBar = getappdata(hFig,'CoordsStBar');
if (~isempty(handsStBar))
    delete(handsStBar)
    rmappdata(hFig,'CoordsStBar');
    set(0,'CurrentFigure',hFig)     % This may be need if another figure (e.g. a warning figure) is the gcf
    pixval_stsbar('exit')
end

if (nargin == 1)
    opt = [];
elseif (~ischar(opt) && length(opt) == 2)            % imSize was transmited in input (e.g. histograms)
    imSize = opt;
    opt = 'fixed_size';
elseif (~ischar(opt) && length(opt) == 1)            % Case of anysotropic dx/dy
    opt = ['adjust_size_' num2str(opt,'%.8f')];
end
if strcmp(opt,'screen_capture')
    set(axHandle,'Visible','off');
    set(get(axHandle,'Title'),'Visible','on');
    delete(findobj(hFig,'Tag','sbAxes'))
elseif strcmp(opt,'after_screen_capture')
    set(axHandle,'Visible','on');
end

Resize1(axHandle, imHandle, imSize, opt);
%set(gca,'DataAspectRatio',[1 1 1])

%--------------------------------------------
% Subfunction ParseInputs
%--------------------------------------------
function [axHandle,imHandle,colorbarHandle,imSize,resizeType,msg] = ParseInputs(varargin)

imSize = [];        colorbarHandle = [];    msg = '';
axHandle = [];      imHandle = [];          figHandle = [];     resizeType = [];

if (nargin == 0);    axHandle = gca;    end
    
if (nargin >= 1)
    if ((nargin == 1) && isequal(size(varargin{1}), [1 2]))        % truesize([M N])
        figHandle = gcf;        imSize = varargin{1};
    else        % truesize(FIG, ...)
        if (~ishandle(varargin{1}) || ~strcmp(get(varargin{1},'type'),'figure'))
            msg = 'FIG must be a valid figure handle';      return;
        else
            figHandle = varargin{1};
        end
    end
    axHandle = get(figHandle, 'CurrentAxes');
    if (isempty(axHandle));     msg = 'Current figure has no axes';     return;    end
end

if (nargin >= 2)
    imSize = varargin{2};
    if (~isequal(size(imSize), [1 2]));     msg = 'REQSIZE must be a 1-by-2 vector';
        return;
    end
end

if (isempty(figHandle))
    figHandle = get(axHandle, 'Parent');
end

% Find all the images and texturemapped surfaces in the current figure.  These are the candidates.
h = [findobj(figHandle, 'Type', 'image') ;
    findobj(figHandle, 'Type', 'surface', 'FaceColor', 'texturemap')];

% If there's a colorbar, ignore it.
colorbarHandle = findobj(figHandle, 'type', 'image', 'Tag', 'TMW_COLORBAR');
if (~isempty(colorbarHandle))
    for k = 1:length(colorbarHandle),       h(h == colorbarHandle(k)) = [];    end
end

if (isempty(h)),    msg = 'No images or texturemapped surfaces in the figure';    return;   end

% Start with the first object on the list as the initial candidate.
% If it's not in the current axes, look for another one that is.
imHandle = h(1);
if (get(imHandle,'Parent') ~= axHandle)
    for k = 2:length(h)
        if (get(h(k),'Parent') == axHandle),    imHandle = h(k);    break;        end
    end
end

figKids = allchild(figHandle);
if ((length(findobj(figKids, 'flat', 'Type', 'axes')) == 1) && ...
            isempty(findobj(figKids, 'flat', 'Type', 'uicontrol', ...
            'Visible', 'on')) && (length(imHandle) == 1))
    % Resize type 1. Figure contains only one axes object, which contains only one image.
    resizeType = 1;
elseif (isempty(colorbarHandle))    % Figure contains other objects and not a colorbar.
    resizeType = 2;
else            % Figure contains a colorbar.  Do we have a one image one colorbar situation?
    if (length(colorbarHandle) > 1)  % No
        resizeType = 2;
    else
        colorbarAxes = get(colorbarHandle, 'Parent');
        uimenuKids = findobj(figKids, 'flat', 'Type', 'uimenu');
        invisUicontrolKids = findobj(figKids, 'flat', 'Type', 'uicontrol', 'Visible', 'off');
        otherKids = setdiff(figKids, [uimenuKids ; colorbarAxes ; axHandle ; invisUicontrolKids]);
        if (length(otherKids) > 1)  % No
            resizeType = 2;
        else                        % Yes, one image and one colorbar
            resizeType = 3;
        end
    end
end

%--------------------------------------------
% Subfunction Resize1
%--------------------------------------------
function Resize1(axHandle, imHandle, imSize, opt)
% Resize figure containing a single axes object with a single image.

figHandle = get(axHandle, 'Parent');

% h_Xlabel = get(axHandle,'Xlabel');      h_Ylabel = get(axHandle,'Ylabel');
% Xlabel_pos = get(h_Xlabel,'pos');       Ylabel_pos = get(h_Ylabel,'pos');
set(axHandle,'Units','normalized','Position',[0 0 1 1]) % Don't realy understand why, but I need this
% set(h_Xlabel,'pos',Xlabel_pos);         set(h_Ylabel,'pos',Ylabel_pos);

if (isempty(imSize))    % How big is the image?
    imageWidth  = size(get(imHandle, 'CData'), 2);
    imageHeight = size(get(imHandle, 'CData'), 1);
else
    imageWidth  = imSize(2);
    imageHeight = imSize(1);
end

if (length(opt) > 12 && strcmp(opt(1:11),'adjust_size'))   % We have an anisotropic dx/dy. OPT is of the form adjust_size_[aniso]
    aniso = str2double(opt(13:end));
    if (aniso > 1)
        imageWidth  = imageWidth * aniso;
    else
        imageHeight = imageHeight / aniso;
    end
end

% Don't try to handle the degenerate case.
if (imageWidth * imageHeight == 0),    return;  end

% What are the screen dimensions
screenSize = get(0, 'ScreenSize');      screenWidth = screenSize(3);    screenHeight = screenSize(4);
if ((screenWidth <= 1) || (screenHeight <= 1))
    screenWidth = Inf;    screenHeight = Inf;
end

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

axUnits = get(axHandle, 'Units');
set(axHandle, 'Units', 'pixels');
axPos = get(axHandle, 'Position');

figUnits = get(figHandle, 'Units');
rootUnits = get(0, 'Units');
set(figHandle, 'Units', 'pixels');
set(0, 'Units', 'pixels');

figLeftBorder = 10;  % assume left figure decorations are 10 pixels
figRightBorder = 10;
figBottomBorder = 30;
figTopBorder = 50;
figTopBorder = figTopBorder + 30;  % scribe hack

minFigWidth = 600;      % don't try to display a figure smaller than this.
minFigHeight = 128;

% What are the gutter sizes?
figPos = get(figHandle, 'Position');
gutterLeft = max(axPos(1) - 1, 0);

nonzeroGutters = (gutterLeft > 0);
if (nonzeroGutters)
    defAxesPos = get(0,'DefaultAxesPosition');
    gutterWidth  = round((1 - defAxesPos(3)) * imageWidth / defAxesPos(3));
    gutterHeight = round((1 - defAxesPos(4)) * imageHeight / defAxesPos(4));
    newFigWidth  = imageWidth + gutterWidth;
    newFigHeight = imageHeight + gutterHeight;
else
    newFigWidth = imageWidth;        newFigHeight = imageHeight;
end
while (((newFigWidth + figLeftBorder + figRightBorder) > screenWidth) || ...
            ((newFigHeight + figBottomBorder + figTopBorder) > (screenHeight - 40)))
    imageWidth  = round(imageWidth * 0.95);
    imageHeight = round(imageHeight * 0.95);
    newFigWidth  = round(newFigWidth * 0.95);
    newFigHeight = round(newFigHeight * 0.95);
end

% ---------------------------------------------
h_Xlabel = get(axHandle,'Xlabel');      h_Ylabel = get(axHandle,'Ylabel');
units_save = get(h_Xlabel,'units');
set(h_Xlabel,'units','pixels');         set(h_Ylabel,'units','pixels');
Xlabel_pos = get(h_Xlabel,'pos');

% One more atempt to make any sense out of this non-sense
%YTickLabel = get(axHandle,'YTickLabel');
%set(get(axHandle,'YLabel'),'String',YTickLabel)
Ylabel_pos = get(h_Ylabel,'Extent');
Ylabel_pos(1) = Ylabel_pos(4) - 5;

if (abs(Ylabel_pos(1)) < 20)    % Stupid hack, but there is a bug somewhere
    Ylabel_pos(1) = 30;
end

if (isempty(opt) || strcmp(opt(1:5),'fixed') || strcmp(opt(1:6),'adjust')) 
    setappdata(figHandle,'Backup_LabelPos',[Xlabel_pos Ylabel_pos])
elseif strcmp(opt,'after_screen_capture')
    lab_tmp = getappdata(figHandle,'Backup_LabelPos');
    Xlabel_pos = [lab_tmp(1) lab_tmp(2) lab_tmp(3)];
    Ylabel_pos = [lab_tmp(4) lab_tmp(5) lab_tmp(6)];
end

if strcmp(opt,'screen_capture'),    stsbr_height = 0;
else                                stsbr_height = 20;    end

y_margin = abs(Xlabel_pos(2))+get(h_Xlabel,'Margin')+stsbr_height;  % Devera conter a altura em pixeis do Xlabel
x_margin = abs(Ylabel_pos(1))+get(h_Ylabel,'Margin');  % Devera conter a largura em pixeis do Ylabel
y_margin = max(y_margin,40);            % Another hack due to the LabelPos non-sense
topMarg = 5;                            % To account for Ylabels exceeding image height
if strcmp(get(axHandle,'Visible'),'off')               % No Labels, give only a 20 pixels margin to account for Status bar
    x_margin = 0;   y_margin = stsbr_height;
    topMarg = 0;
elseif (minFigWidth - x_margin > imageWidth + x_margin)% Image + x_margin still fits inside minFigWidth
    x_margin = 0;
end
set(h_Xlabel,'units',units_save);     set(h_Ylabel,'units',units_save);

newFigWidth = max(newFigWidth, minFigWidth) + x_margin;
newFigHeight = max(newFigHeight, minFigHeight) + y_margin + topMarg;

figPos(1) = max(1, figPos(1) - floor((newFigWidth - figPos(3))/2));
figPos(2) = max(1, figPos(2) - floor((newFigHeight - figPos(4))/2));
figPos(3) = newFigWidth;
figPos(4) = newFigHeight;

% Figure out where to place the axes object in the resized figure
gutterWidth = figPos(3) - imageWidth;
gutterHeight = figPos(4) - imageHeight;
gutterLeft = floor(gutterWidth/2) + x_margin/2 - 1;
gutterBottom = floor(gutterHeight/2) + y_margin/2 - 1;

axPos(1) = gutterLeft + 1;  axPos(2) = gutterBottom + 1;
axPos(3) = imageWidth;      axPos(4) = imageHeight;

% Force the window to be in the "north" position. 75 is the height of the blue Bar + ...
figPos(2) = screenHeight - figPos(4) - 75;
set(figHandle, 'Position', figPos);     set(axHandle, 'Position', axPos);

if ~strcmp(opt,'screen_capture')
    %-------------- This section simulates a box at the bottom of the figure
    H = 22;
    sbPos(1) = 1;               sbPos(2) = 2;
    sbPos(3) = figPos(3)-2;     sbPos(4) = H-1;
    h = axes('Parent',figHandle,'Box','off','Visible','off','Tag','sbAxes','Units','Pixels',...
        'Position',sbPos,'XLim',[0 sbPos(3)],'YLim',[0 H-1]);
    hFieldFrame = createframe(h,[1 (figPos(3) - 1)],H);
    setappdata(figHandle,'CoordsStBar',[h hFieldFrame]);  % Save it for use in ...
    set(hFieldFrame,'Visible','on')
    set(h,'HandleVisibility','off')
end
%------------------------------------

% Restore the units
drawnow;  % necessary to work around HG bug   -SLE
set(figHandle, 'Units', figUnits);
set(axHandle, 'Units', axUnits);
set(0, 'Units', rootUnits);

if ~strcmp(opt,'screen_capture'),   pixval_stsbar(figHandle);  end

% Set the nextplot property of the figure so that the
% axes object gets deleted and replaced for the next plot.
% That way, the new plot gets drawn in the default position.

% ISTO DA MERDA
%set(figHandle, 'NextPlot', 'replacechildren');
%set(axHandle, 'NextPlot', 'replacechildren');

% Warn if the display is not truesize (suspended)
% if (scale < 100)
%     try,    led = getappdata(gcf,'Semaforo');   % Don't know why but the message below prevented led -> green
%             set(led,'CData',cat(3,zeros(12),ones(12),zeros(12)));   pause(.01)
%             setappdata(gcf,'Semaforo',led);
%     end
%     %if isempty(croped)  % It's false giving this message with croped images (they were previously expanded)
%         %message = ['Image is too big to fit on screen;', sprintf('displaying at %d%% scale.', floor(scale))];
%         %msgbox(message,'Warning')
%     %end
% end

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

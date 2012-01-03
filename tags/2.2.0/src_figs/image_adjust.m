function varargout = image_adjust(varargin)
% command line arguments to image_adjust

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

	if (isempty(varargin))		return,		end

	hObject = figure('Vis','off');
	image_adjust_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hObject,'right')

	handles.h_mirone_fig = varargin{1};
	handMir = guidata(handles.h_mirone_fig);	% Get Mirone handles
	handles.OriginalImage = handMir.origFig;	% Get Mirone original Image
	if (isempty(handMir.origFig))
		errordlg('ERROR: Mirone original image is not stored in memory. Increase "Grid max size".','Error')
		delete(hObject);    return
	end
	handles.hMironeAxes = findobj(handles.h_mirone_fig,'Type','axes');
	handles.h_mirone_img = findobj(handles.h_mirone_fig,'Type','image');
	img = get(handles.h_mirone_img,'CData');
	if (isempty(img))
		errordlg('C''mon load a file first. It''s logic, isn''t it?','Error')
		delete(hObject);    return
	elseif (ndims(img) > 2)
		errordlg('ERROR: This tool doesn''t work with RGB images.','Error')
		delete(hObject);    return
	end
	%set(handles.h_mirone_img,'CDataMapping','scaled')
	handles.OriginalImage_low  = double(min(min(img))) / 255;  % E SE NAO UINT8?
	handles.OriginalImage_high = double(max(max(img))) / 255;
	handles.frameAx1_pos = get(handles.frame_ax1,'pos');
	handles.frameAx2_pos = get(handles.frame_ax2,'pos');

	% determines how to scale original image so that it fits in the defined preview area.
	% Constraining dimension in "image panel"
	[m,n] = size(img);
	frameImageConstr = min(handles.frameAx1_pos(3),handles.frameAx1_pos(4));
	dImage = ceil(sqrt(m^2 + n^2)) + 5;     %Length of Image diagonal in pixels
											% +5 added to allow for some panel border effects
	reductionFactor = (frameImageConstr/dImage);
	handles.previewImage = img_fun('imresize',img,reductionFactor);

	% Put the image in the right axes
	set(handles.figure1,'ColorMap',get(handles.h_mirone_fig,'ColorMap'))
	[m,n] = size(handles.previewImage);
	bg_aspect = m / n;
	handles.hPreviewImage = image(handles.previewImage,'Parent',handles.axes1);
	set(handles.axes1,'PlotBoxAspectRatio',[1 bg_aspect 1],'Visible','off')
	set(handles.axes1,'YDir',get(handles.hMironeAxes,'Ydir'))
	handles.hAdjustedImage = image(handles.previewImage,'Parent',handles.axes2);
	set(handles.axes2,'PlotBoxAspectRatio',[1 bg_aspect 1],'Visible','off')
	set(handles.axes2,'YDir',get(handles.hMironeAxes,'Ydir'))

	% Compute and plot the histograms
	% axes(handles.axes3);        localImhist(img,handles);
	% axes(handles.axes4);        localImhist(img,handles);
	set(hObject,'CurrentAxes',handles.axes3);        localImhist(img,handles);
	set(hObject,'CurrentAxes',handles.axes4);        localImhist(img,handles);

	% Line object for intensity transform of Histogram Eq Operation.
	Std.Interruptible = 'off';
	Std.BusyAction = 'queue';
	handles.hT = line(Std, 'Xdata',0:(1/255):1,'YData',0:(1/255):1,...
   'Parent',handles.axes5, 'Visible','off', 'Color', [0 0 1]);

	% Shit the axes so it thinks the display range is [0 255]
	set(handles.axes5,'XTick',[0 100 200]/255,'YTick',[0 100 200]/255)
	set(handles.axes5,'XTickLabel',['  0'; '100'; '200'],'YTickLabel',['  0'; '100'; '200'])

	handles.Computer = computer;

	% Line objects for Intensity Adjustment
	CtlEraseMode = 'xor';
	handles.hAdjLineCtl = line(Std, 'Tag', 'linectl', 'Parent', handles.axes5, ...
		'ButtonDownFcn', {@BtnDown,handles,1}, 'Xdata',[0 .15 .85 1], ...
		'Ydata', [.15 .15 .85 .85], 'EraseMode', CtlEraseMode, 'Color', [0 0 1]);
	handles.hAdjTopCtl = line(Std, 'Tag','topctl', 'Parent', handles.axes5, ...
		'ButtonDownFcn', {@BtnDown,handles,2}, 'LineStyle', 'none', ...
		'Xdata',.85, 'Ydata', .85, 'Marker', 'square', 'MarkerFaceCol', [.8 0 0], ...
		'EraseMode', CtlEraseMode, 'Color', [0 0 0]);
	handles.hAdjGammaCtl = line(Std, 'Tag', 'gammactl', 'Parent', handles.axes5, ...
		'ButtonDownFcn', {@BtnDown,handles,3}, 'LineStyle', 'none', ...
		'Xdata',.5, 'Ydata', .5, 'Marker', 'o', 'MarkerFaceCol', [1 1 0], ...
		'EraseMode', CtlEraseMode, 'Color', [0 0 0]);
	handles.hAdjBotCtl = line(Std, 'Tag', 'botctl', 'Parent', handles.axes5, ...
		'ButtonDownFcn', {@BtnDown,handles,4}, 'LineStyle', 'none', ...
		'Xdata',.15, 'Ydata', .15, 'Marker', 'square', 'MarkerFaceCol', [0 1 0], ...
		'EraseMode', CtlEraseMode, 'Color', [0 0 0]);

	handles = InitializeAdjustmentTool(handles);

	%------------ Give a Pro look (3D) to the frame boxes  -------------------------------
	new_frame3D(hObject, NaN)
	%------------- END Pro look (3D) -----------------------------------------------------

	guidata(hObject, handles);
	set(hObject,'Visible','on');
	if (nargout),	varargout{1} = hObject;		end

% ----------------------------------------------------------------------------
function popup_operations_CB(hObject, handles)
	UpdateOperations(handles)

% ----------------------------------------------------------------------------
function edit_gamma_CB(hObject, handles)
	gamma = str2double(get(hObject, 'String'));
	if (isempty(gamma) || gamma(1) < 0)
		set(hObject,'String','1')
		handles.Gamma = 1;
		guidata(handles.figure1,handles)
	else   
		handles.Gamma = gamma(1);
		guidata(handles.figure1,handles)
		BtnMotion([],[],handles,5);  % Redraw
		DoAdjust(handles);
	end

% ----------------------------------------------------------------------------
function push_Brighten_CB(hObject, handles)
	low = handles.LowHiBotTop(1);       high = handles.LowHiBotTop(2);
	bot = handles.LowHiBotTop(3);       top = handles.LowHiBotTop(4);
	change = .1;
	% Don't shift the whole line out of the axes range:
	if top >= bot
		if (bot+change > 1),     return;     end     % Don't make it any brighter
	else
		if (top+change > 1),     return;     end     % Don't make it any brighter
	end   
	top = top + change;     bot = bot + change;
	BtnMotion([],[],handles,5, low, high, bot, top);      % Redraw
	handles.LowHiBotTop = [low high bot top];
	guidata(handles.figure1,handles)
	DoAdjust(handles);

% ----------------------------------------------------------------------------
function push_Darken_CB(hObject, handles)
	low = handles.LowHiBotTop(1);        high = handles.LowHiBotTop(2);
	bot = handles.LowHiBotTop(3);        top = handles.LowHiBotTop(4);
	change = .1;
	% Don't shift the whole line out of the axes range:
	if (top >= bot)
		if (top-change < 0),     return;     end     % Don't make it any darker
	else
		if (bot-change < 0),     return;     end     % Don't make it any darker
	end
	top = top - change;     bot = bot - change;
	BtnMotion([],[],handles,5, low, high, bot, top);      % Redraw
	handles.LowHiBotTop = [low high bot top];
	guidata(handles.figure1,handles)
	DoAdjust(handles);

% ----------------------------------------------------------------------------
function push_incrContrast_CB(hObject, handles)
	high = get(handles.hAdjTopCtl, 'Xdata');    low = get(handles.hAdjBotCtl, 'Xdata');
	top = get(handles.hAdjTopCtl, 'Ydata');     bot = get(handles.hAdjBotCtl, 'Ydata');
	change = .1*(high-low);
	low = low + change;         high = high - change;
	BtnMotion([],[],handles,5, low, high, bot, top);      % Redraw
	DoAdjust(handles);

% ----------------------------------------------------------------------------
function push_decrContrast_CB(hObject, handles)
	high = get(handles.hAdjTopCtl, 'Xdata');    low = get(handles.hAdjBotCtl, 'Xdata');
	top = get(handles.hAdjTopCtl, 'Ydata');     bot = get(handles.hAdjBotCtl, 'Ydata');
	change = .1*(high-low);
	low = max(low - change,0);  high = min(high + change,1);
	BtnMotion([],[],handles,5, low, high, bot, top);      % Redraw
	DoAdjust(handles);

% ----------------------------------------------------------------------------
function push_decrGamma_CB(hObject, handles)
	handles.Gamma = handles.Gamma * 0.8;
	set(handles.edit_gamma, 'String', num2str(handles.Gamma)); %drawnow
	guidata(handles.figure1,handles)
	BtnMotion([],[],handles,5);  % Redraw
	DoAdjust(handles);

% ----------------------------------------------------------------------------
function push_incrGamma_CB(hObject, handles)
	handles.Gamma = handles.Gamma * 1.2;
	set(handles.edit_gamma, 'String', num2str(handles.Gamma)); %drawnow
	guidata(handles.figure1,handles)
	BtnMotion([],[],handles,5);  % Redraw
	DoAdjust(handles);

% --------------------------------------------------------------------------
function [y,x] = localImhist(varargin)

[a, n, isScaled, top] = parse_inputs(varargin{1:end-1});
if islogical(a)
    y(2) = sum(a(:));
    y(1) = numel(a) - y(2);
    y = y';
else
    y = imhistc(a, n, isScaled, top);       % Call MEX file to do work.
end

classin = class(a);
switch classin
    case 'uint8',   x = linspace(0,255,n)';
    case 'uint16',  x = linspace(0,65535,n)';
    case 'double',  x = linspace(0,1,n)';
    case 'logical', x = [0,1]';
otherwise
    errordlg('The input image must be uint8, uint16, double, or logical.','ERROR');
    return
end

if (nargout == 0)
    hAxes = gca;
    localStem(x,y)
    set(hAxes,'XLim',[min(x) max(x)],'YLim',[min(y) max(y)], 'Xtick',[], 'Ytick',[]);
end

% --------------------------------------------------------------------------
function [a, n, isScaled, top] = parse_inputs(varargin)

a = varargin{1};    n = 256;

if (isa(a,'double'))
    isScaled = 1;   top = 1;
elseif (isa(a,'uint8'))
    isScaled = 1;   top = 255;
elseif (islogical(a))
    n = 2;      isScaled = 1;    top = 1;
else % a must be uint16
    isScaled = 1;   top = 65535;
end

if (nargin ==2)
    if (numel(varargin{2}) == 1)        % IMHIST(I, N)
        n = varargin{2};
    elseif (size(varargin{2},2) == 3)   % IMHIST(X,MAP) or invalid second argument
        n = size(varargin{2},1);
        isScaled = 0;
        top = n;
    else
        messageId = sprintf('Images:%s:invalidSecondArgument', mfilename);
        error(messageId, '%s', 'Second argument must be a colormap or a positive integer.'); 
    end  
end

% --------------------------------------------------------------------------
function localStem(x,y)
%   This is a hacked code from the stem function
%   Only strictly necessary code was kept

% Set up data using fancing indexing
[m,n] = size(x);
xx = zeros(3*m,n);      xx(1:3:3*m,:) = x;
xx(2:3:3*m,:) = x;      xx(3:3:3*m,:) = NaN;

[m,n] = size(y);        yy = zeros(3*m,n);
yy(2:3:3*m,:) = y;      yy(3:3:3*m,:) = NaN;

plot(xx,yy);

% --------------------------------------------------------------------------
function handles = InitializeAdjustmentTool(handles)

high = handles.OriginalImage_high;
low = handles.OriginalImage_low;
handles.Gamma = 1.0;
set(handles.edit_gamma, 'String', 1.0);
set(handles.hAdjLineCtl,'Xdata',[0 low high 1],'Ydata',[0 0 1 1]);
set(handles.hAdjTopCtl,'Xdata',high, 'Ydata', 1);
set(handles.hAdjBotCtl,'Xdata',low, 'Ydata', 0);
set(handles.hAdjGammaCtl,'Xdata',(high+low)/2, 'Ydata', 0.5);
handles.LowHiBotTop = [low high 0 1];  % Actual Low,high,bot,top
guidata(handles.figure1,handles)

% --------------------------------------------------------------------------
function BtnDown(obj,event,handles,control)
% The input argument is the control which called the BtnDown function

handles = guidata(handles.figure1);
if control==1  % Remember out-of-limits control points
   low  = handles.LowHiBotTop(1);  
   bot  = handles.LowHiBotTop(3);
   high = handles.LowHiBotTop(2);
   top  = handles.LowHiBotTop(4);
else   % Use control points within the axes
   low  = get(handles.hAdjBotCtl, 'Xdata');
   bot  = get(handles.hAdjBotCtl, 'Ydata');
   high = get(handles.hAdjTopCtl, 'Xdata');
   top  = get(handles.hAdjTopCtl, 'Ydata');
end   

if isequal(handles.Computer(1:2), 'MA')  % For the Macintosh
   set(handles.hAdjLineCtl , 'EraseMode', 'Xor');
   set(handles.hAdjTopCtl , 'EraseMode', 'Xor');
   set(handles.hAdjBotCtl , 'EraseMode', 'Xor');
   set(handles.hAdjGammaCtl , 'EraseMode', 'Xor');
end

switch control
case 1  % The line
   pt = get(handles.axes5, 'CurrentPoint');
   handles.StartingPoint = pt(1,1:2);
   setptr(handles.figure1, 'closedhand');
case 3  % Gamma
   pt = get(handles.axes5, 'CurrentPoint');
   handles.StartingPoint = pt(1,1:2);
   setptr(handles.figure1, 'uddrag');
case {2,4}  % The Top and Bottom Controls
   setptr(handles.figure1, 'fleur');
end

handles.StartingLowHiBotTop = [low high bot top];
guidata(handles.figure1,handles)

motionfcn = {@BtnMotion,handles,control};
set(handles.figure1, 'WindowButtonMotionFcn',motionfcn,'WindowButtonUpFcn', {@BtnUp,handles});

% --------------------------------------------------------------------------
function BtnUp(obj,event,handles)

set(handles.figure1, 'WindowButtonMotionFcn','','WindowButtonUpFcn','');
DoAdjust(handles);
if isequal(handles.Computer(1:2), 'MA')  % For the Macintosh
   set(handles.hAdjLineCtl, 'EraseMode', 'normal');
   set(handles.hAdjTopCtl, 'EraseMode', 'normal');
   set(handles.hAdjBotCtl, 'EraseMode', 'normal');
   set(handles.hAdjGammaCtl, 'EraseMode', 'normal');
end
setptr(handles.figure1, 'arrow');

% --------------------------------------------------------------------------
function BtnMotion(obj,event,handles,control, low, high, bot, top)
% The first input argument is the control which called the BtnDown function
% If there are two input arguments, we are triggering a redraw from a callback
% and the second argument contains [low high bot top].

handles = guidata(handles.figure1);
% Parameters for imadjust
pt = get(handles.axes5, 'CurrentPoint');
pt = pt(1,1:2);

if (control == 5)  % Just do a redraw from a button callback function
   if (nargin == 4)
      low = handles.LowHiBotTop(1);
      high = handles.LowHiBotTop(2);
      bot = handles.LowHiBotTop(3);
      top = handles.LowHiBotTop(4);
   end   
   handles.StartingLowHiBotTop = [low high bot top];
else
   low = handles.StartingLowHiBotTop(1);
   high = handles.StartingLowHiBotTop(2);
   bot = handles.StartingLowHiBotTop(3);
   top = handles.StartingLowHiBotTop(4);
end

pt = Constrain(control, pt, handles);

% Change parameters for the line redraw
switch control
case 1                 % the user grabbed the line
   start = handles.StartingPoint;
   dx = pt(1)-start(1);     dy = pt(2)-start(2);
   low = low+dx;            high = high+dx;
   bot = bot+dy;            top = top+dy;
case 2                 % the user grabbed the top box
   high = pt(1);            top = pt(2);
case 3                 % The user grabbed the Gamma control (circle)
   % Figure out Gamma for redraw
   x = pt(1);                y = pt(2);
   gamma_old = handles.Gamma;         % Save this in case we get negative Gamma
   if ((y-bot) ~= 0) && ((x-low) ~= 0) 
      denominator = log(abs((x-low)/(high-low)));
      if denominator ~= 0
         handles.Gamma = log(abs((y-bot)/(top-bot)))/denominator;
      end
   end  
   if (handles.Gamma < 0),  handles.Gamma = gamma_old;  end
   set(handles.edit_gamma, 'String', num2str(handles.Gamma));
case 4                 % The user grabbed the Bottom control
    low = pt(1);    bot = pt(2);
end

handles.LowHiBotTop = [low high bot top];
guidata(handles.figure1,handles)

if (control == 1 || control == 5)
   if top > 1  % We dragged it off the axes
      high = (1-bot)*(high-low)/(top-bot)+low;
      top = 1;
   elseif top < 0
      high = low-bot*(high-low)/(top-bot);
      top = 0;
   elseif bot < 0
      low = high - top*(high-low)/(top-bot);
      bot = 0;
   elseif bot > 1
      low = high+(top-1)*(high-low)/(bot-top);
      bot = 1;
   end
end

% Now just draw the lines
if (abs(handles.Gamma-1) < .01)      % if Gamma is close to 1, don't worry
   set(handles.hAdjLineCtl,'Xdata',[0 low high 1], 'Ydata',[bot bot top top]);
   set(handles.hAdjBotCtl,'Xdata',low, 'Ydata', bot);
   set(handles.hAdjTopCtl,'Xdata',high, 'Ydata', top);
   set(handles.hAdjGammaCtl,'Xdata',(high+low)/2, 'Ydata', (bot+top)/2);
else
   axpos = get(handles.axes5, 'Position');
   xres = 2/axpos(3);
   curveX = low:xres:high;  % MARKED LINE
   if (length(curveX) == 1),    curveX = [low high];   end
   curveY = ( (curveX-low)/(high-low) ) .^ handles.Gamma;
   curveY = curveY*(top-bot) + bot;
   curveY(end) = top;   % This is needed because of the line marked above - HG Bug?
   set(handles.hAdjLineCtl, 'Xdata',[0 (curveX) 1], 'Ydata',[bot (curveY) top]);
   set(handles.hAdjBotCtl,'Xdata',low,'Ydata',bot);
   set(handles.hAdjTopCtl,'Xdata',high,'Ydata', top);
   npts = length(curveX);
   if (rem(npts,2) == 0)
      set(handles.hAdjGammaCtl,...
         'Xdata', (curveX(npts/2) + curveX(1+npts/2))/2, ...
         'Ydata', (curveY(npts/2) + curveY(1+npts/2))/2);
   else
      set(handles.hAdjGammaCtl, 'Xdata', curveX((npts+1)/2),'Ydata', curveY((npts+1)/2));
   end
end

% --------------------------------------------------------------------------
function out = Constrain(control, in, handles)
% Make sure the following conditions are all met:
% 1. The line still represents a function (1 to 1 mapping)
% 2. No part of the line extends beyond the range of the axes

%handles = guidata(handles.figure1);
low = handles.StartingLowHiBotTop(1);
high = handles.StartingLowHiBotTop(2);
bot = handles.StartingLowHiBotTop(3);
top = handles.StartingLowHiBotTop(4);

out = in;

out(in > 1) = 1;  % Make sure the point stays in the axes
out(in < 0) = 0;

axpos = get(handles.axes5, 'Position');
xres = 1/axpos(3);
yres = 1/axpos(4);
switch control
case 1                 % the user grabbed the line
   start = handles.StartingPoint;
   dx = out(1)-start(1);
   if high+dx > 1
      out(1) = 1-high+start(1);
   end
   if low+dx < 0
      out(1) = start(1)-low;
   end
case 2                 % the user grabbed the top box
   if out(1) <= low,   
      out(1) = low + xres;  % One pixel
   end
case 3                 % The user grabbed the Gamma control (circle)
   start = handles.StartingPoint;
   out(1) = start(1);  % Only move vertically
   if top > bot  
      if out(2) >= top
         out(2) = top-yres;  % One pixel
      elseif out(2) <= bot
         out(2) = bot+yres;
      end
   elseif top <= bot
      if out(2) > bot
         out(2) = bot-yres;
      elseif out(2) <= top
         out(2) = top+yres;
      end
   end
case 4                 % The user grabbed the Bottom control
   if out(1) >= high,  
      out(1) = high - xres; 
   end
end

% --------------------------------------------------------------------------
function DoAdjust(handles)

	setptr(handles.figure1, 'watch');

	% The add/subtracts in Constrain can introduce eps level errors which put us outside of [0 1]
	low  = max(0.0, get(handles.hAdjBotCtl, 'Xdata'));
	bot  = max(0.0, get(handles.hAdjBotCtl, 'Ydata'));
	high = min(1.0, get(handles.hAdjTopCtl, 'Xdata'));
	top  = min(1.0, get(handles.hAdjTopCtl, 'Ydata'));

	if( abs(high-low)<eps )  % Protect imadjust against divide by 0
	   high = low+.0001;
	end

	% Adjust first the preview image
	img = get(handles.hPreviewImage, 'Cdata');
	imgAd = img_fun('imadjust_j',img, [low;high],[bot;top],handles.Gamma);
	set(handles.hAdjustedImage, 'Cdata', imgAd);

	% Adjust the original image
	imgAd = img_fun('imadjust_j',handles.OriginalImage, [low;high],[bot;top],handles.Gamma);
	set(handles.h_mirone_img, 'Cdata', imgAd);

	axes(handles.axes4);        localImhist(imgAd,handles);
	setptr(handles.figure1, 'arrow');

% --------------------------------------------------------------------------
function UpdateOperations(handles)

op = get(handles.popup_operations, 'Value');
switch op
	case 1   % Intensity Adjustment
       set(handles.hT, 'Visible', 'off');   % Turn off histeq part
       % Turn on imadjust part
       set([handles.push_Brighten handles.push_Darken handles.push_incrContrast handles.push_decrContrast handles.push_incrGamma ...
             handles.push_decrGamma handles.text_gammaLabel handles.edit_gamma], 'Enable', 'on');
       set([handles.hAdjLineCtl handles.hAdjTopCtl handles.hAdjGammaCtl handles.hAdjBotCtl], 'Visible', 'on');
       InitializeAdjustmentTool(handles);
       DoAdjust(handles);
	case 2   % Histogram Equalization
       % Turn off imadjust part
       set([handles.push_Brighten handles.push_Darken handles.push_incrContrast handles.push_decrContrast handles.push_incrGamma ...
             handles.push_decrGamma handles.text_gammaLabel handles.edit_gamma],'Enable', 'off');
       set([handles.hAdjLineCtl handles.hAdjTopCtl handles.hAdjGammaCtl handles.hAdjBotCtl], 'Visible', 'off');
       set(handles.hT, 'Visible', 'on');    % Turn on histeq part
       EqualizeHistogram(handles);
	case 3   % Adaptative Histogram Equalization
       % Turn off imadjust part
       set([handles.push_Brighten handles.push_Darken handles.push_incrContrast handles.push_decrContrast handles.push_incrGamma ...
             handles.push_decrGamma handles.text_gammaLabel handles.edit_gamma],'Enable', 'off');
       set([handles.hAdjLineCtl handles.hAdjTopCtl handles.hAdjGammaCtl handles.hAdjBotCtl], 'Visible', 'off');
       set(handles.hT, 'Visible', 'on');    % Turn on histeq part
       AdaptiveHisteq(handles)
	otherwise
       error('Undefined Operation');
	end

% --------------------------------------------------------------------------
function EqualizeHistogram(handles)
	setptr(handles.figure1, 'watch');

	% Equalize first the preview image
	img = get(handles.hPreviewImage, 'CData');
	imgEq = img_fun('histeq_j',img,256);
	set(handles.hAdjustedImage, 'CData', imgEq)

	% Now the original image
	[imgEq,T] = img_fun('histeq_j',handles.OriginalImage,256);
	set(handles.h_mirone_img, 'Cdata', imgEq);

	set(handles.hT, 'Xdata',0:(1/255):1,'YData',T);
	axes(handles.axes4);    localImhist(imgEq,handles);
	setptr(handles.figure1, 'arrow');

% --------------------------------------------------------------------------
function AdaptiveHisteq(handles)
setptr(handles.figure1, 'watch');

% Adjust first the preview image
img = get(handles.hPreviewImage, 'CData');
imgEq = img_fun('adapthisteq',img,'NumTiles',[8 8],'ClipLimit',0.02,'Distribution','rayleigh');
set(handles.hAdjustedImage, 'CData', imgEq)

% Now the original image
imgEq = img_fun('adapthisteq',handles.OriginalImage,'NumTiles',[8 8],'ClipLimit',0.02,'Distribution','rayleigh');
set(handles.h_mirone_img, 'Cdata', imgEq);

axes(handles.axes4);    localImhist(imgEq,handles);
setptr(handles.figure1, 'arrow');


% --- Creates and returns a handle to the GUI figure. 
function image_adjust_LayoutFcn(h1)

set(h1,...
'PaperUnits',get(0,'defaultfigurePaperUnits'),...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'MenuBar','none',...
'Name','Image Adjust',...
'NumberTitle','off',...
'Position',[520 439 530 361],...
'Resize','off',...
'Tag','figure1');

axes('Parent',h1,'Units','pixels','Position',[10 200 161 151],'Tag','axes1');
axes('Parent',h1,'Units','pixels','Position',[180 200 161 151],'Tag','axes2');
axes('Parent',h1,'Units','pixels','Position',[10 20 151 141],'Tag','axes3');
axes('Parent',h1,'Units','pixels','Position',[180 20 151 141],'Tag','axes4');

uicontrol('Parent',h1,'FontName','Helvetica','FontSize',9,...
'Position',[370 344 144 16],...
'String','Output vs. Input Intensity',...
'Style','text','Tag','text5');

axes('Parent',h1,'Units','pixels',...
'Position',[370 200 151 141],'Tag','axes5');

uicontrol('Parent',h1,'Position',[340 20 181 141],'Style','frame','Tag','frame1');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@main_uiCB,...
'Position',[350 128 161 22],...
'String',{'Intensity adjustment'; 'Histogram Equalization'; 'Adaptative Histogram Eq'},...
'Style','popupmenu','Value',1,'Tag','popup_operations');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@main_uiCB,...
'HorizontalAlignment','right',...
'Position',[451 165 49 21],...
'String','1','Style','edit','Tag','edit_gamma');

uicontrol('Parent',h1,...
'Call',@main_uiCB,...
'Position',[350 88 75 19],...
'String','+ Brightness','Tag','push_Brighten');

uicontrol('Parent',h1,...
'Call',@main_uiCB,...
'Position',[437 88 75 19],...
'String','- Brightness','Tag','push_Darken');

uicontrol('Parent',h1,...
'Call',@main_uiCB,...
'Position',[350 61 75 19],...
'String','+ Contrast','Tag','push_incrContrast');

uicontrol('Parent',h1,...
'Call',@main_uiCB,...
'Position',[437 61 75 19],...
'String','- Contrast','Tag','push_decrContrast');

uicontrol('Parent',h1,...
'Call',@main_uiCB,...
'Position',[350 31 75 19],...
'String','+ Gamma','Tag','push_incrGamma');

uicontrol('Parent',h1,...
'Call',@main_uiCB,...
'Position',[437 31 75 19],...
'String','- Gamma','Tag','push_decrGamma');

uicontrol('Parent',h1,'FontName','Helvetica','FontSize',9,...
'Position',[220 166 71 15],'String','Histogram',...
'Style','text','Tag','text2');

uicontrol('Parent',h1,'FontName','Helvetica','FontSize',9,...
'Position',[50 166 71 15],'String','Histogram',...
'Style','text','Tag','text3');

uicontrol('Parent',h1,'Position',[403 168 45 15],...
'String','Gamma:','Style','text','Tag','text_gammaLabel');

uicontrol('Parent',h1,'Position',[10 200 161 151],...
'Style','frame','Tag','frame_ax1');

uicontrol('Parent',h1,'Position',[180 200 161 151],...
'Style','frame','Tag','frame_ax2');

function main_uiCB(hObject, eventdata)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));

function imscroll_j(hAxes,opt)
% IMSCROLL_J updates the slider values of a figure when called from the
% zoom or pan functions. It contains also functions to update the image
% when the slider are changed.
%
% The sliders MUST have tags with HOR and VER, respectively for horizontal
% and vertical sliders. Both or only one of the sliders may exist.
%
% The sliders MUST be registered on the axes where they work on via appdata.
% EG: setappdata(hAxes,'SliderAxes',[hSlider_hor hSlider_ver])
%
% An example of the Slider settings is given bellow
% set(hSlider_hor,'Min',0,'Max',1,'Value',0,'Tag','HOR','Callback',{@slider_Cb,hAxes,'SetSliderHor'})
% set(hSlider_ver,'Min',0,'Max',1,'Value',0,'Tag','VER','Callback',{@slider_Cb,hAxes,'SetSliderVer'})
%
% The "slider_Cb" Callback is a simple function with this form
% function slider_Cb(obj,evt,hAxes,opt)
% imscroll_j(ax,opt)
%
% See the SCROLL_DEMO function for a working example of using IMSCROLL_J
%
% NOTE1: the zoom_j.m function is a R13 modified version that interfaces with
% the imscroll_j fucntion. It may be used by R13 or R14 releases.
%
% NOTE2: The PAN function was introduced only on R14 and in order that it knows
% about IMSCROLL_J you must do the following change to the pan.m function
%
% pan.m         (R14 only)
% On the locWindowButtonMotionFcn function at the end of the 
% "if ok2pan" test (but inside this IF test)  ad
% imscroll_j(ax,'PanSetSliders')

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

	hSliders = getappdata(hAxes,'SliderAxes');
	if (isempty(hSliders)),     return;     end

	switch opt
		case 'PanSetSliders'
			PanZoomSetSliders(hAxes,hSliders)
		case 'ZoomSetSliders'
			PanZoomSetSliders(hAxes,hSliders)
		case 'SetSliderHor'
			if (numel(hSliders) == 2)
				hSliders = findobj(hSliders,'Tag','HOR');
			end
			SetSliderHor(hAxes,hSliders)
		case 'SetSliderVer'
			if (numel(hSliders) == 2)
				hSliders = findobj(hSliders,'Tag','VER');
			end
			SetSliderVer(hAxes,hSliders)
	end

%------------------------------------------------------------------------------------
function PanZoomSetSliders(ax,hSliders)
% Find out which slider function to call

	limits = getappdata(ax,'ThisImageLims');
	if (isempty(limits)),     limits = get_full_limits(ax);   end

	x = limits(1:2);                  y = limits(3:4);
	xVisibleLim = get(ax,'XLim');   yVisibleLim = get(ax,'YLim');

	for (i = 1:numel(hSliders))
		ud = get(hSliders(i),'Tag');           % NEED TO TEST IF THIS MAY GIVE ERROR
		if (strcmpi(ud(1),'h'))
			if (abs(diff([x; xVisibleLim])) < 1e-5)
				set(hSliders(i),'Visible','off')
			else
				set(hSliders(i),'Visible','on')
				ZoomSetSliderHOR(hSliders(i),x,xVisibleLim)
			end
		elseif (strcmpi(ud(1),'v'))
			if (abs(diff([y; yVisibleLim])) < 1e-5)
				set(hSliders(i),'Visible','off')
			else
				set(hSliders(i),'Visible','on')
				ZoomSetSliderVER(ax,hSliders(i),y,yVisibleLim)
			end
		end
	end

%------------------------------------------------------------------------------------
function ZoomSetSliderHOR(hSlider,x,xVisibleLim)
% Set the Horizontal slider value from a zoom command

	VisWinWidth = diff(xVisibleLim);		% Width of the visible window
	xVisWinLast = xVisibleLim(1) - x(1);	% X coords of the last visible point

	f = VisWinWidth/diff(x);
	maxStep = 1/(1/f - 1 + eps);			% +eps to avoid Divide by zero warnings
	maxStep = max(maxStep,0);
	minStep = min(1, maxStep/10);			% must be between 0 and 1.
	if (~isequal([minStep maxStep],[0 0]))
		set(hSlider,'SliderStep',[minStep maxStep],'Visible','on')
	else
		set(hSlider,'Visible','off')		% Nothing more to unzoom, so we may hide the slider
		return
	end

	to_scale = VisWinWidth / diff(x);		% Factor to scale the xSlider value in the [0 1] interval
	xSlider =  xVisWinLast / diff(x);		% somewhere in the interval [0 <1]
	xSlider = xSlider / (1 - to_scale);		% Now it should range from 0 to 1
	xSlider = min(max(xSlider,0),1);		% Make sure it is so.
	set(hSlider,'Value',xSlider),	refresh

%------------------------------------------------------------------------------------
function ZoomSetSliderVER(ax,hSlider,y,yVisibleLim)
% Set the Vertical slider value from a zoom command

	VisWinHeight = diff(yVisibleLim);		% Height of the visible window
	yVisWinLast = yVisibleLim(2) - y(1);	% Y relative coords of the last visible point

	f = VisWinHeight/diff(y);
	maxStep = 1/(1/f - 1 + eps);
	maxStep = max(maxStep,0);
	minStep = min(1, maxStep/10);			% must be between 0 and 1.
	if (~isequal([minStep maxStep],[0 0]))
		set(hSlider,'SliderStep',[minStep maxStep],'Visible','on')
	else
		set(hSlider,'Visible','off')		% Nothing more to unzoom, so we may hide the slider
		return
	end

	to_scale = VisWinHeight / diff(y);		% Factor to scale the ySlider value in the [0 1] interval
	ySlider =  yVisWinLast / diff(y);				% somewhere in the interval [>0 1]
	ySlider = (ySlider - to_scale) / (1 - to_scale);% Now it should range from 0 to 1
	ySlider = min(max(ySlider,0),1);				% Make sure it is so.
	if (strcmp(get(ax,'YDir'),'reverse'))
		ySlider = 1 - ySlider;
	end
	set(hSlider,'Value',ySlider),	refresh

%------------------------------------------------------------------------------------
function SetSliderHor(ax,hSlider)
% Update the image according to de Horizontal slider value

	limits = getappdata(ax,'ThisImageLims');
	if (isempty(limits)),     limits = get_full_limits(ax);   end

	x = limits(1:2);
	xVisibleLim = get(ax,'XLim');
	VisWinWidth = diff(xVisibleLim);				% Width of the visible window
	%xVisWinLast = xVisibleLim(1) + VisWinWidth;	% X coords of the last visible point

	val = get(hSlider,'Value');
	to_scale = VisWinWidth / diff(x);				% Factor to scale the xSlider value in the [0 <1] interval
	val = val * (1-to_scale);						% somewhere in the interval [0 <1]
	new_xLeft = x(1) + val*diff(x);
	new_xLim = [new_xLeft new_xLeft+VisWinWidth];
	if (new_xLim(2) > x(2))
		new_xLim = [x(2)-VisWinWidth x(2)];
	end
	set(ax,'XLim',new_xLim);   refresh;

%------------------------------------------------------------------------------------
function SetSliderVer(ax,hSlider)
% Update the image according to de Vertical slider value

	limits = getappdata(ax,'ThisImageLims');
	if (isempty(limits)),     limits = get_full_limits(ax);   end

	y = limits(3:4);
	yVisibleLim = get(ax,'YLim');
	VisWinHeight = diff(yVisibleLim);				% Height of the visible window
	%yVisWinLast = yVisibleLim(1) + VisWinHeight;	% Y coords of the last visible point

	val = 1 - get(hSlider,'Value');
	to_scale = VisWinHeight / diff(y);				% Factor to scale the ySlider value in the [0 <1] interval

	if (strcmp(get(ax,'YDir'),'normal'))
		val = val * (1-to_scale);					% somewhere in the interval [0 <1]
		new_yTop = y(2) - val*diff(y);
		new_yLim = [new_yTop-VisWinHeight new_yTop];
		if (new_yLim(1) < y(1))
			new_yLim = [y(1) y(1)+VisWinHeight];
		end
	else
		val = val * (1-to_scale);					% somewhere in the interval [0 <1]
		new_yTop = y(1) + val*diff(y);
		new_yLim = [new_yTop new_yTop+VisWinHeight];
		if (new_yLim(1) > y(2))
			new_yLim = [y(2)-VisWinHeight y(2)];
		end
	end

	set(ax,'YLim',new_yLim),	refresh

%------------------------------------------------------------------------------------
function limits = get_full_limits(hAx)
%Get the full limits of the contents of axes hAx.
%   If the axes has no children or they do not contribute to the limits, 
%	the limits of the axes hAx are returned.   

	xmin = nan;		xmax = nan;
	ymin = nan;		ymax = nan;

	hObj = findobj(hAx,'-depth',1);		% Get childrens

	for (i = 1:numel(hObj))
		type = get(hObj(i),'type');
		switch type
			case 'axes'
				x = get(hObj(i), 'xlim');
				y = get(hObj(i), 'ylim');
				xmin = x(1);	xmax = x(2);
				ymin = y(1);	ymax = y(2);
			case {'image', 'line', 'surface'}
				if (strcmp(type, 'image'))
					% Make sure X|Ydata are in pixel registered mode. If they were already than we have
					% a problem, but that shouldn't happen because 'ThisImageLims' has been filled before
					x = get(hObj(i), 'XData');
					y = get(hObj(i), 'YData');
					[m, n, k] = size(get(hObj(i),'CData'));
					if (numel(x) == 2)		% An image
						dx2 = ( (x(end) - x(1)) / (n - 1) ) / 2;
						dy2 = ( (y(end) - y(1)) / (m - 1) ) / 2;
					else
						dx2 = (x(2) - x(1)) / 2;
						dy2 = (y(2) - y(1)) / 2;
					end
					x = [x(1) x(end)] + [-dx2 dx2];
					y = [y(1) y(end)] + [-dy2 dy2];
				end

			case 'patch'
				verts = get(hObj(i), 'vertices');
				faces = get(hObj(i), 'faces');
				verts = verts(faces(isfinite(faces)),:);
				x = verts(:,1);
				y = verts(:,2);
			end

			xmin = min(xmin, min(x));
			xmax = max(xmax, max(x));
			ymin = min(ymin, min(y));
			ymax = max(ymax, max(y));
	end

	limits = [xmin xmax ymin ymax];

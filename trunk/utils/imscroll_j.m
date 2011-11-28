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

	lims = getappdata(ax,'ThisImageLims');
	if (isempty(lims)),     lims = get_full_limits(ax);   end

	x = lims(1:2);                  y = lims(3:4);
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

	lims = getappdata(ax,'ThisImageLims');
	if (isempty(lims)),     lims = get_full_limits(ax);   end

	x = lims(1:2);
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

	lims = getappdata(ax,'ThisImageLims');
	if (isempty(lims)),     lims = get_full_limits(ax);   end

	y = lims(3:4);
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
function lims = get_full_limits(h)
%GET_FULL_LIMITS 3D object limits.
%   LIMS=GET_FULL_LIMITS(H)  limits of the objects in vector H.
%   LIMS=GET_FULL_LIMITS(AX) limits of the objects that are children 
%                      of axes AX. If the axes has no children or none of
%                      children contribute to the limits, the limits of the
%                      axes AX are returned.    
%   LIMS=GET_FULL_LIMITS     limits of the objects that are children 
%                      of the current axes. If the axes has no children or none of
%                      children contribute to the limits, the limits of the
%                      current axes are returned.
%   GET_FULL_LIMITS calculates the 3D limits of the objects specified. The
%   limits are returned in the form [xmin xmax ymin ymax zmin zmax].
%
%  If the limits are invalid, an empty array is returned.

%   Copyright 1984-2009 The MathWorks, Inc.

	if (nargin == 0),	h = gca;	end

	h = h(ishandle(h));		% A cose des mouches
	%Save a reference to the input object
	hInput = h;

	xmin = nan;		xmax = nan;
	ymin = nan;		ymax = nan;
	zmin = nan;		zmax = nan;

	if (strcmp(get(h, 'type'), 'axes'))
		ch = [];
		for (i = 1:numel(h))		% Get the axes children.
			newCh = findobj(h(i));
			ch = [ch;newCh(2:end)]; %#ok<AGROW>
		end
		h = ch;
	end

	for (i = 1:numel(h))
		validtype = false;
		type = get(h(i),'type');
		switch type
			case {'image', 'line', 'surface'}
				xd = get(h(i), 'xdata');
				yd = get(h(i), 'ydata');
				validtype = ~isempty(xd) && ~isempty(yd) && any(isfinite(xd(:))) && any(isfinite(yd(:)));
				if (strcmp(type, 'image'))
					if ( validtype )
						zd = 0;
						[xd,yd] = localGetImageBounds(h(i));
					end
				else
					zd = get(h(i), 'zdata');
					if isempty(zd)
						if (strcmp(type, 'line')),	zd = 0; % a line can have empty zdata
						else						validtype = false;
						end
					end
				end

			case 'patch'
				v = get(h(i), 'vertices');
				validtype = ~isempty(v) && any(isfinite(v(:)));
				if validtype
					f = get(h(i), 'faces');
					v = v(f(isfinite(f)),:);
					xd = v(:,1);
					yd = v(:,2);
					if (size(v,2) == 2),	zd = 0;
					else					zd = v(:,3);
					end
				end
			end

		if validtype
			xmin = min(xmin, min(xd(:)));
			xmax = max(xmax, max(xd(:)));
			ymin = min(ymin, min(yd(:)));
			ymax = max(ymax, max(yd(:)));
			zmin = min(zmin, min(zd(:)));
			zmax = max(zmax, max(zd(:)));
		end
	end

	lims = [xmin xmax ymin ymax zmin zmax];
	if any(isnan(lims)),	lims = [];	end
	% If the lims are empty and the input argument was
	% an AXES object, return the lims of the AXES.
	% This can happen if either the axes has no children
	% or none of the children contributed to the limits.
	if( isempty(lims) && (numel(findobj(hInput,'flat','type','axes')) == numel(hInput)))
		xl = xlim(hInput);
		yl = ylim(hInput);
		zl = zlim(hInput);
		xmin = xl(1);  xmax = xl(2); 
		ymin = yl(1);  ymax = yl(2); 
		zmin = zl(1);  zmax = zl(2); 
		lims = [xmin xmax ymin ymax zmin zmax];
	end

%----------------------------------
function [xd,yd] = localGetImageBounds(h)
% Determine the bounds of the image

xdata = get(h,'XData');
ydata = get(h,'YData');
cdata = get(h,'CData');
m = size(cdata,1);
n = size(cdata,2);

[xd(1), xd(2)] = localComputeImageEdges(xdata,n);
[yd(1), yd(2)] = localComputeImageEdges(ydata,m);

%----------------------------------
function [min,max]= localComputeImageEdges(xdata,num)
% Determine the bounds of an image edge

% This algorithm is an exact duplication of the image HG c-code 
% Reference: src/hg/gs_obj/image.cpp, ComputeImageEdges(...)

offset = .5;
nreals = length(xdata);
old_nreals = nreals;

if (old_nreals>1 && isequal(xdata(1),xdata(end)))
    nreals = 1;
end

first_last(1) = 1;
first_last(2) = num;

if (num==0) 
	min = nan;
	max = nan;
else
    first_last(1) = xdata(1);
    if (nreals>1) 
        first_last(2) = xdata(end);
    else
        first_last(2) = first_last(1) + num - 1;
    end
    
    % Data should be monotonically increasing
    if (first_last(2) < first_last(1)) 
        first_last = fliplr(first_last);
    end
    
    if (num > 1) 
		offset = (first_last(2) - first_last(1)) / (2 * (num-1));
    elseif (nreals > 1)
		offset = xdata(end) - xdata(1);
    end
    min = first_last(1) - offset;
    max = first_last(2) + offset;
end

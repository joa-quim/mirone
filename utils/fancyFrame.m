function fancyFrame(handles, opt)
% Create a FANCY frame with 'frame' objects if OPT == 'set' or delete it if OPT == 'unset'

%	Copyright (c) 2004-2013 by J. Luis
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

% $Id: fancyFrame.m 7921 2016-06-01 01:42:56Z j $

	if (opt(1) == 'p')		% Temporary solution for printing (create a frame made of patches)
		frame_patch(handles, opt(2:end))
		return
	end

	xLim = get(handles.axes1, 'XLim');		yLim = get(handles.axes1, 'YLim');

	hFrancyFrames = findobj(handles.figure1, '-depth',1, 'Style','frame', 'Tag', 'FancyFrame');
	delete(hFrancyFrames);		% For now we just delete but in future we will recicle them
	if (strcmp(opt,'unset'))
		set(handles.axes1,'TickDir','in', 'XLimMode','auto','XTickMode','auto', 'YLimMode','auto','YTickMode','auto')
		set(handles.axes1, 'XLim', xLim, 'YLim', yLim);
		return
	end

	fW = 4;		% Frame width in points
	axLW = get(handles.axes1,'LineWidth');
	units = get(handles.axes1,'Units');			set(handles.axes1,'Units','points')
	pos = get(handles.axes1,'Pos');				set(handles.axes1,'Units',units)
	units = get(handles.figure1,'Units');		set(handles.figure1,'Units','points')
	posFigPts = get(handles.figure1,'Pos');		set(handles.figure1,'Units',units)
	posAxNorm = get(handles.axes1,'Position');	% Axes position in normalized units

	h = uicontrol('Parent',handles.figure1,'Style','frame','Units','points', 'Tag','FancyFrame', ...
		'Pos',[pos(1:3) 0]+[-fW -fW 2*fW fW+axLW]);		% Bot
	set(h,'Units','normalized')
	p = get(h,'Pos');	BarThickX = p(4);	% X thickness our frame bars
	h = uicontrol('Parent',handles.figure1,'Style','frame','Units','points', 'Tag','FancyFrame', ...
		'Pos',[pos(1)+pos(3)-axLW pos(2)-fW fW pos(4)+2*fW]);	% Right
	set(h,'Units','normalized')
	p = get(h,'Pos');	BarThickY = p(3);
	h = uicontrol('Parent',handles.figure1,'Style','frame','Units','points', 'Tag','FancyFrame', ...
		'Pos',[pos(1)-fW pos(2)+pos(4)-axLW pos(3)+2*fW fW]);	% Top
	set(h,'Units','normalized')
	h = uicontrol('Parent',handles.figure1,'Style','frame','Units','points', 'Tag','FancyFrame', ...
		'Pos',[pos(1)-fW+axLW pos(2)-fW fW pos(4)+2*fW]);	% Left
	set(h,'Units','normalized')

	ddx = dtick(diff(xLim), handles.geog > 0);
	ddy = dtick(diff(yLim), handles.geog > 0);

	xtick = (ddx*ceil(xLim(1)/ddx)):ddx:xLim(2);	% XTicks in Data units
	XDataPerUnit = diff(xLim) / posAxNorm(3);
	xtFU = (xtick - xLim(1)) / XDataPerUnit;		% xtick in Figure's units but still referenced to axes (0,0)
	if (xtFU(1) ~= 0),	xtFU = [0 xtFU];	end		% At this point first element must always be at the axes origin
	if (abs(xtFU(end)-posAxNorm(3)) > 1e-6), xtFU = [xtFU posAxNorm(3)];	end		% Ensure that last black bar is closed
	ytick = (ddy*ceil(yLim(1)/ddy)):ddy:yLim(2);
	YDataPerUnit = diff(yLim) / posAxNorm(4);
	ytFU = (ytick - yLim(1)) / YDataPerUnit;		% ytick in Figure's units (normalized)
	if (ytFU(1) ~= 0),	ytFU = [0 ytFU];	end		% At this point first element must always be at the axes origin
	if (abs(ytFU(end)-posAxNorm(4)) > 1e-6), ytFU = [ytFU posAxNorm(4)];	end
	xtFU = xtFU + posAxNorm(1);						% Make them refer to Figure's origin (0,0)
	ytFU = ytFU + posAxNorm(2);

	tinny = axLW / posFigPts(4);					% Axes line width in normalized units

	for (k = 1:2:numel(xtFU)-rem(numel(xtFU),2))	% Make sure loop variable k is always < numel(xtFU)
		uicontrol('Parent',handles.figure1,'Style','frame','BackgroundColor','k','Units','normalized', ...	% Bot
			'Pos',[xtFU(k) posAxNorm(2)-BarThickX+tinny xtFU(k+1)-xtFU(k) BarThickX], 'Tag','FancyFrame');
		uicontrol('Parent',handles.figure1,'Style','frame','BackgroundColor','k','Units','normalized', ...	% Top
			'Pos',[xtFU(k) posAxNorm(2)+posAxNorm(4)-tinny xtFU(k+1)-xtFU(k) BarThickX], 'Tag','FancyFrame');
	end

	for (k = 1:2:numel(ytFU)-rem(numel(ytFU),2))
		uicontrol('Parent',handles.figure1,'Style','frame','BackgroundColor','k','Units','normalized', ...	% Left
			'Pos',[xtFU(1)-BarThickY+tinny ytFU(k) BarThickY ytFU(k+1)-ytFU(k)], 'Tag','FancyFrame');
		uicontrol('Parent',handles.figure1,'Style','frame','BackgroundColor','k','Units','normalized', ...	% Right
			'Pos',[posAxNorm(1)+posAxNorm(3) ytFU(k) BarThickY ytFU(k+1)-ytFU(k)], 'Tag','FancyFrame');
	end

	% Count number of decimals and force to use the same numbel for ALL labels.
	XTickLabel = num2str(xtick(:));			% Rely on default value by num2str to decide on the number of decimals.
	ind = strfind(XTickLabel(1,:), '.');
	if (isempty(ind)),	n_dec = 0;
	else				n_dec  = size(XTickLabel, 2) - ind;
	end

	XTickLabel = num2str(xtick(:), sprintf('%%.%df', n_dec));
	
	YTickLabel = num2str(ytick(:));
	ind = strfind(YTickLabel(1,:), '.');
	if (isempty(ind)),	n_dec = 0;
	else				n_dec  = size(YTickLabel, 2) - ind;
	end
	YTickLabel = num2str(ytick(:), sprintf('%%.%df', n_dec));

	set(handles.axes1,'TickDir','out','XTick',xtick, 'XTickLabel',XTickLabel, 'YTick',ytick, 'YTickLabel',YTickLabel)
	setappdata(handles.axes1, 'XTickOrig', XTickLabel),		setappdata(handles.axes1, 'XTickOrigNum', xtick)
	setappdata(handles.axes1, 'YTickOrig', YTickLabel),		setappdata(handles.axes1, 'YTickOrigNum', ytick)

% -----------------------------------------------------------------------------------------
function frame_patch(handles, opt)
% This is solution used by François Beauducel in his 'DEM' function (BSD) but the problem is
% that the patches get zoomed and move with zooming, which is an unacceptable behavior for axes.
% But my 'frame' solution does not print, so before printing I call this function to create a
% temporary FANCY frame made of patches, which are axes children and therefore print.
% It's still far from perfect, namely the two FANCY frames are not exactly equal meaning that
% what we see on screen is not exactly what we get on printing (but not that different).
% http://www.ipgp.fr/~beaudu/matlab/dem.m

	if (strcmp(opt,'unset'))
		hFrancyFrames = findobj(handles.axes1, '-depth',1, 'Type','patch', 'Tag', 'PatchFrame');
		delete(hFrancyFrames),		return
	end

	xLim = get(handles.axes1, 'XLim');		yLim = get(handles.axes1, 'YLim');
	DAR = get(handles.axes1, 'DataAspectRatio');
	bwy = 0.009*diff(yLim);			% Y border width = 1%
	bwx = bwy/DAR(2);				% border width (in degree of longitude)

	patch([xLim(1)-bwx,xLim(2)+bwx,xLim(2)+bwx,xLim(1)-bwx],yLim(1) - bwy*[0,0,1,1],'k','FaceColor','none','clipping','off','Tag','PatchFrame')
	patch([xLim(1)-bwx,xLim(2)+bwx,xLim(2)+bwx,xLim(1)-bwx],yLim(2) + bwy*[0,0,1,1],'k','FaceColor','none','clipping','off','Tag','PatchFrame')
	patch(xLim(1) - bwx*[0,0,1,1],[yLim(1)-bwy,yLim(2)+bwy,yLim(2)+bwy,yLim(1)-bwy],'k','FaceColor','none','clipping','off','Tag','PatchFrame')
	patch(xLim(2) + bwx*[0,0,1,1],[yLim(1)-bwy,yLim(2)+bwy,yLim(2)+bwy,yLim(1)-bwy],'k','FaceColor','none','clipping','off','Tag','PatchFrame')

	ddx = dtick(diff(xLim), handles.geog > 0);
	ddy = dtick(diff(yLim), handles.geog > 0);

	xtick = get(handles.axes1,'XTick');
	for xt = xtick(1:2:end)
		dt = ddx - max(0,xt + ddx - xLim(2));
		patch(repmat(xt + dt*[0,1,1,0]',[1,2]),[yLim(1) - bwy*[0,0,1,1];yLim(2) + bwy*[0,0,1,1]]','k','clipping','off','Tag','PatchFrame')
	end

	ytick = get(handles.axes1,'YTick');
	for yt = ytick(1:2:end)
		dt = ddy - max(0,yt + ddy - yLim(2));
		patch([xLim(1) - bwx*[0,0,1,1];xLim(2) + bwx*[0,0,1,1]]',repmat(yt + dt*[0,1,1,0]',[1,2]),'k','clipping','off','Tag','PatchFrame')
	end

% -----------------------------------------------------------------------------------------
function dd = dtick(dlim, deg)
% DTICK Tick intervals
% From François Beauducel's 'DEM' function (BSD)

	if (nargin < 2),	deg = false;	end

	if (deg && dlim <= 2/60)			% less than 2 minutes: base 36
		m = 10^floor(log10(dlim*36))/36;
	elseif (deg && dlim <= 2)			% less than 2 degrees: base 6
		m = 10^floor(log10(dlim*6))/6;
	else								% more than few degrees or not degrees: decimal rules
		m = 10^floor(log10(dlim));
	end
	p = ceil(dlim/m);
	if (p <= 1),		dd = 0.1*m;
	elseif (p == 2),	dd = 0.2*m;
	elseif (p <= 5),	dd = 0.5*m;
	else				dd = m;
	end

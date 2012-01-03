function varargout = image_histo(varargin)
% Show Histogram of indexed or RGB images

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

	if (isempty(varargin)),		return,		end
	hMirHand = varargin{1};
	img = get(hMirHand.hImg, 'CData');
	if (isempty(img)),			return,		end

	hObject = figure('Vis','off');
	image_histo_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hObject,'right')

	hPolyg = [];
	if (nargin == 2)
		hPolyg = varargin{2};
	end
	handles.hMirFig = hMirHand.figure1;
	handles.is_polyg = false;					% For ROI histos
	handles.validGrid = hMirHand.validGrid;		% Grids have NaN color as first

	if (ndims(img) == 3)
		handles.isRGB = 1;
	elseif (ndims(img) == 2)
		handles.isRGB = 0;
		% OK, here we must reshape the GUI, but we can delete some uis first
		pos_fig = get(hObject,'pos');
		pos_a1 = get(handles.axes1,'pos');
		new_height = pos_fig(4) - pos_a1(2) + 40;
		all_childs = get(hObject,'Children');
		set(hObject,'pos',[pos_fig(1:3) new_height])
		for i = 1:length(all_childs)
			set(all_childs(i),'pos', get(all_childs(i),'pos') + [0 -(pos_fig(4) - new_height) 0 0])
		end
		set(handles.axes2,'Vis','off')
		handles.currAxes = zeros(1,3);
	end

	if (~handles.isRGB)
		handles.currAxes = 1;
		set(handles.figure1,'CurrentAxes', handles.axes1);
		if (isempty(hPolyg))
			localImhist(handles, img);
		else
			localROIImhist(handles, hPolyg);
		end
	else
		axes__ = [handles.axes1 handles.axes2 handles.axes3];
		for (k = 3:-1:1)
			set(handles.figure1,'CurrentAxes', axes__(k));
			handles.currAxes = k;
			if (isempty(hPolyg))
				localImhist(handles, img(:,:,k));
			else
				localROIImhist(handles, hPolyg);	% This fun will loop over RGB, so call it only once
				break
			end
		end
	end

	guidata(hObject, handles);
	set(hObject,'Vis','on');
	if (nargout),   varargout{1} = hObject;     end

% --------------------------------------------------------------------------
function [a, n, isScaled, top] = parse_inputs(varargin)

	a = varargin{1};        n = 256;

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
function img = localImresize(img,m)
% Specialized version of IMRESIZE to reduce the size of IMG by factor M
% Define old and new image sizes, and actual scaling
	so = [size(img,1) size(img,2)];		% old image size
	sn = max(floor(m*so(1:2)),1);		% new image size=(integer>0)
	sc = so./sn + eps;					% actual scaling between images, eps is needed
	img = img(floor(sc(1)/2:sc(1):end)+1,floor(sc(2)/2:sc(2):end)+1,:);

% --------------------------------------------------------------------------
function localROIImhist(handles, hPolyg)
% Call Mirone to do the ROI cropping and next call localImhist to do the histos
	img = mirone('ImageCrop_CB', guidata(handles.hMirFig), hPolyg);
	
	x = get(hPolyg, 'XData');		y = get(hPolyg, 'YData');
	if ( ~(numel(x) == 5 && (x(1) == x(2)) && (x(3) == x(4)) && (y(1) == y(4)) && (y(2) == y(3))) )
		handles.is_polyg = true;	
	end  

	if (ndims(img) == 3)
		axes__ = [handles.axes1 handles.axes2 handles.axes3];
		for (k = 3:-1:1)
			set(handles.figure1,'CurrentAxes', axes__(k));
			handles.currAxes = k;
			localImhist(handles, img(:,:,k));
		end
	else
		handles.currAxes = 1;
		set(handles.figure1,'CurrentAxes', handles.axes1);
		localImhist(handles, img);
	end

% --------------------------------------------------------------------------
function [y,x] = localImhist(handles, varargin)

	[a, n, isScaled, top] = parse_inputs(varargin{:});

	if islogical(a)
		y(2) = sum(a(:));
		y(1) = numel(a) - y(2);
		y = y';
	else
		y = imhistc(a, n, isScaled, top);       % Call MEX file to do work.
	end
	
	if (handles.is_polyg)		% In ROI cases we may have lots of whites (kida of NaNs)
		if (handles.validGrid),		y(1) = 0;
		else						y(end) = 0;
		end
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

	if (nargout == 0),    plot_result(x, y, handles);     end

% --------------------------------------------------------------------------
function plot_result(x, y, handles)

	hAxes = get(handles.figure1,'CurrentAxes');

	% Create a patch that encloses (x,y) limits
	x_max = max(x);     y_max = max(y);
	if (handles.isRGB)
		switch handles.currAxes
			case 1,		cor = [255 124 117] / 255;
			case 2,		cor = [120 255 114] / 255;
			case 3,		cor = [119 119 255] / 255;
		end
	else
		cor = [120 255 114] / 255;
	end
	patch([0 0 x_max x_max],[0 y_max y_max 0],cor);

 	set(hAxes,'XLim',[0 x_max],'YLim',[0 y_max]);
	localStem(hAxes,x,y)

	if (handles.isRGB && handles.currAxes < 3)  % We don't want XTicks on R & G axes of a RGB display
		set(hAxes, 'Xtick',[])
		label = get(hAxes,'YTickLabel');
		label(1) = ' ';     % We also don't want the '0' because is overlaps the next axes YTickLabel(end,:)
		set(hAxes,'YTickLabel',label)
	end

% --------------------------------------------------------------------------
function localStem(hAxes,x,y)
%   This is a hacked code from the stem function
%   Only strictly necessary code was kept

	% Set up data using fancing indexing
	[m,n] = size(x);
	xx = zeros(3*m,n);      yy = xx;
	xx(1:3:3*m,:) = x;      yy(2:3:3*m,:) = y;
	xx(2:3:3*m,:) = x;      yy(3:3:3*m,:) = NaN;
	xx(3:3:3*m,:) = NaN;

	line('XData',xx,'YData',yy,'color','k','Parent',hAxes);


% --- Creates and returns a handle to the GUI figure. 
function image_histo_LayoutFcn(h1)

set(h1, 'Color',get(0,'factoryUicontrolBackgroundColor'),...
'Position', [872 318 536 400], ...
'MenuBar','none',...
'Name','Image histogram',...
'NumberTitle','off',...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1');

axes('Parent',h1,...
'Units','pixels',...
'Position',[50 275 481 120],...
'Color',[0.9 1 1],...
'Tag','axes1');

axes('Parent',h1,...
'Units','pixels',...
'Position',[50 152 481 120],...
'Color',[0.9 1 1],...
'Tag','axes2');

axes('Parent',h1,...
'Units','pixels',...
'Position',[50 29 481 120],...
'Color',[0.9 1 1],...
'Tag','axes3');

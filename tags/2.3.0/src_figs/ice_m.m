function varargout = ice_m(varargin)
%ICE_EXPORT Interactive Color Editor.
%
%   OUT = ICE_EXPORT('Property Name', 'Property Value', ...) transforms an
%   image's color components based on interactively specified mapping
%   functions. Inputs are Property Name/Property Value pairs:  
%
%     Name            Value
%     ------------    -----------------------------------------------
%     'image'         An RGB or monochrome input image to be
%                     transformed by interactively specified mappings.
%     'space'         The color space of the components to be
%                     modified. Possible values are 'rgb', 'cmy',
%                     'hsi', 'hsv', 'ntsc' (or 'yiq'), 'ycrcb'. When
%                     omitted, the RGB color space is assumed.
%     'wait'          If 'on' (the default), OUT is the mapped input
%                     image and ICE_EXPORT returns to the calling function
%                     or workspace when closed. If 'off', OUT is the
%                     handle of the mapped input image and ICE_EXPORT
%                     returns immediately. 
%
%   EXAMPLES:
%     ice_m OR ice_m('wait', 'off')				% Demo user interface
%     ice_m('image', f)						% Map RGB or mono image
%     ice_m('image', f, 'space', 'hsv')		% Map HSV of RGB image
%     g = ice_m('image', f)					% Return mapped image
%     g = ice_m('image', f, 'wait', 'off');	% Return its handle
%
%   ICE_EXPORT displays one popup menu selectable mapping function at a time.
%   Each image component is mapped by a dedicated curve (e.g.,
%   R, G, or B) and then by an all-component curve (e.g., RGB). Each
%   curve's control points are depicted as circles that can be moved,
%   added, or deleted with a two- or three-button mouse:
%
%     Mouse Button    Editing Operation
%     ------------    -----------------------------------------------
%     Left            Move control point by pressing and dragging.
%     Middle          Add and position a control point by pressing
%                     and dragging. (Optionally Shift-Left)
%     Right           Delete a control point. (Optionally
%                     Control-Left) 
%
%   Checkboxes determine how mapping functions are computed, whether
%   the input image and reference pseudo- and full-color bars are
%   mapped, and the displayed reference curve information (e.g., PDF):
%
%     Checkbox        Function
%     ------------    -----------------------------------------------
%     Smooth          Checked for cubic spline (smooth curve)
%                     interpolation. If unchecked, piecewise linear.
%     Clamp Ends      Checked to force the starting and ending curve
%                     slopes in cubic spline interpolation to 0. No
%                     effect on piecewise linear.
%     Show PDF        Display probability density function(s) [i.e.,
%                     histogram(s)] of the image components affected
%                     by the mapping function.
%     Show CDF        Display cumulative distributions function(s)
%                     instead of PDFs.
%                     <Note: Show PDF/CDF are mutually exclusive.>
%     Map Image       If checked, image mapping is enabled; else
%                     not. 
%     Map Bars        If checked, pseudo- and full-color bar mapping
%                     is enabled; else display the unmapped bars (a
%                     gray wedge and hue wedge, respectively).
%
%   Mapping functions can be initialized via pushbuttons:
%
%     Button          Function
%     ------------    -----------------------------------------------
%     Reset           Init the currently displayed mapping function
%                     and uncheck all curve parameters.
%     Reset All       Initialize all mapping functions.

%   Copyright 2002-2004 R. C. Gonzalez, R. E. Woods, & S. L. Eddins
%   Digital Image Processing Using MATLAB, Prentice-Hall, 2004
%   $Revision: 1.5 $  $Date: 2003/11/21 14:16:10 $

%	WARNING. Some of the above may no longer be true since I strongly 
%	modified things to use cvlib_mex instead of IPT toolbox and did other
%	improvements, like: 	M-File changed by desGUIDE 
%
%	Joaquim Luis
 
	hObject = figure('Vis','off');
	ice_m_LayoutFcn(hObject);
	handles = guihandles(hObject);
 
	%  When ICE_EXPORT is opened, perform basic initialization (e.g., setup
	%  globals, ...) before it is made visible.

	% Set ICE_EXPORT globals to defaults.
	handles.updown = 'none';         % Mouse updown state
	handles.plotbox = [0 0 1 1];     % Plot area parameters in pixels
	handles.set1 = [0 0; 1 1];       % Curve 1 control points
	handles.set2 = [0 0; 1 1];       % Curve 2 control points
	handles.set3 = [0 0; 1 1];       % Curve 3 control points
	handles.set4 = [0 0; 1 1];       % Curve 4 control points
	handles.curve = 'set1';          % Structure name of selected curve
	handles.cindex = 1;              % Index of selected curve
	handles.node = 0;                % Index of selected control point
	handles.below = 1;               % Index of node below control point
	handles.above = 2;				% Index of node above control point
	handles.smooth = [0; 0; 0; 0];	% Curve smoothing states
	handles.slope = [0; 0; 0; 0];	% Curve end slope control states
	handles.cdf = [0; 0; 0; 0];		% Curve CDF states
	handles.pdf = [0; 0; 0; 0];		% Curve PDF states
	handles.output = [];			% Output image handle
	handles.df = [];				% Input PDFs and CDFs
	handles.colortype = 'rgb';		% Input image color space
	handles.input = [];				% Input image data
	handles.imagemap = 1;			% Image map enable
	handles.barmap = 1;				% Bar map enable
	handles.graybar = [];			% Pseudo (gray) bar image
	handles.colorbar = [];			% Color (hue) bar image
	handles.hMirFig = [];
	handles.first_bars = true;		% For creating colorbars
	handles.is2D = false;

	% Process Property Name/Property Value input argument pairs.
	wait = 'on';
	n_arg = numel(varargin);

	% If first arg is a Mirone handle, don't check later for the 'image' prop
	if (n_arg >= 1 && ishandle(varargin{1}) && strcmp(get(varargin{1}, 'type'), 'figure'))
		handles.hMirFig = varargin{1};
		hMirHand = guidata(handles.hMirFig);
		handles.input = get(hMirHand.hImg,'CData');
		handles.hMirImg = hMirHand.hImg;
		varargin(1) = [];
		n_arg = numel(varargin);
	end

	if (n_arg)
		for (i = 1:2:n_arg)
			if (n_arg == i),		break,		end
			switch lower(varargin{i})
				case 'image'
					if (isempty(handles.input)) 	% Otherwise ignore this because, ... Mirone feeded
						handles.input = varargin{i + 1};
					end

				case 'space'
					handles.colortype = lower(varargin{i + 1});
					list = cell(7,4);
					list(1,:) = {'CMY' 'Cyan' 'Magenta' 'Yellow'};
					list(2,:) = {'YCrCb' 'Luminance' 'Red Difference'  'Blue Difference'};
					list(3,:) = {'HSV' 'Hue' 'Saturation' 'Value'};
					list(4,:) = {'HSI' 'Hue' 'Saturation' 'Intensity'};
					list(5,:) = {'YIQ' 'Luminance' 'Hue' 'Saturation'};
					list(6,:) = {'La*b*' 'Luminance' 'a*' 'b*'};
					list(7,:) = {'RGB' 'Red' 'Green' 'Blue'};
					switch handles.colortype
						case 'cmy',			ind = 1;
						case 'ycrcb',		ind = 2;		handles.colortype = 'YCrCb';
						case 'hsv',			ind = 3;
						case 'hsi',			ind = 4;		handles.colortype = 'hls';
						case {'ntsc','yiq'},ind = 5;		handles.colortype = 'yiq';
						case 'lab',			ind = 6;		handles.colortype = 'Lab';
						otherwise,			ind = 7;		handles.colortype = 'rgb';
					end
					handles.space_list = list;
					set(handles.popup_component, 'String', list(ind,:));
					set(handles.popup_colorSpace,'String', list(:,1), 'Val', ind);

				case 'wait' 
					wait = lower(varargin{i + 1});
			end	
		end
	end

	if (ndims(handles.input) == 2)
		handles.input = cat(3, handles.input, handles.input, handles.input);
		handles.is2D = true;
	end

	if (isempty(handles.input))			% Just give the user some trash to play with.
		handles.input = uint8(rand(512,512,3) * 255);
		hImg = image(handles.input);
		handles.hMirImg = hImg;
		handles.hMirFig = get(get(hImg,'Parent'), 'Parent');
		hMirHand = guihandles(handles.hMirFig);
		hMirHand.origFig = handles.input;		% Copy of the dumb image
		guidata(handles.hMirFig, hMirHand)
	end

	% Create pseudo- and full-color mapping bars (grays and hues). Store
	% a color space converted 1x128x3 line of each bar for mapping.
	handles = set_colorSpace(handles);

	% Compute ICE's screen position and display image/graph.
	set(0, 'Units', 'pixels');      ssz = get(0, 'Screensize');
	uisz = get(handles.figure1, 'Position');
	fsz = get(handles.output, 'Position');
	bc = (fsz(4) - uisz(4)) / 3;
	if (bc > 0),		bc = bc + fsz(2);
	else				bc = fsz(2) + fsz(4) - uisz(4) - 10;  
	end
	lc = fsz(1) + (size(handles.input, 2) / 4) + (3 * fsz(3) / 4);
	lc = min(lc, ssz(3) - uisz(3)) + 5;
	set(handles.figure1, 'Position', [lc bc uisz(3:4)]);
	graph(handles);
	render(handles);

	% Add this figure handle to the carraças list
	plugedWin = getappdata(handles.hMirFig,'dependentFigs');
	plugedWin = [plugedWin hObject];
	setappdata(handles.hMirFig,'dependentFigs',plugedWin);

	% Update handles and make ICE wait before exit if required.
	guidata(hObject, handles);
	%if (strcmpi(wait, 'on')),	uiwait(handles.figure1);	end
	set(hObject,'Visible','on');
	if (nargout),   varargout{1} = hObject;     end

% ---------------------------------------------------------------------------------------------
function handles = set_colorSpace(handles)
% Do color space transforms, clamp to [0, 255], compute histograms
% and cumulative distribution functions, and create output figure.
% Create pseudo- and full-color mapping bars (grays and hues). Store
% a color space converted 1x128x3 line of each bar for mapping.
	xi = 0:1/127:1;		x = (0:1/6:1)';
	y = [1 1 0 0 0 1 1; 0 1 1 1 0 0 0; 0 0 0 1 1 1 0]';
 	gb = single(repmat(xi, [1 1 3]));		cb = reshape(single(interp1q(x, y, xi')), [1 128 3]);
	ctype = handles.colortype;
	if ~strcmp(ctype, 'rgb')
		if (strcmp(ctype ,'cmy'))
			gb = rgb2cmy(gb);		cb = rgb2cmy(cb);		handles.input = rgb2cmy(handles.input);
		elseif (strcmp(ctype ,'yiq'))
			gb = rgb2yiq(gb);		cb = rgb2yiq(cb);		handles.input = rgb2yiq(handles.input);
		else
			gb = cvlib_mex('color', gb, ['rgb2' ctype]);
			cb = cvlib_mex('color', cb, ['rgb2' ctype]);
			handles.input = cvlib_mex('color', handles.input, ['rgb2' ctype]);
		end
		
		cb = max(0, cb);		gb = max(0, gb);
		if (size(cb,3) == 1)
			gb = cat(3, gb, gb, gb);		cb = cat(3, cb, cb, cb);
		end
	end
	if (ctype(1) == 'h')
		cb(:,:,1) = single(double(cb(:,:,1)) / 360);
	end
	handles.graybar = round(255 * double(gb));     handles.colorbar = round(255 * double(cb));

	for (i = 1:3)
		color = handles.input(:, :, i);
		df = double(histc(color(:), 0:255))';
		%df = hist(double(color(:)), 0:255);
		handles.df = [handles.df; df / max(df(:))];
		df = df / sum(df(:));   df = cumsum(df);
		handles.df = [handles.df; df];
	end
	if (~isempty(handles.hMirFig)),		handles.output = handles.hMirFig;
	else								handles.output = figure;
	end
	%guidata(handles.figure1, handles)

%-------------------------------------------------------------------%
function ice_m_WindowButtonDownFcn(hObject, eventdata)
%  Start mapping function control point editing. Do move, add, or
%  delete for left, middle, and right button mouse clicks ('normal',
%  'extend', and 'alt' cases) over plot area.

handles = guidata(hObject);
handles.plotbox = get(handles.curve_axes, 'Position');
[inplot, x, y] = cursor(hObject, handles);
if inplot
	nodes = handles.(handles.curve);
	i = find(x >= nodes(:, 1));      below = max(i);
	above = min(below + 1, size(nodes, 1));
	if (x - nodes(below, 1)) > (nodes(above, 1) - x)    
		node = above;
	else  
		node = below;    
	end
	deletednode = 0;
        
	switch get(hObject, 'SelectionType')
		case 'normal'
			if (node == above),				above = min(above + 1, size(nodes, 1));
			elseif (node == below)			below = max(below - 1, 1);      
			end
			if (node == size(nodes, 1)),	below = above;
			elseif (node == 1)				above = below;      
			end
			if (x > nodes(above, 1)),		x = nodes(above, 1);
			elseif (x < nodes(below, 1))	x = nodes(below, 1);    
			end
			handles.node = node;     handles.updown = 'down';
			handles.below = below;   handles.above = above;
			nodes(node, :) = [x y];
		case 'extend'
			if isempty(find(nodes(:, 1) == x))
				nodes = [nodes(1:below, :); [x y]; nodes(above:end, :)];
				handles.node = above;   handles.updown = 'down';
				handles.below = below;  handles.above = above + 1;
			end
		case 'alt'
			if (node ~= 1) && (node ~= size(nodes, 1))
				nodes(node, :) = [];   deletednode = 1;    
			end
			handles.node = 0;
			set(handles.input_text, 'String', '');
			set(handles.output_text, 'String', '');
	end
    
   handles = setfield(handles, handles.curve, nodes);
   guidata(hObject, handles);
   graph(handles);
   if (deletednode),	render(handles);	end
end

%-------------------------------------------------------------------%
function ice_m_WindowButtonMotionFcn(hObject, eventdata)
%  Do nothing unless a mouse 'down' event has occurred. If it has,
%  modify control point and make new mapping function. 
	handles = guidata(hObject);
	if (~strcmpi(handles.updown, 'down')),		return,		end
	[inplot, x, y] = cursor(hObject, handles);
	if inplot
		nodes = handles.(handles.curve);
		nudge = handles.smooth(handles.cindex) / 256;
		if (handles.node ~= 1) && (handles.node ~= size(nodes, 1))
			if (x >= nodes(handles.above, 1)),			x = nodes(handles.above, 1) - nudge;
			elseif (x <= nodes(handles.below, 1))		x = nodes(handles.below, 1) + nudge;    
			end
		else
			if (x > nodes(handles.above, 1)),			x = nodes(handles.above, 1);
			elseif (x < nodes(handles.below, 1))		x = nodes(handles.below, 1);    
			end
		end
		nodes(handles.node, :) = [x y];
		handles = setfield(handles, handles.curve, nodes);
		guidata(hObject, handles);
		graph(handles);
	end

%-------------------------------------------------------------------%
function ice_m_WindowButtonUpFcn(hObject, eventdata)
%  Terminate ongoing control point move or add operation. Clear
%  coordinate text below plot and update display.
	handles = guidata(hObject);
	update = strcmpi(handles.updown, 'down');
	handles.updown = 'up';      handles.node = 0;
	guidata(hObject, handles);
	if update
		set(handles.input_text, 'String', '');
		set(handles.output_text, 'String', '');
		render(handles);
	end

%-------------------------------------------------------------------%
function popup_component_CB(hObject, handles, noRedraw)
%  Accept color component selection, update component specific
%  parameters on GUI, and draw the selected mapping function.

	if (nargin == 2),	noRedraw = false;	end
	c = get(hObject, 'Value');
	handles.cindex = c;
	handles.curve = strcat('set', num2str(c));
	guidata(hObject, handles);
	set(handles.check_smooth, 'Value', handles.smooth(c));
	set(handles.check_slope, 'Value', handles.slope(c));
	set(handles.check_pdf, 'Value', handles.pdf(c));
	set(handles.check_cdf, 'Value', handles.cdf(c));
	if (~noRedraw),		graph(handles),		end

%-------------------------------------------------------------------%
function popup_colorSpace_CB(hObject, handles)
% Reinitiate with a new color space

	set(handles.popup_component,'Val',1)
	popup_component_CB(handles.popup_component, handles, true)		% true e Do  not call graph
	handles = guidata(handles.figure1);
	push_resetall_CB(handles.push_resetall, handles, true)			% true e noRedraw
	handles = guidata(handles.figure1);
	%push_reset_CB(handles.push_reset, [], handles, true)
	val = get(hObject,'Value');
	set(handles.popup_component, 'String', handles.space_list(val,:))
	handles.colortype = lower(handles.space_list{val,1});
	if (numel(handles.colortype) == 5),		handles.colortype = 'YCrCb';	end
	if (strcmp(handles.colortype, 'hsi')),	handles.colortype = 'hls';		end
	hMirHand = guidata(handles.hMirImg);
	handles.input = hMirHand.origFig;		% Start with the original data
	if (handles.is2D)						% Orig file is 2D. Need the expand
		handles.input = cat(3, handles.input, handles.input, handles.input);
	end
	handles = set_colorSpace(handles);
	graph(handles);
	render(handles);
	guidata(handles.figure1, handles)

%-------------------------------------------------------------------%
function check_smooth_CB(hObject, handles)
%  Accept smoothing parameter for currently selected color
%  component and redraw mapping function.
	if get(hObject, 'Value')
		handles.smooth(handles.cindex) = 1;
		nodes = handles.(handles.curve);
		nodes = spreadout(nodes);
		handles = setfield(handles, handles.curve, nodes);
	else  
		handles.smooth(handles.cindex) = 0;    
	end
	guidata(hObject, handles);
	set(handles.figure1, 'Pointer', 'watch');
	graph(handles);      render(handles);
	set(handles.figure1, 'Pointer', 'arrow');

%-------------------------------------------------------------------%
function push_reset_CB(hObject, handles)
% Init all display parameters for currently selected color component and redraw it.
	handles = setfield(handles, handles.curve, [0 0; 1 1]);
	c = handles.cindex;
	handles.smooth(c) = 0;		handles.slope(c) = 0;
	handles.pdf(c) = 0;			handles.cdf(c) = 0;
	set([handles.check_smooth handles.check_slope handles.check_pdf handles.check_cdf], 'Val', 0)
	guidata(hObject, handles);
	set(handles.figure1, 'Pointer', 'watch');
	graph(handles);		render(handles);
	set(handles.figure1, 'Pointer', 'arrow');

%-------------------------------------------------------------------%
function push_resetall_CB(hObject, handles, noRedraw)
% Init display parameters for color components and redraw display.
% Third arg, 'noRedraw', is used when one changes the color space and
% the redrawing is requested somewhere else.

	if (nargin == 2),	noRedraw = false;	end
	for c = 1:4
		handles.smooth(c) = 0;		handles.slope(c) = 0;
		handles.pdf(c) = 0;			handles.cdf(c) = 0;
		handles = setfield(handles, ['set' num2str(c)], [0 0; 1 1]);
	end
	set([handles.check_smooth handles.check_slope handles.check_pdf handles.check_cdf], 'Val', 0)
	guidata(hObject, handles);
	graph(handles);
	if (~noRedraw)
		set(handles.figure1, 'Pointer', 'watch');
		render(handles);
		set(handles.figure1, 'Pointer', 'arrow');
	end

%-------------------------------------------------------------------%
function check_slope_CB(hObject, handles)
%  Accept slope clamp for currently selected color component and
%  draw function if smoothing is on.

	if get(hObject, 'Value'),	handles.slope(handles.cindex) = 1;
	else						handles.slope(handles.cindex) = 0;    
	end
	guidata(hObject, handles);
	if handles.smooth(handles.cindex)
		set(handles.figure1, 'Pointer', 'watch');
		graph(handles);			render(handles);
		set(handles.figure1, 'Pointer', 'arrow');   
	end

%-------------------------------------------------------------------%
function check_pdf_CB(hObject, handles)
%  Accept PDF (probability density function or histogram) display
%  parameter for currently selected color component and redraw
%  mapping function if smoothing is on. If set, clear CDF display.
	if get(hObject, 'Value')
		handles.pdf(handles.cindex) = 1;
		set(handles.check_cdf, 'Value', 0);
		handles.cdf(handles.cindex) = 0;
	else  
		handles.pdf(handles.cindex) = 0;    
	end
	guidata(handles.figure1, handles);
	graph(handles);

%-------------------------------------------------------------------%
function check_cdf_CB(hObject, handles)
%  Accept CDF (cumulative distribution function) display parameter
%  for selected color component and redraw mapping function if
%  smoothing is on. If set, clear CDF display.
	if get(hObject, 'Value')
		handles.cdf(handles.cindex) = 1;
		set(handles.check_pdf, 'Value', 0);
		handles.pdf(handles.cindex) = 0;
	else   
		handles.cdf(handles.cindex) = 0;    
	end
	guidata(handles.figure1, handles);      graph(handles);

%-------------------------------------------------------------------%
function check_mapbar_CB(hObject, handles)
%  Accept changes to bar map enable state and redraw bars.
	handles.barmap = get(hObject, 'Value');
	guidata(handles.figure1, handles);      render(handles);

%-------------------------------------------------------------------%
function check_mapimage_CB(hObject, handles)
%  Accept changes to the image map state and redraw image.
	handles.imagemap = get(hObject, 'Value');
	guidata(handles.figure1, handles);      render(handles);

%-------------------------------------------------------------------%
function graph(handles)
%  Interpolate and plot mapping functions and optional reference PDF(s) or CDF(s).

nodes = handles.(handles.curve);
c = handles.cindex;     dfx = 0:1/255:1;
colors = ['k' 'r' 'g' 'b'];
    
% For piecewise linear interpolation, plot a map, map + PDF/CDF, or map + 3 PDFs/CDFs.
if ~handles.smooth(handles.cindex)
	if (~handles.pdf(c) && ~handles.cdf(c)) || (size(handles.df, 2) == 0)
		plot(nodes(:, 1), nodes(:, 2), 'b-', nodes(:, 1),  nodes(:, 2), 'ko', 'Parent', handles.curve_axes);
	elseif c > 1
		i = 2 * c - 2 - handles.pdf(c);
		plot(dfx, handles.df(i, :), [colors(c) '-'], nodes(:, 1), nodes(:, 2), 'k-', nodes(:, 1), nodes(:, 2), 'ko', 'Parent', handles.curve_axes);
	elseif c == 1
		i = handles.cdf(c);
		plot(dfx, handles.df(i + 1, :), 'r-', dfx, handles.df(i + 3, :), 'g-', dfx, handles.df(i + 5, :), 'b-', ...
			nodes(:, 1), nodes(:, 2), 'k-', nodes(:, 1), nodes(:, 2), 'ko', 'Parent', handles.curve_axes);
	end

% Do the same for smooth (cubic spline) interpolations.
else
	x = 0:0.01:1;
	if ~handles.slope(handles.cindex)
		y = spline(nodes(:, 1), nodes(:, 2), x);
	else    
		y = spline(nodes(:, 1), [0; nodes(:, 2); 0], x);    
	end
	y(y > 1) = 1;
	y(y < 0) = 0;
	
	if (~handles.pdf(c) && ~handles.cdf(c)) || (size(handles.df, 2) == 0)
		plot(nodes(:, 1), nodes(:, 2), 'ko',  x, y, 'b-', 'Parent', handles.curve_axes);
	elseif c > 1
		i = 2 * c - 2 - handles.pdf(c);
		plot(dfx, handles.df(i, :), [colors(c) '-'], nodes(:, 1), nodes(:, 2), 'ko', x, y, 'k-', 'Parent', handles.curve_axes);
	elseif c == 1
		i = handles.cdf(c);
		plot(dfx, handles.df(i + 1, :), 'r-', dfx, handles.df(i + 3, :), 'g-', dfx, handles.df(i + 5, :), 'b-', ...
			nodes(:, 1), nodes(:, 2), 'ko', x, y, 'k-', 'Parent', handles.curve_axes);
	end
end

% Put legend if more than two curves are shown.
s = handles.colortype;
if (strcmp(s, 'ntsc')),		s = 'yiq';		end
if (c == 1) && (handles.pdf(c) || handles.cdf(c))
	s1 = ['-- ' upper(s(1))];
	if length(s) == 3
       s2 = ['-- ' upper(s(2))];		s3 = ['-- ' upper(s(3))];
	else
       s2 = ['-- ' upper(s(2)) s(3)];	s3 = ['-- ' upper(s(4)) s(5)];
	end
else    
	s1 = '';    s2 = '';    s3 = '';    
end
set(handles.red_text, 'String', s1);
set(handles.green_text, 'String', s2);
set(handles.blue_text, 'String', s3);
    
%-------------------------------------------------------------------%
function [inplot, x, y] = cursor(h, handles)
%  Translate the mouse position to a coordinate with respect to
%  the current plot area, check for the mouse in the area and if so
%  save the location and write the coordinates below the plot.

	p = get(h, 'CurrentPoint');
	x = (p(1, 1) - handles.plotbox(1)) / handles.plotbox(3);
	y = (p(1, 2) - handles.plotbox(2)) / handles.plotbox(4);
	if (x > 1.05 || x < -0.05 || y > 1.05 || y < -0.05)
		inplot = 0;
	else
		x = min(x, 1);      x = max(x, 0);
		y = min(y, 1);      y = max(y, 0);
		x = round(256 * x) / 256;
		inplot = 1;
		set(handles.input_text, 'String', num2str(x, 3));
		set(handles.output_text, 'String', num2str(y, 3));
	end

%-----------------------------------------------------------------------------
function render(handles)
%  Map the input image and bar components and convert them to RGB (if needed) and display.

	set(handles.figure1, 'Interruptible', 'off', 'Pointer', 'watch');
	ygb = handles.graybar;		ycb = handles.colorbar;
	yi = handles.input;			mapon = handles.barmap;
	imageon = handles.imagemap;

	for i = 2:4
		nodes = getfield(handles, ['set' num2str(i)]);
		t = lut(nodes, handles.smooth(i), handles.slope(i));
		if imageon  
			yi(:, :, i - 1) = t(cvlib_mex('addS', yi(:, :, i - 1), 1));   
		end
		if mapon
			ygb(:, :, i - 1) = t(ygb(:, :, i - 1) + 1); 
			ycb(:, :, i - 1) = t(ycb(:, :, i - 1) + 1);     
		end
	end
	t = uint8( lut(handles.set1, handles.smooth(1), handles.slope(1)) );
	if imageon
		for (i = 1:size(yi,3))		% Use a loop due to a BUG in OpenCV
			yi(:,:,i) = cvlib_mex('addS', yi(:,:,i), 1);
		end
		yi = t(yi);
	end
	if mapon    
		ygb = t(ygb + 1);     ycb = t(ycb + 1);
		ygb = uint8(ygb);     ycb = uint8(ycb);
	end

	ctype = handles.colortype;
	if ~strcmp(ctype, 'rgb')
		if (strcmp(ctype ,'cmy')),		yi = cmy2rgb(yi);
		elseif (strcmp(ctype ,'yiq')),	yi = yiq2rgb(yi);
		else							yi = cvlib_mex('color', yi, [ctype '2rgb']);
		end

		if mapon
			if (strcmp(ctype ,'cmy')),		ygb = cmy2rgb(ygb);		ycb = cmy2rgb(ycb);
			elseif (strcmp(ctype ,'yiq')),	ygb = yiq2rgb(ygb);		ycb = yiq2rgb(ycb);
			else							ygb = cvlib_mex('color', ygb, [ctype '2rgb']);
											ycb = cvlib_mex('color', ycb, [ctype '2rgb']);
			end
		end
	end

	set(handles.hMirImg,'CData',yi)
	if mapon
		ygb = repmat(ygb, [32 1 1]);    ycb = repmat(ycb, [32 1 1]);
		if (handles.first_bars)
			handles.hImg_gray_axes  = image(ygb, 'Parent',handles.gray_axes);
			handles.hImg_color_axes = image(ycb, 'Parent',handles.color_axes);
			handles.first_bars = false;
			guidata(handles.figure1, handles)
			set(handles.gray_axes,  'XTick',[], 'YTick',[]);
			set(handles.color_axes, 'XTick',[], 'YTick',[]);
		else
			set(handles.hImg_gray_axes,  'CData',ygb);
			set(handles.hImg_color_axes, 'CData',ycb);
		end
	end
	set(handles.figure1, 'Interruptible', 'on', 'Pointer', 'arrow');

%-------------------------------------------------------------------%
function t = lut(nodes, smooth, slope)
%  Create a 256 element mapping function from a set of control points.
%  The output values are integers in the interval [0, 255]. Use piecewise
%  linear or cubic spline with or without zero end slope interpolation.

	t = 255 * nodes;    i = 0:255;
	if ~smooth    
		t = [t; 256 256];   t = interp1q(t(:, 1), t(:, 2), i');
	else
		if ~slope   
			t = spline(t(:, 1), t(:, 2), i);
		else        
			t = spline(t(:, 1), [0; t(:, 2); 0], i);    
		end
	end
	t = round(t);		t = max(0, t);			t = min(255, t);

%-------------------------------------------------------------------%
function out = spreadout(in)
%  Make all x values unique.

% Scan forward for non-unique x's and bump the higher indexed x--
% but don't exceed 1. Scan the entire range.
nudge = 1 / 256;
for i = 2:size(in, 1) - 1
	if in(i, 1) <= in(i - 1, 1)
       in(i, 1) = min(in(i - 1, 1) + nudge, 1);  
	end
end

% Scan in reverse for non-unique x's and decrease the lower indexed
% x -- but don't go below 0. Stop on the first non-unique pair.
if in(end, 1) == in(end - 1, 1)
	for i = size(in, 1):-1:2
		if in(i, 1) <= in(i - 1, 1)
			in(i - 1, 1) = max(in(i, 1) - nudge, 0);
		else 
			break
		end
	end
end

% If the first two x's are now the same, init the curve.
if in(1, 1) == in(2, 1)     
	in = [0 0; 1 1];    
end
out = in;

%-------------------------------------------------------------------%
function g = rgb2cmy(f)
% Convert RGB to CMY, using a kind of imcomplement
	switch class(f)
		case 'logical'
			g = ~f;
		case 'double'
			g = 1 - f;
		case {'int8' 'int16' 'single' 'int32'}
			if (ndims(f) > 2)
				for (k = 1:size(f,3))
					cvlib_mex('subRS', f(:,:,k), 1);		% Due to a BUG in OpenCV
				end
				g = f;
			else
				g = cvlib_mex('subRS', f, 1);
			end
		case 'uint8'
			lut = uint8(255:-1:0);
			g = intlutc(f,lut);
		case 'uint16'
			lut = uint16(65535:-1:0);
			g = intlutc(f,lut);
		case 'uint32'
			error('UInt32 is currently not supported image class.')
		otherwise
			error('Invalid input image class.')
	end

%-------------------------------------------------------------------%
function g = cmy2rgb(f)
%   Convert CMY to RGB
	g = rgb2cmy(f);

% --------------------------------------------------------------------
function g = rgb2yiq(f)
%RGB2YIQ Convert RGB values to NTSC (YIQ) colorspace.
%   YIQ = RGB2YIQ(RGB) converts the truecolor image RGB to the
%   equivalent NTSC image YIQ that contains the NTSC luminance (Y)
%	and chrominance (I and Q) color components as columns that are
%	equivalent to the colors in the RGB colormap.
%
%   Class Support
%   -------------
%   RGB can be uint8,uint16 or double; The output is double.

	class_orig = class(f);
	f = local_im2double(f);
	T = [1.0 0.956 0.621; 1.0 -0.272 -0.647; 1.0 -1.106 1.703].';
	[m n thirdD] = size(f);
	g = reshape(reshape(f,m*n,thirdD)/T,m,n,thirdD);
	if (strcmp(class_orig, 'uint8'))
		g = uint8(g * 255);
	end

% --------------------------------------------------------------------
function [r,g,b] = yiq2rgb(y,i,q)
%   RGB = YIQ2RGB(YIQ) converts the NTSC image YIQ to the equivalent truecolor image RGB.
%
%   Class Support
%   -------------
%   The input image must be of class double. The output is of class double.

	T = [1.0 0.956 0.621; 1.0 -0.272 -0.647; 1.0 -1.106 1.703];

	threeD = (ndims(y)==3); % Determine if input includes a 3-D array.

	if threeD
		y = local_im2double(y);
		rgb = reshape(y(:),size(y,1)*size(y,2),3)*T';
	else
		y = local_im2double(y);		i = local_im2double(i);		q = local_im2double(q);
		rgb = [y(:) i(:) q(:)]*T';
	end

	% Make sure the rgb values are between 0.0 and 1.0
	rgb = max(0,rgb);
	d = find(any(rgb' > 1));
	rgb(d,:) = rgb(d,:)./(max(rgb(d,:)')'*ones(1,3));

	if (nargout==0 || nargout==1)
		r = reshape(rgb,size(y,1),size(y,2),3);
	else
		[m,n] = size(y);
		r = reshape(rgb(:,1),m,n);
		g = reshape(rgb(:,2),m,n);
		b = reshape(rgb(:,3),m,n);
	end

% -------------------------------------------------------------------
function d = local_im2double(img)
% Convert image to double precision. Local version without any error check
	if isa(img, 'double')
		d = img;
	elseif (isa(img, 'logical') || isa(img, 'single'))
		d = double(img);
	elseif isa(img, 'uint8') || isa(img, 'uint16')
		if isa(img, 'uint8'),	d = double(img)/255;	
		else					d = double(img)/65535;
		end
	end

% ----------------------------------------------------------------
function ice_m_LayoutFcn(h1)

set(h1,...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'DoubleBuffer','on',...
'MenuBar','none',...
'Name','ICE - Interactive Color Editor',...
'NumberTitle','off',...
'Position',[5 868 439 365],...
'RendererMode','manual',...
'WindowButtonDownFcn',@ice_m_WindowButtonDownFcn,...
'WindowButtonMotionFcn',@ice_m_WindowButtonMotionFcn,...
'WindowButtonUpFcn',@ice_m_WindowButtonUpFcn,...
'HandleVisibility','callback',...
'Tag','figure1');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'FontSize',12,...
'FontWeight','bold',...
'HorizontalAlignment','left',...
'Position',[10 342 100 20],...
'String','Color Space',...
'Style','text');

uicontrol('Parent',h1, 'Position',[309 179 120 152], 'Style','frame', 'Tag','frame1');

axes('Parent',h1,...
'Units','pixels',...
'Box','on',...
'CameraPosition',[0.5 0.5 9.16025403784439],...
'FontSize',8,...
'NextPlot','replacechildren',...
'Position',[27 74 274 261],...
'XColor',get(0,'defaultaxesXColor'),...
'XGrid','on',...
'XLim',get(0,'defaultaxesXLim'),...
'XLimMode','manual',...
'XTick',[0 0.25 0.5 0.75 1],...
'XTickMode','manual',...
'YColor',get(0,'defaultaxesYColor'),...
'YGrid','on',...
'YLim',get(0,'defaultaxesYLim'),...
'YLimMode','manual',...
'YTick',[0 0.25 0.5 0.75 1],...
'YTickMode','manual',...
'Tag','curve_axes');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'FontSize',12,...
'FontWeight','bold',...
'HorizontalAlignment','left',...
'Position',[207 342 100 20],...
'String','Component:',...
'Style','text',...
'Tag','text1');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'FontSize',10,...
'HorizontalAlignment','left',...
'Position',[42 24 50 25],...
'String','Input:',...
'Style','text');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'FontSize',10,...
'HorizontalAlignment','left',...
'Position',[33 6 50 20],...
'String','Output:',...
'Style','text');

uicontrol('Parent',h1,...
'Call',@ice_m_uiCB,...
'FontName','Helvetica',...
'FontSize',10,...
'HorizontalAlignment','left',...
'Position',[323 296 87 24],...
'String','Smooth',...
'Style','checkbox',...
'Tag','check_smooth');

uicontrol('Parent',h1,...
'Call',@ice_m_uiCB,...
'FontName','Helvetica',...
'FontSize',10,...
'Position',[320 191 98 30],...
'String','Reset',...
'Tag','push_reset');

uicontrol('Parent',h1,...
'FontAngle','italic',...
'FontName','Helvetica',...
'FontSize',12,...
'Position',[333 320 67 22],...
'String','Curve',...
'Style','text');

uicontrol('Parent',h1,...
'BackgroundColor',[0.831372559070587 0.815686285495758 0.7843137383461],...
'FontSize',12,...
'ForegroundColor',[1 0 0],...
'Position',[83 32 50 16],...
'Style','text',...
'Tag','input_text');

uicontrol('Parent',h1,...
'BackgroundColor',[0.831372559070587 0.815686285495758 0.7843137383461],...
'FontSize',12,...
'ForegroundColor',[1 0 0],...
'Position',[83 9 50 16],...
'Style','text',...
'Tag','output_text');

uicontrol('Parent',h1,...
'Call',@ice_m_uiCB,...
'FontName','Helvetica',...
'FontSize',10,...
'HorizontalAlignment','left',...
'Position',[323 273 99 24],...
'String','Clamp Ends',...
'Style','checkbox',...
'Tag','check_slope');

uicontrol('Parent',h1,...
'Call',@ice_m_uiCB,...
'FontName','Helvetica',...
'FontSize',10,...
'Position',[311 14 120 30],...
'String','Reset All',...
'Tag','push_resetall');

uicontrol('Parent',h1,...
'Call',@ice_m_uiCB,...
'FontName','Helvetica',...
'FontSize',10,...
'HorizontalAlignment','left',...
'Position',[324 250 87 24],...
'String','Show PDF',...
'Style','checkbox',...
'Tag','check_pdf');

uicontrol('Parent',h1,...
'Call',@ice_m_uiCB,...
'FontName','Helvetica',...
'FontSize',10,...
'HorizontalAlignment','left',...
'Position',[324 227 87 24],...
'String','Show CDF',...
'Style','checkbox',...
'Tag','check_cdf');

axes('Parent',h1,...
'Units','pixels',...
'CameraPosition',[0.5 0.5 9.16025403784439],...
'Color',get(0,'defaultaxesColor'),...
'ColorOrder',get(0,'defaultaxesColorOrder'),...
'Position',[305 136 128 32],...
'XTick',[],...
'XTickLabel','',...
'XTickLabelMode','manual',...
'XTickMode','manual',...
'YTick',[],...
'YTickMode','manual',...
'Tag','gray_axes');

axes('Parent',h1,...
'Units','pixels',...
'CameraPosition',[0.5 0.5 9.1602540378444],...
'Color',get(0,'defaultaxesColor'),...
'ColorOrder',get(0,'defaultaxesColorOrder'),...
'Position',[305 78 128 32],...
'XTick',[],...
'XTickMode','manual',...
'YTick',[],...
'YTickMode','manual',...
'Tag','color_axes');

uicontrol('Parent',h1,...
'FontAngle','italic',...
'FontName','Helvetica',...
'FontSize',10,...
'Position',[309 117 120 16],...
'String','Pseudo-color Bar',...
'Style','text');

uicontrol('Parent',h1,...
'FontAngle','italic',...
'FontName','Helvetica',...
'FontSize',10,...
'Position',[324 59 89 16],...
'String','Full-color Bar',...
'Style','text');

uicontrol('Parent',h1,...
'Call',@ice_m_uiCB,...
'FontName','Helvetica',...
'FontSize',10,...
'Position',[193 28 80 24],...
'String','Map Bars',...
'Style','checkbox',...
'Value',1,...
'Tag','check_mapbar');

uicontrol('Parent',h1,...
'Call',@ice_m_uiCB,...
'FontName','Helvetica',...
'FontSize',10,...
'Position',[193 5 90 24],...
'String','Map Image',...
'Style','checkbox',...
'Value',1,...
'Tag','check_mapimage');

uicontrol('Parent',h1,...
'FontSize',12,...
'ForegroundColor',[0 0 1],...
'HorizontalAlignment','left',...
'Position',[152 0 40 20],...
'Style','text',...
'Tag','blue_text');

uicontrol('Parent',h1,...
'FontSize',12,...
'ForegroundColor',[0 1 0],...
'HorizontalAlignment','left',...
'Position',[152 17 40 20],...
'Style','text',...
'Tag','green_text');

uicontrol('Parent',h1,...
'FontSize',12,...
'ForegroundColor',[1 0 0],...
'HorizontalAlignment','left',...
'Position',[152 34 40 20],...
'Style','text',...
'Tag','red_text');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@ice_m_uiCB,...
'FontSize',10,...
'Position',[107 341 80 23],...
'String','RGB',...
'Style','popupmenu',...
'Value',1,...
'Tag','popup_colorSpace');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@ice_m_uiCB,...
'FontSize',10,...
'Position',[305 340 125 23],...
'String',{'RGB'; 'Red'; 'Green'; 'Blue'},...
'Style','popupmenu',...
'Value',1,...
'Tag','popup_component');

function ice_m_uiCB(hObject, eventdata)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));
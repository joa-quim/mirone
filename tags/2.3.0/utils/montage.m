function montage(I, varargin)
%MONTAGE  Display multiple images as a montage of subplots
%
% Examples:
%   montage
%   montage(I)
%   montage(I, map)
%   montage(..., param1, value1, param2, value2, ...)
%
% This function displays multiple images, in a stack or cell array or
% defined by filenames in a cell array or simply in the current directory,
% in a grid of subplots in the current figure.
%
% The size of grid is calculated or user defined. Images not fitting in the
% grid can be scrolled through using the scroll keys. This allows fast
% scrolling through a movie or image stack, e.g.
%    montage(imstack, 'Size', 1)
% creates a single frame montage with images scrolled through using arrow keys
%
% This function is designed to replace the MONTAGE function provided with
% the Image Processing Toolbox (IPT). It is syntax compatible, but operates
% differently, and while results are generally comparable given identical
% syntaxes, this is not guaranteed.
%
% Differences from the IPT version are:
%    - The IPT is not required!
%    - Images are placed in subplots, so can be zoomed separately.
%    - Small images are properly enlarged on screen.
%    - Gaps can be placed between images.
%    - Images can be viewed on a grid smaller than the number of images.
%    - Surplus images can be viewed by scrolling through pages.
%    - A directory of images can be viewed easily.
%    - It cannot return an image handle (as there are multiple images)
%
% Keys:
%    Up - Back a row.
%    Down - Forward a row.
%    Left - Back a page (or column if there is only one row).
%    Right - Forward a page (or column if there is only one row).
%    Shift - 2 x speed.
%    Ctrl - 4 x speed.
%    Shift + Ctrl - 8 x speed.
%
% IN:
%   I - MxNxCxP array of images, or 1xP cell array. C is 1 for indexed
%       images or 3 for RGB images. P is the number of images. If I is a
%       cell array then each cell must contain an image or image filename.
%       If I is empty then all the images in the current directory are
%       used. Default: [].
%   map - Kx3 colormap to be used with indexed images. Default: gray(256).
%   Optional parameters - name, value parameter pairs for the following:
%      'Size' - [H W] size of grid to display image on. If only H is given
%               then W = H. If either H or W is NaN then the number of rows
%               or columns is chosen such that all images fit. If both H
%               and W are NaN or the array is empty then the size of grid
%               is chosen to fit all images in as large as possible.
%               Default: [].
%      'Indices' - 1xL list of indices of images to display. Default: 1:P.
%      'Border' - [B R] borders to give each image top and bottom (B) and
%                 left and right (R), to space out images. Borders are
%                 normalized to the subplot size, i.e. B = 0.01 gives a border
%                 1% of the height of each subplot. If only B is given, R =
%                 B. Default: 0.01.
%      'DisplayRange' - [LOW HIGH] display range for indexed images.
%                       Default: [min(I(:)) max(I(:))].
%      'Map' - Kx3 colormap or (additionally from above) name of MATLAB
%              colormap, for use with indexed images. Default: gray(256).

% $Id: montage.m,v 1.7 2009/02/25 16:39:01 ojw Exp $

[map layout gap indices lims] = parse_inputs(varargin);

if (nargin == 0 || isempty(I))		% Read in all the images in the directory
    I = get_im_names;
    if isempty(I),	return,		end	% No images found
end

if (isnumeric(I))
    [y x c n] = size(I);
    if isempty(lims)
        lims = [min(reshape(I, numel(I), 1)) max(reshape(I, numel(I), 1))];
    elseif isequal(0, lims)
        lims = default_limits(I);
    end
    if ( c == 3 && (isa(I,'single') || isa(I,'double')) )
        I = uint8(I * 256 - 0.5);
        lims = round(lims * 256 - 0.5);
    end
    I = squeeze(num2cell(I, [1 2 3]));

elseif iscell(I)
    A = I{1};
    if ischar(A)
        A = imread_rgb(A);
        I{1} = imread_rgb(A);
    end
    n = numel(I);
    % Assume all images are the same size and type as the first
    [y x c] = size(A);
    if (isempty(lims) || isequal(0, lims)),		lims = default_limits(A);	end

elseif (isa(I,'struct'))
	if (~isfield(I,'cdata') || ~isfield(I,'colormap'))
		error('Input must be a MATLAB movie.');
	end
    if (isempty(lims) || isequal(0, lims)),		lims = default_limits(I(1).cdata);	end
	map = I(1).colormap;
    [y x c] = size(I(1).cdata);
	if (c == 1),	I = squeeze(num2cell(cat(3,I.cdata), [1 2]));		% Indexed
	else			I = squeeze(num2cell(cat(4,I.cdata), [1 2 3]));		% RGB
	end
    n = numel(I);
else
	error('I not of recognized type.');
end

% Select indexed images
if ~isequal(indices, -1)
    I = I(indices);
    n = numel(I);
end

layout = choose_layout(n, y, x, layout);		% Compute a good layout

% Create a data structure to store the data in
num = prod(layout);
state.n = num * ceil(n / num);
state.h = zeros(layout);
I = [I(:); cell(state.n-n, 1)];

% Get and clear the figure
%fig = gcf;		clf(fig);
fig = figure('NumberTitle','off', 'Name', 'Montage');

% Set the figure size well
MonSz = get(0, 'ScreenSize');
MaxSz = MonSz(3:4) - [20 120];
ImSz = layout([2 1]) .* [x y] ./ (1 - 2 * gap([end 1]));
RescaleFactor = min(MaxSz ./ ImSz);
if (RescaleFactor > 1)				% Integer scale for enlarging, but don't make too big
	MaxSz = min(MaxSz, [1000 680]);
	RescaleFactor = max(floor(min(MaxSz ./ ImSz)), 1);
end
figPosNew = ceil(ImSz * RescaleFactor);
% Don't move the figure if the size isn't changing
figPosCur = get(fig, 'Position');
if ~isequal(figPosCur(3:4), figPosNew)
    % Keep the centre of the figure stationary
    figPosNew = [max(1, floor(figPosCur(1:2)+(figPosCur(3:4)-figPosNew)/2)) figPosNew];
    % Ensure the figure bar is in bounds
    figPosNew(1:2) = min(figPosNew(1:2), MonSz(1:2)+MonSz(3:4)-[6 101]-figPosNew(3:4));
    set(fig, 'Position', figPosNew);
end

colormap(map);			% Set the colourmap

% Set the first lot of images
index = mod(0:num-1, state.n) + 1;
hw = 1 ./ layout;
gap = gap ./ layout;
dims = hw - 2 * gap;
dims = dims([2 1]);
for a = 1:layout(1)
    for b = 1:layout(2)
        c = index(b + (layout(1) - a) * layout(2));
        A = I{c};
        if ischar(A)
            A = imread_rgb(A);
            I{c} = A;
        end
        subplot('Position', [(b-1)*hw(2)+gap(2) (a-1)*hw(1)+gap(1) dims]);
        if isempty(A)
            state.h(a,b) = imagesc(zeros(1, 1, 3), lims);
            axis image off;
            set(state.h(a,b), 'CData', []);
        else
            state.h(a,b) = imagesc(A, lims);
            axis image off;
        end
    end
end
drawnow;
figure(fig); % Bring the montage into view

% Check if we need to be able to scroll through images
if n > num
    % Intialize rest of data structure
    state.index = 1;
    state.layout = layout;
    state.I = I;
    % Set the callback for image navigation, and save the image data in the figure
    set(fig, 'KeyPressFcn', @keypress_callback, 'Interruptible', 'off', 'UserData', state);
end

% ------------------------------------------------------------------------
function keypress_callback(fig, event_data)
% The function which does all the display stuff
% Check what key was pressed and update the image index as necessary
switch event_data.Character
    case 28,	up = -1;	% Left - Back a page
    case 29,	up = 1;		% Right - Forward a page
    case 30,	up = -0.1;	% Up - Back a row
    case 31,	up = 0.1;	% Down - Forward a row
    otherwise				% Another key was pressed - ignore it
		return
end
% Use control and shift for faster scrolling
if ~isempty(event_data.Modifier)
    up = up * (2 ^ (strcmpi(event_data.Modifier, {'shift', 'control'}) * [1; 2]));
end
state = get(fig, 'UserData');		% Get the state data, if not given
index = state.index;				% Get the current index
n = prod(state.layout);				% Get number of images
% Generate 12 valid indices
if (abs(up) < 1)					% Increment by row
	index = index + state.layout(2) * (up * 10) - 1;
else
	if (state.layout(1) == 1),		index = index + up - 1;		% Increment by column
	else							index = index + n * up - 1;	% Increment by page
	end
end
index = mod(index:index+n, state.n) + 1;
% Plot the images
figure(fig);
for a = 1:state.layout(1)
    for b = 1:state.layout(2)
        c = index(b + (state.layout(1) - a) * state.layout(2));
        A = state.I{c};
        if ischar(A)
            A = imread_rgb(A);
            state.I{c} = A;
        end
        set(state.h(a,b), 'CData', A);
    end
end
drawnow;
% Save the current index
state.index = index(1);
set(fig, 'UserData', state);

% ------------------------------------------------------------------------
function layout = choose_layout(n, y, x, layout)
%% Choose a good layout for the images
layout = reshape(layout, 1, min(numel(layout), 2));
v = numel(layout);
N = isnan(layout);
if v == 0 || all(N)
    sz = get(0, 'ScreenSize');
    sz = sz(3:4) ./ [x y];
    layout = ceil(sz([2 1]) ./ sqrt(prod(sz) / n));
    switch ([prod(layout - [1 0]) prod(layout - [0 1])] >= n) * [2; 1]
        case 0
        case 1
            layout = layout - [0 1];
        case 2
            layout = layout - [1 0];
        case 3
            if min(sz .* (layout - [0 1])) > min(sz .* (layout - [1 0]))
                layout = layout - [0 1];
            else
                layout = layout - [1 0];
            end
    end
elseif v == 1
    layout = layout([1 1]);
elseif any(N)
    layout(N) = ceil(n / layout(~N));
end

% ------------------------------------------------------------------------
function A = imread_rgb(name)
%% Read image to uint8 rgb array
try
    [A map] = imread(name);
catch
    % Format not recognized by imread, so create a red cross (along diagonals)
    A = eye(101) | diag(ones(100, 1), 1) | diag(ones(100, 1), -1);
    A = uint8(255 * (1 - (A | flipud(A))));
    A = cat(3, zeros(size(A), 'uint8')+uint8(255), A, A);
    return
end
if ~isempty(map)
    map = uint8(map * 256 - 0.5);
    A = reshape(map(A,:), [size(A) size(map, 2)]);
elseif size(A, 3) == 4
    ll = lower(name(end));
    if (ll == 'f')			% TIFF in CMYK colourspace - convert to RGB
        error('CMYK image files not yet supported - please fix.');
    elseif (ll == 's')		% RAS in RGBA colourspace - convert to RGB
        error('RGBA image files not yet supported - please fix.');
    end
end

% ------------------------------------------------------------------------
function L = get_im_names
%% Get the names of all images in a directory
D = dir;
n = 0;
L = cell(size(D));
% Go through the directory list
for a = 1:numel(D)
    % Check if file is a supported image type
    if numel(D(a).name) > 4 && ~D(a).isdir && (any(strcmpi(D(a).name(end-3:end), {'.png', '.tif', '.jpg', '.bmp', '.ppm', '.pgm', '.pbm', '.gif', '.ras'})) || any(strcmpi(D(a).name(end-4:end), {'.tiff', '.jpeg'})))
        n = n + 1;
        L{n} = D(a).name;
    end
end
L = L(1:n);

% ------------------------------------------------------------------------
function [map, layout, gap, indices, lims] = parse_inputs(inputs)

% Set defaults
map = gray(256);
layout = [];
gap = 0.01;
indices = -1;
lims = 0;

% Check for map
if numel(inputs) && isnumeric(inputs{1}) && size(inputs{1}, 2) == 3
    map = inputs{1};
    inputs = inputs(2:end);
end

% Go through option pairs
for a = 1:2:numel(inputs)
    switch lower(inputs{a})
        case 'map'
            map = inputs{a+1};
            if (ischar(map)),	map = eval([map '(256)']);	end
        case {'size', 'grid'}
            layout = inputs{a+1};
        case {'gap', 'border'}
            gap = inputs{a+1};
        case 'indices'
            indices = inputs{a+1};
        case {'lims', 'displayrange'}
            lims = inputs{a+1};
        otherwise
            error('Input option %s not recognized', inputs{a});
    end
end

% ------------------------------------------------------------------------
function lims = default_limits(A)
%% Return default limits for the image type
if size(A, 3) == 1
    lims = [min(reshape(A, numel(A), 1)) max(reshape(A, numel(A), 1))];
else
    lims = [0 1];
    if ~(isa(A,'single') || isa(A,'double'))
        lims = lims * double(loc_intmax(class(A)));
    end
end

% ------------------------------------------------------------------------
function imax = loc_intmax(classname)
%INTMAX Largest positive integer value.
%   LOC_INTMAX(CLASSNAME) is the largest positive value in the integer class
%   CLASSNAME. Valid values of CLASSNAME are 'int8', 'uint8', 'int16',
%   'uint16', 'int32', 'uint32'.

	switch (classname)
		case 'int8'
			imax = int8(127);
		case 'uint8'
			imax = uint8(255);
		case 'int16'
			imax = int16(32767);
		case 'uint16'
			imax = uint16(65535);
		case 'int32'
			imax = int32(2147483647);
		case 'uint32'
			imax = uint32(4294967295);
		otherwise
			error('montage:loc_intmax:invalidClassName','Invalid class name.')
	end


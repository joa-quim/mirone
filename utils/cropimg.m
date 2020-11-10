function varargout = cropimg(varargin)
%cropimg Crop image.
%   cropimg crops an image to a specified rectangle. In the syntaxes below, cropimg
%   displays the input image and waits for you to specify the crop rectangle with the mouse:
%
%        I2 = cropimg(I)
%        X2 = cropimg(X,MAP)
%        RGB2 = cropimg(RGB)
%
%   If you omit the input arguments, cropimg operates on the image in the current axes.
%
%   To specify the rectangle:
%
%   - For a single-button mouse, press the mouse button and drag to
%     define the crop rectangle. Finish by releasing the mouse button.
%
%   - For a 2- or 3-button mouse, press the left mouse button and drag to
%     define the crop rectangle. Finish by releasing the mouse button.
%
%   If you hold down the SHIFT key while dragging, or if you
%   press the right mouse button on a 2- or 3-button mouse,
%   cropimg constrains the bounding rectangle to be a square.
%
%   When you release the mouse button, cropimg returns the cropped
%   image in the supplied output argument. If you do not supply an
%   output argument, cropimg displays the output image in a new figure.
%
%   You can also specify the cropping rectangle noninteractively, using these syntaxes:
%
%        I2 = cropimg(I,RECT)
%        X2 = cropimg(X,MAP,RECT)
%        RGB2 = cropimg(RGB,RECT)
%
%   RECT is a 4-element vector with the form [XMIN YMIN WIDTH
%   HEIGHT]; these values are specified in spatial coordinates.
%
%   To specify a nondefault spatial coordinate system for the
%   input image, precede the other input arguments with two
%   2-element vectors specifying the XData and YData:
%
%        [...] = cropimg(X,Y,...)
%
%   If you supply additional output arguments, cropimg returns
%   information about the selected rectangle and the coordinate
%   system of the input image:
%
%        [A,RECT] = cropimg(...)
%        [X,Y,A,RECT] = cropimg(...)
%
%   A is the output image. X and Y are the XData and YData of the input image.
%
%        rect_ind = cropimg(X,Y,I,RECT,OPT)
%   returns a vector with [start_row end_row start_col end_col]
%   OPT must be either: "out" or "out_ind" to indicate that row and col are integers
%                       corresponding image rows and columns.
%                       "out_precise" to indicate that row and col are real numbers
%                       consisting on decimal image coordinates on the n_row,n_col system
%                       (This is used for (geo)registering images.)
%        [A,rect_ind] = cropimg(X,Y,I,RECT,'out_grid')
%                       Outputs the croped grid as well as the start and end row and column indices
%                       in "rect_ind". That is rect_ind = [start_row end_row start_col end_col]
%
%   Remarks
%   -------
%   Because RECT is specified in terms of spatial coordinates, the WIDTH and HEIGHT
%   of RECT do not always correspond exactly with the size of the output image.
%   For example, suppose RECT is [20 20 40 30], using the default spatial coordinate
%   system. The upper left corner of the specified rectangle is the center of the pixel
%   (20,20) and the lower right corner is the center of the pixel (50,60). The resulting
%   output image is 31-by-41, not 30-by-40, because the output image includes all pixels
%   in the input that are completely or partially enclosed by the rectangle.
%
%   Class Support
%   -------------
%   A must be a real image of class logical, uint8, uint16, or double. The output
%   image B is of the same class as A. RECT is always of class double.

%   (Most) Copyright 1993-2003 The MathWorks, Inc.  

n_out = nargout;

if (~isempty(varargin) && ( strcmp(varargin{end},'out') || strcmp(varargin{end},'out_ind') ))
	out = 1;    precise = 0;    varargin(end) = [];
elseif (~isempty(varargin) && strcmp(varargin{end},'out_precise'))
	out = 1;    precise = 1;    varargin(end) = [];
elseif (~isempty(varargin) && strcmp(varargin{end},'out_grid'))
	out = 1;    precise = 0;    varargin(end) = [];
	if (n_out == 2)		% OUT == 2 will signal that the 2 outputs are [crop, [r1,r2,c1,c2]]
		out = 2;
	end
else
	out = 0;    precise = 0;
end

[x,y,a,cm,rect] = ParseInputs(varargin{:});
[m,n,o] = size(a);
xmin = min(x(:));   ymin = min(y(:));
xmax = max(x(:));   ymax = max(y(:));
% Transform rectangle into row and column indices.
if (m == 1)
    pixelsPerVerticalUnit = 1;
else
    pixelsPerVerticalUnit = (m - 1) / (ymax - ymin);
end
if (n == 1)
    pixelsPerHorizUnit = 1;
else
    pixelsPerHorizUnit = (n - 1) / (xmax - xmin);
end
pixelHeight = rect(4) * pixelsPerVerticalUnit;
pixelWidth = rect(3) * pixelsPerHorizUnit;
r1 = (rect(2) - ymin) * pixelsPerVerticalUnit + 1;
c1 = (rect(1) - xmin) * pixelsPerHorizUnit + 1;
r2 = r1 + pixelHeight;
c2 = c1 + pixelWidth;
if ~precise         % row & col indexes
    r2 = round(r2);     c2 = round(c2);
    r1 = round(r1);     c1 = round(c1);
end
% Check for selected rectangle completely outside the image
if ((r1 > m) || (r2 < 1) || (c1 > n) || (c2 < 1))
    b = [];
elseif (out == 0 || out == 2)
    r1 = max(r1, 1);    r2 = min(r2, m);
    c1 = max(c1, 1);    c2 = min(c2, n);
    b = a(r1:r2, c1:c2, :);
else
    r1 = max(r1, 1);    r2 = min(r2, m);
    c1 = max(c1, 1);    c2 = min(c2, n);
    n_out = 5;          rect = [r1 r2 c1 c2];
end

switch n_out 
    case 0
        if (isempty(b))
            errordlg('The crop rectangle does not intersect the image','Error');
        end
    case 1
        varargout{1} = b;
    case 2
        if ~(out)               % varargout = [B,RECT]
            varargout{1} = b;    varargout{2} = rect;
        else                    % varargout = [B,[r1 r2 c1 c2]]
            varargout{1} = b;    varargout{2} = [r1 r2 c1 c2];
        end
    case 4
        varargout{1} = x;    varargout{2} = y;
        varargout{3} = b;    varargout{4} = rect;
    case 5
        varargout{1} = rect;
    otherwise
        errordlg('Too many output arguments','Error');
end

% --------------------------------------------------------------------------------------------------
function [x,y,a,cm,rect] = ParseInputs(varargin)
x = [];     y = [];     a = [];     flag = 0;   cm = [];    rect = [];

switch nargin
case 0
    % Cropimg. Get information from current figure
    [x,y,a,flag] = getimage;
    if (flag == 0)
        errordlg('No image found in the current axes','Error');
    end
    if (flag == 1)        % input image is indexed; get its colormap;
        cm = colormap;
    end
    %rect = getrect(gcf);
    [p1,p2]=rubberbandbox;    rect = [p1(1) p1(2) p2(1)-p1(1) p2(2)-p1(2)];
case 1                   % cropimg(I) or cropimg(RGB)
    a = varargin{1};
    x = [1 size(a,2)];    y = [1 size(a,1)];
    CheckCData(a);
    %rect = getrect(gcf);
    [p1,p2]=rubberbandbox;    rect = [p1(1) p1(2) p2(1)-p1(1) p2(2)-p1(2)];
case 2                   % cropimg(X,MAP) or cropimg(I,RECT) or cropimg(RGB,RECT)
    a = varargin{1};
    x = [1 size(a,2)];    y = [1 size(a,1)];
    if (size(varargin{2},2) == 3)        % cropimg(X,MAP)
        cm = varargin{2};
        CheckCData(a);
        %rect = getrect(gcf);
        [p1,p2]=rubberbandbox;    rect = [p1(1) p1(2) p2(1)-p1(1) p2(2)-p1(2)];
    else
        rect = varargin{2};
    end
case 3                    % cropimg(x,y,RGB) or cropimg(X,MAP,RECT)
    if (size(varargin{3},3) == 3)        % cropimg(x,y,RGB)
        x = varargin{1};        y = varargin{2};
        a = varargin{3};
        CheckCData(a);
        %rect = getrect(gcf);
        [p1,p2]=rubberbandbox;    rect = [p1(1) p1(2) p2(1)-p1(1) p2(2)-p1(2)];
    else                  % cropimg(X,MAP,RECT)
        a = varargin{1};        cm = varargin{2};
        rect = varargin{3};
        x = [1 size(a,2)];        y = [1 size(a,1)];
    end
case 4                    % cropimg(x,y,I,RECT) or cropimg(x,y,RGB,RECT)
    x = varargin{1};    y = varargin{2};
    a = varargin{3};
    rect = varargin{4};
case 5                    % cropimg(x,y,X,MAP,RECT)
    x = varargin{1};    y = varargin{2};
    a = varargin{3};    cm = varargin{4};
    rect = varargin{5};
end

% In some cases CheckCData gets called twice.  This could be avoided with more
% complex logic, but the check is quick and the code is simpler this way.  -sle
CheckCData(a);

% --------------------------------------------------------------------------------------------------
function CheckCData(cdata)
right_type = (isnumeric(cdata) | islogical(cdata)) & isreal(cdata) & ~issparse(cdata);

is_2d  = ndims(cdata) == 2;
is_rgb = (ndims(cdata) == 3) & (size(cdata,3) == 3);

if ~right_type || ~(is_2d || is_rgb)
    errordlg('CROPIMG: Invalid input image.', 'Error');
end

% --------------------------------------------------------------------------------------------------
function varargout = getimage(varargin)
%GETIMAGE Get image data from axes.
%   A = GETIMAGE(H) returns the first image data contained in
%   the Handle Graphics object H.  H can be a figure, axes,
%   image, or texture-mapped surface.  A is identical to the
%   image CData; it contains the same values and is of the same
%   class (uint8, uint16, double, etc.) as the image CData. 
%   If H is not an image or does not contain an image or 
%   texture-mapped surface, A is empty.
%
%   [X,Y,A] = GETIMAGE(H) returns the image XData in X and the
%   YData in Y. XData and YData are two-element vectors that
%   indicate the range of the x-axis and y-axis.
%
%   [...,A,FLAG] = GETIMAGE(H) returns an integer flag that
%   indicates the type of image H contains. FLAG is one of these
%   values:
%   
%       0   not an image; A is returned as an empty matrix
%       1   indexed image
%       2   intensity image with values in standard range ([0,1]
%           for double arrays, [0,255] for uint8 arrays,
%           [0,65535] for uint16 arrays)
%       3   intensity data, but not in standard range
%       4   RGB image
%
%   [...] = GETIMAGE returns information for the current axes. It
%   is equivalent to [...] = GETIMAGE(GCA).
%
%   Class Support
%   -------------
%   The output array A is of the same class as the image
%   CData. All other inputs and outputs are of class double.
%

%   Copyright 1993-2003 The MathWorks, Inc.  

him = findim(varargin{:});

[x,y,A,state] = get_image_info(him);

switch nargout
case 0                  % GETIMAGE(...)
    varargout{1} = A;
case 1                  % A = GETIMAGE(...)
    varargout{1} = A;
case 2                  % [A,FLAG] = GETIMAGE(...)
    varargout{1} = A;
    varargout{2} = state;
case 3                  % [x,y,A] = GETIMAGE(...)
    varargout{1} = x;    varargout{2} = y;
    varargout{3} = A;
case 4                  % [x,y,A,FLAG] = GETIMAGE(...)
    varargout{1} = x;    varargout{2} = y;
    varargout{3} = A;
    varargout{4} = state;
otherwise
    eid = sprintf('Images:%s:tooManyOutputArgs',mfilename);
    msg = 'Too many output arguments.';
    error(eid,msg);
end

%----------------------------------------------------------------------
% Local Function: FINDIM
%----------------------------------------------------------------------
function him = findim(varargin)
%FINDIM Find image object.
%   HIM = FINDIM(H) searches for a valid Handle Graphics Image
%   object starting from the handle H and returns its handle in
%   HIM. H may be the handle of a Figure, Axes, or Image object.
%
%   If H is an Image object, FINDIM returns it.
%
%   If H is an Axes object, FINDIM searches H for Image objects.
%   If more than one Image object is found in the Axes, FINDIM
%   looks to see if one of the Images is the current object. If
%   so, FINDIM returns that Image. Otherwise, FINDIM returns the
%   highest Image in the stacking order.
%
%   If H is a Figure object, FINDIM searches H's current Axes.
%
%   HIM = FINDIM searchs the current Figure.

him = [];
if (nargin == 0)
    rootKids = get(0,'Children');
    if (~isempty(rootKids))
        figHandle = get(0,'CurrentFigure');
        figAxes = findobj(get(figHandle, 'Children'), 'flat', 'Type', 'axes');
        if (~isempty(figAxes))
            axHandle = get(figHandle, 'CurrentAxes');
            him = findim_in_axes(axHandle);
        end
    end
else                    % User specified a handle.
    h = varargin{1};    h = h(1);
    if (~ishandle(h))
        eid = sprintf('Images:%s:invalidHandleH',mfilename);
        msg = sprintf('%s: Invalid handle H.',upper(mfilename));
        error(eid,msg);
    end
    switch get(varargin{1},'Type')
    case 'figure'
        figHandle = varargin{1};
        figAxes = findobj(get(figHandle, 'Children'), 'flat', 'Type', 'axes');
        if (~isempty(figAxes))
            axHandle = get(figHandle, 'CurrentAxes');
            him = findim_in_axes(axHandle);
        end
    case 'axes'
        axHandle = varargin{1};
        him = findim_in_axes(axHandle);
    case 'image'
        him = h;
    otherwise
        eid = sprintf('Images:%s:handleHMustBeFigAxesOrImage',mfilename);
        msg = sprintf('%s: Input handle H must be a figure, axes, or image.',upper(mfilename));
        error(eid,msg);
    end
end

%----------------------------------------------------------------------
% Local Function: FINDIM_IN_AXES
%----------------------------------------------------------------------
function him = findim_in_axes(axHandle)

figHandle = get(axHandle, 'Parent');
% If the current object is a texture-mapped surface, use that.
currentObj = get(figHandle, 'CurrentObject');
if (~isempty(currentObj) && strcmp(get(currentObj,'type'),'surface') && ...
            strcmp(get(currentObj,'FaceColor'),'texturemap'))
    him = currentObj;
else
    him = findobj(axHandle, 'Type', 'image');
    if (length(him) > 1)
        % Found more than one image in the axes. If one of the images is the current
        % object, use it. Otherwise, use the first image in the stacking order.
        if (isempty(currentObj))
            % No current object; use the one on top.
            him = him(1);
        else            % If the current object is one of the images we found, use it.
            idx = find(him == currentObj);
            if (isempty(idx))
                him = him(1);
            else
                him = him(idx);
            end
        end
    end
end
if (isempty(him))        % Didn't find an image.  Is there a texturemapped surface we can use?
    him = findobj(axHandle, 'Type', 'surface', 'FaceColor', 'texturemap');
    if (~isempty(him)),     him = him(1);    end
end
            
%----------------------------------------------------------------------
% Local Function: GET_IMAGE_INFO
%----------------------------------------------------------------------
function [x,y,A,state] = get_image_info(him)

if (isempty(him))                               % We didn't find an image.
    x = [];    y = [];    A = [];    state = 0;
elseif (strcmp(get(him, 'Type'), 'surface'))    % We found a texturemapped surface object.
    A = get(him, 'CData');
    x = get(him, 'XData');    y = get(him, 'YData');
    state = 2;
else                                            % We did find an image.  Find out about it.
    userdata = get(him, 'UserData');
    cdatamapping = get(him, 'CDataMapping');
    x = get(him, 'XData');    y = get(him, 'YData');
    A = get(him, 'CData');
    if ((ndims(A) == 3) && (size(A,3) == 3))      % We have an RGB image
        state = 4;
    else        % Not an RGB image
        if (isequal(cdatamapping,'direct'))
            % Do we have an indexed image or an old-style intensity or scaled image?
            if (isequal(size(userdata), [1 2]))
                % We have an old-style intensity or scaled image. How long is the colormap?
                N = size(get(get(get(him,'Parent'),'Parent'),'Colormap'),1);
                if (isequal(userdata, [0 1]))   % We have an old-style intensity image.
                    A = (A-1)/(N-1);
                    state = 2;
                else                            % We have an old-style scaled image.
                    A = (A-1)*((userdata(2)-userdata(1))/(N-1))+userdata(1);
                    state = 3;
                end
            else                % We have an indexed image.
                state = 1;
            end
        else            % CDataMapping is 'scaled'
            hax = get(him, 'Parent');
            clim = get(hax, 'CLim');
            if ((isa(A,'double') && isequal(clim,[0 1])) || (isa(A,'uint8') && isequal(clim,[0 255])) || ...
                  (isa(A,'uint16') && isequal(clim,[0 65535])))      % We have an intensity image.
                state = 2;
            else                % We have a scaled image.
                state = 3;
            end
        end
    end
end
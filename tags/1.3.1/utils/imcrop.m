function varargout = imcrop(varargin)
%IMCROP Crop image.
%   IMCROP crops an image to a specified rectangle. In the
%   syntaxes below, IMCROP displays the input image and waits for
%   you to specify the crop rectangle with the mouse:
%
%        I2 = IMCROP(I)
%        X2 = IMCROP(X,MAP)
%        RGB2 = IMCROP(RGB)
%
%   If you omit the input arguments, IMCROP operates on the image
%   in the current axes.
%
%   To specify the rectangle:
%
%   - For a single-button mouse, press the mouse button and drag
%     to define the crop rectangle. Finish by releasing the mouse
%     button.
%
%   - For a 2- or 3-button mouse, press the left mouse button and
%     drag to define the crop rectangle. Finish by releasing the
%     mouse button.
%
%   If you hold down the SHIFT key while dragging, or if you
%   press the right mouse button on a 2- or 3-button mouse,
%   IMCROP constrains the bounding rectangle to be a square.
%
%   When you release the mouse button, IMCROP returns the cropped
%   image in the supplied output argument. If you do not supply
%   an output argument, IMCROP displays the output image in a new
%   figure.
%
%   You can also specify the cropping rectangle noninteractively,
%   using these syntaxes:
%
%        I2 = IMCROP(I,RECT)
%        X2 = IMCROP(X,MAP,RECT)
%        RGB2 = IMCROP(RGB,RECT)
%
%   RECT is a 4-element vector with the form [XMIN YMIN WIDTH
%   HEIGHT]; these values are specified in spatial coordinates.
%
%   To specify a nondefault spatial coordinate system for the
%   input image, precede the other input arguments with two
%   2-element vectors specifying the XData and YData:
%
%        [...] = IMCROP(X,Y,...)
%
%   If you supply additional output arguments, IMCROP returns
%   information about the selected rectangle and the coordinate
%   system of the input image:
%
%        [A,RECT] = IMCROP(...)
%        [X,Y,A,RECT] = IMCROP(...)
%
%   A is the output image. X and Y are the XData and YData of the
%   input image.
%
%   Remarks
%   -------
%   Because RECT is specified in terms of spatial coordinates,
%   the WIDTH and HEIGHT of RECT do not always correspond exactly
%   with the size of the output image. For example, suppose RECT
%   is [20 20 40 30], using the default spatial coordinate
%   system. The upper left corner of the specified rectangle is
%   the center of the pixel (20,20) and the lower right corner is
%   the center of the pixel (50,60). The resulting output image
%   is 31-by-41, not 30-by-40, because the output image includes
%   all pixels in the input that are completely or partially
%   enclosed by the rectangle.
%
%   Class Support
%   -------------    
%   If you specify RECT as an input argument, then the input image can be
%   logical or numeric, and must be real and nonsparse. RECT is double.
%
%   If you do not specify RECT as an input argument, then IMCROP calls
%   IMSHOW. IMSHOW expects I to be logical, uint8, uint16, int16, double, or
%   single. RGB can be uint8, int16, uint16, double, or single. X can be
%   logical, uint8, uint16, double, or single. The input image must be real
%   and nonsparse.
% 
%   If you specify the image as an input argument, then the output image has
%   the same class as the input image.
% 
%   If you don't specify the image as an input argument, i.e., you call IMCROP
%   with 0 input arguments, then the output image has the same class as the
%   input image except for the int16 or single data type.  The output image is
%   double if the input image is int16 or single.
%
%   Example
%   -------
%        I = imread('circuit.tif');
%        I2 = imcrop(I,[60 40 100 90]);
%        figure, imshow(I), figure, imshow(I2)
%
%   See also ZOOM.

%   Copyright 1993-2004 The MathWorks, Inc.
%   $Revision: 5.21.4.5 $  $Date: 2004/08/10 01:40:07 $

[x,y,a,cm,rect] = ParseInputs(varargin{:});

m = size(a,1);
n = size(a,2);
xmin = min(x(:)); ymin = min(y(:));
xmax = max(x(:)); ymax = max(y(:));
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
r2 = round(r1 + pixelHeight);
c2 = round(c1 + pixelWidth);
r1 = round(r1);
c1 = round(c1);
% Check for selected rectangle completely outside the image
if ((r1 > m) || (r2 < 1) || (c1 > n) || (c2 < 1))
    b = [];
else
    r1 = max(r1, 1);
    r2 = min(r2, m);
    c1 = max(c1, 1);
    c2 = min(c2, n);
    b = a(r1:r2, c1:c2, :);
end

switch nargout 
case 0
    if (isempty(b))
        msg = 'The crop rectangle does not intersect the image';
        wid = sprintf('Images:%s:cropRectDoesNotIntersectImage',mfilename);
        warning(wid,msg);
    end
    figure;
    if (~isempty(cm))
        imshow(b,cm)
    else
        imshow(b);
    end
    
case 1
    varargout{1} = b;
    
case 2
    varargout{1} = b;
    varargout{2} = rect;
    
case 4
    varargout{1} = x;
    varargout{2} = y;
    varargout{3} = b;
    varargout{4} = rect;

otherwise
    msg = 'Too many output arguments';
    eid = sprintf('Images:%s:tooManyOutputArguments',mfilename);
    error(eid,msg);
end

%%%
%%% Subfunction ParseInputs
%%%
function [x,y,a,cm,rect] = ParseInputs(varargin)

x = [];
y = [];
a = [];
cm = [];
rect = [];

iptchecknargin(0,5,nargin,mfilename);

switch nargin
case 0
    % IMCROP
    % Get information from current figure
    [x,y,a,flag] = getimage;
    if (flag == 0)
        msg = 'No image found in the current axes';
        eid = sprintf('Images:%s:noImageFoundInCurrentAxes',mfilename);
        error(eid,msg);
    end
    if (flag == 1)
        % input image is indexed; get its colormap;
        cm = colormap;
    end
    rect = getrect(gcf);
    
case 1
    % IMCROP(I) or IMCROP(RGB)
    a = varargin{1};
    x = [1 size(a,2)];
    y = [1 size(a,1)];
    iptcheckinput(a,{'logical','int16','single','double','uint16', ...
                     'uint8'},{'real','nonsparse'},mfilename,'I or RGB',1);
    imshow(a);
    rect = getrect(gcf);
    
case 2
    % IMCROP(X,MAP) or IMCROP(I,RECT) or IMCROP(RGB,RECT)
    a = varargin{1};
    x = [1 size(a,2)];
    y = [1 size(a,1)];
    if numel(varargin{2}) ~= 4
        % IMCROP(X,MAP)
        cm = varargin{2};
        iptcheckinput(a,{'logical','single','double','uint16', 'uint8'}, ...
                      {'real','nonsparse'},mfilename,'X',1);
        imshow(a,cm);
        rect = getrect(gcf);
    else
        % IMCROP(I,RECT) or IMCROP(RGB,RECT)
        CheckCData(a);
        rect = varargin{2};
    end
    
case 3
    % IMCROP(x,y,RGB), IMCROP(x,y,I), or IMCROP(X,MAP,RECT)
    if (size(varargin{3},3) == 3)
        % IMCROP(x,y,RGB)
        x = varargin{1};
        y = varargin{2};
        a = varargin{3};
        iptcheckinput(a,{ 'int16','single','double','uint16', 'uint8'}, ...
                      {'real','nonsparse'},mfilename,'RGB',1);
        imshow(a,'XData',x,'YData',y);
        rect = getrect(gcf);
    elseif isvector(varargin{3}) 
        % This logic has some holes but it is less hole-ly than the previous
        % version. Furthermore, it is very unlikely that a user
        % would use IMCROP(x,y,I) if I was a vector.
        % IMCROP(X,MAP,RECT)
        a = varargin{1};
        CheckCData(a);
        cm = varargin{2};
        rect = varargin{3};
        x = [1 size(a,2)];
        y = [1 size(a,1)];
    else
        % IMCROP(x,y,I)
        x = varargin{1};
        y = varargin{2};
        a = varargin{3};
        iptcheckinput(a,{'int16','logical','single','double','uint16', 'uint8'}, ...
                      {'real','nonsparse'},mfilename,'I',1);
        imshow(a,'XData',x,'YData',y);
        rect = getrect(gcf);
    end
    
case 4
    % IMCROP(x,y,I,RECT) or IMCROP(x,y,RGB,RECT)
    x = varargin{1};
    y = varargin{2};
    a = varargin{3};
    CheckCData(a);
    rect = varargin{4};
    
case 5
    % IMCROP(x,y,X,MAP,RECT)
    x = varargin{1};
    y = varargin{2};
    a = varargin{3};
    CheckCData(a);
    cm = varargin{4};
    rect = varargin{5};
end

%%%
%%% Subfunction CheckCData
%%%
function CheckCData(cdata)

right_type = (isnumeric(cdata) || islogical(cdata)) && isreal(cdata) && ...
    ~issparse(cdata);

is_2d = ndims(cdata) == 2;
is_rgb = (ndims(cdata) == 3) && (size(cdata,3) == 3);

if ~right_type || ~(is_2d || is_rgb)
    msg = 'Invalid input image.';
    eid = sprintf('Images:%s:invalidInputImage',mfilename);
    error(eid, msg);
end

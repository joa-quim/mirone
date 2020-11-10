function pixelx = axes2pix(dim, x, axesx)
%AXES2PIX Convert axes coordinates to pixel coordinates.
%   PIXELX = AXES2PIX(DIM, X, AXESX) converts axes coordinates
%   (as returned by get(gca, 'CurrentPoint'), for example) into
%   pixel coordinates.  X should be the vector returned by
%   X = get(image_handle, 'XData') (or 'YData').  DIM is the
%   number of image columns for the x coordinate, or the number
%   of image rows for the y coordinate.

%   Copyright 1993-2002 The MathWorks, Inc.  
%   $Revision: 5.12 $  $Date: 2002/03/15 15:57:01 $

% Error checking on input arguments.
error_str = nargchk(3, 3, nargin);
if (~isempty(error_str)),   error('There must be 3 input arguments.');  end
if (max(size(dim)) ~= 1),   error('First argument must be a scalar.');  end
if (min(size(x)) > 1),      error('X must be a vector.');               end

xfirst = x(1);      xlast = x(max(size(x)));

if (dim == 1),      pixelx = axesx - xfirst + 1;    return;     end
xslope = (dim - 1) / (xlast - xfirst);
if ((xslope == 1) && (xfirst == 1)),     pixelx = axesx;
else    pixelx = xslope * (axesx - xfirst) + 1;
end

function pixval_stsbar(arg1)
%   This is based on the Matlab's pixval, but hacked in several ways.
%   For example, it shows the DEM Z value instead of it's color and it doesn't create a new figure
%
%	Coofeewright 2004 by J. Luis

if (nargin==0),  arg1 = [];    end

if ischar(arg1) % execute a callback
    switch lower(arg1)
        case 'pixvalmotionfcn'
            PixvalMotionFcn;
        case 'buttondownonimage'
            ButtonDownOnImage;
        case 'backtonormalpixvaldisplay'
            BackToNormalPixvalDisplay;
        case 'exit'
            displayBar = findobj(gcf, 'Tag', 'pixel value status bar');
            dbud = get(displayBar, 'UserData');
            uirestore_j(dbud.oldFigureUIState);
            if ishandle(displayBar),    delete(displayBar);     end
            return
    end 
else
    if (nargin == 1 && ishandle(arg1))
        figureHandle = arg1;
    else
        figureHandle = gcf;
    end
    dbud.figureHandle = figureHandle;
    dbud.displayMode = 'normal';
    % Save the interactive state of the figure. UICLEARMODE now does everything that UISUSPEND
    % used to do. It calls SCRIBECLEARMODE to take care of SCRIBE.  We WILL be restoring
    % the state - this is a change.  UICLEARMODE also registers a way to disable pixval_stsbar
    dbud.oldFigureUIState = uiclearmode_j(figureHandle,'pixval_stsbar');
    % Create position vectors for Text
    handsBar = getappdata(figureHandle,'CoordsStBar');
    tmp1 = get(handsBar(3),'YData');
    tmp2 = get(handsBar(2),'XData');
    pos = [5 tmp1(1)+2 tmp2(2)-5 abs(tmp1(2)-tmp1(1))-4];
                       
    % Create the display bar
    dbud.buttonDownImage = 0;  % Image 0 never exists      
    DisplayBar = uicontrol('Parent', figureHandle, ...
                       'Style','text', ...
                       'Units','pixels', ...
                       'Position',pos, ...
                       'Foreground', [0 0 0], ...
                       'Horiz','left', ...
                       'Tag', 'pixel value status bar', ...
                       'fontname', 'FixedWidth', ...
                       'BusyAction', 'queue', ...
                       'enable', 'inactive', ...
                       'Interruptible', 'off');
            
    % Register the motion function and button up function
    set(figureHandle, 'WindowButtonMotionFcn', 'pixval_stsbar(''PixvalMotionFcn'')')
    set(DisplayBar, 'UserData', dbud);
    PixvalMotionFcn(DisplayBar);
end

%-----------------------------------------------------------------------------------------
function PixvalMotionFcn(displayBar)
if (nargin<1),  displayBar = findobj(gcbf, 'Tag', 'pixel value status bar');   end

if isempty(displayBar)
   % Pixval is in a half-broken state.  Another function (like zoom) has
   % restored pixval callbacks and PIXVAL has already been uninstalled.
   % Call uisuspend to gracefully get rid of the callbacks.

   % Note 7/21/98 - Since we are now using UICLEARMODE instead of
   % UISUSPEND, I think that we should never get into this
   % state.  It will only happen if a user writes a function
   % which calls UIRESTORE without ever calling UICLEARMDOE or
   % SCRIBECLEARMODE.  We'll leave this code here just to be safe.
   uisuspend_j(gcbf);
   return
end

dbud = get(displayBar, 'UserData');

if strcmp(dbud.displayMode, 'normal')           % See if we're above the displayBar
      [imageHandle,imageType,img,x,y] = OverImage(dbud.figureHandle);
      if imageHandle ~= 0       
         % Update the Pixel Value display
         UpdatePixelValues(dbud.figureHandle,imageHandle, imageType, displayBar,img,x,y);
      end
elseif strcmp(dbud.displayMode, 'distance')     % If we're in distance mode and in another image, clean up a bit.
   [imageHandle,imageType,img,x,y]= OverImage(dbud.figureHandle);
   if imageHandle ~= 0 
      if imageHandle == dbud.buttonDownImage
         UpdatePixelValues(dbud.figureHandle,imageHandle, imageType, displayBar,img,x,y);
      end
   end
end
set(displayBar, 'UserData', dbud);

%-----------------------------------------------------------------------------------------
function [imageHandle,imageType,img,x,y] = OverImage(figHandle)
% Return the index of which image we are over, and return a 0 if we aren't above an image.
images = findobj(figHandle, 'type', 'image');
if isempty(images)
   imageHandle=0;   imageType='';   img=[];     x=0;    y=0;   return
end
% Make sure that the Image's Button Down & Up functions will queue
set(images, 'ButtonDownFcn', 'pixval_stsbar(''ButtonDownOnImage'')', ...
   'Interruptible', 'off', 'BusyAction', 'Queue');
axHandles = get(images, {'Parent'});
%axPositions = get([axHandles{:}], {'Position'});
axCurPt = get([axHandles{:}], {'CurrentPoint'});

% Loop over the axes, see if we are above any of them
imageHandle = 0;  
for k=1:length(axHandles)
   XLim = get(axHandles{k}, 'XLim');   YLim = get(axHandles{k}, 'YLim');
   pt = axCurPt{k};
   x = pt(1,1); y = pt(1,2);
   if x>=XLim(1) && x<=XLim(2) && y>=YLim(1) && y<=YLim(2)
      imageHandle = images(k);
      break;
   end
end

% Figure out image type
if imageHandle ~= 0
   [img,flag] = get_image_info(imageHandle);
   switch flag
   case 1
      imageType = 'indexed';
   case 2     %Grayscale in standard range
      imageType = 'intensity'; 
   case 3     %Logical or Grayscale in nonstandard range
      if islogical(img)
        imageType = 'logical';
      else
        imageType = 'intensity';
      end
   case 4
      imageType = 'rgb';
   otherwise
      msg = ['Invalid image, GETIMAGE returned flag = ' flag '.'];
      error('Images:pixval_stsbar:invalidImage', '%s', msg);
   end
else
   imageHandle=0; imageType=''; img=[]; x=0; y=0;
end

%-----------------------------------------------------------------------------------------
function UpdatePixelValues(figHandle,imageHandle, imageType, displayBar,img,x,y)
%   This is the motion callback for when we are displaying pixels.
%   Either we are in automatic display mode and the mouse pointer is
%   moving around or we are in normal mode and there has been a button-down
%   but not yet a button up. I get the current point and update the string.
Zlim = getappdata(figHandle,'Zmin_max');
if isempty(Zlim),                       no_Zlim = 1;
else                                    no_Zlim = 0;    end
if isempty(getappdata(figHandle,'dem_x')),    haveGrid = 0;
else                                    haveGrid = 1;   end

dbud = get(displayBar, 'UserData');
if (haveGrid == 1)
    if (getappdata(figHandle,'PixelMode'))
        % Inside each grid cell, which is a pixel in the screen, display only the grid node value
        [rows,cols,colors] = size(img);
        rp = axes2pix(rows, get(imageHandle,'YData'),y);
        cp = axes2pix(cols, get(imageHandle,'XData'),x);
        r = min(rows, max(1, round(rp)));   c = min(cols, max(1, round(cp)));
        Z = getappdata(figHandle,'dem_z');
        pixel(1) = double(Z(r,c));
    else
        pixel(1) = bi_linear(getappdata(figHandle,'dem_x'),getappdata(figHandle,'dem_y'),getappdata(figHandle,'dem_z'),x,y);
    end
else                    % work on a image type
    [rows,cols,colors] = size(img);
    rp = axes2pix(rows, get(imageHandle,'YData'),y);
    cp = axes2pix(cols, get(imageHandle,'XData'),x);
    r = min(rows, max(1, round(rp)));   c = min(cols, max(1, round(cp)));
    if strcmp(imageType,'indexed')
        map=get(dbud.figureHandle, 'Colormap');
        idx = img(r,c);
        if (~isa(idx,'double')), idx = double(idx)+1;   end
        idx = round(idx);
        if (idx <= size(map,1)), pixel = map(idx,:)*255;
        else                     pixel = map(end,:)*255;   end
        if (no_Zlim == 0)        % This trick is limited to 255 intervals
            dz = Zlim(2) - Zlim(1);
            pixel(1) = Zlim(1) + idx/size(map,1) * dz;
            pixel(2:3) = [];
        end
    else
        pixel = double(img(r,c,:));
    end
end

% vamos la a fazer o teste da projeccao
% ESTOIRA COM: C:\programs\gmt4\share\pslib\ISOLatin1+.ps: Too many open files
% xy_p = mapproject_m([x y],'-JG0/0/1', '-R0/360/-90/90', '-F', '-C', '-I');
% x = xy_p(1);    y = xy_p(2);
% xy_p = ogrproj([x y],'+proj=ortho +datum=WGS84');
% aaaa = ogrproj([x y],'+proj=ortho +datum=WGS84');
% x = xy_p(1);    y = xy_p(2);

% axHandle = findobj(figHandle,'Type','Axes');
% wkt = getappdata(axHandle,'ProjWKT');
% if (~isempty(wkt))
%     projStruc.SrcProjWKT = wkt;
%     xy_p = ogrproj([x y],projStruc);
%     x = xy_p(1);    y = xy_p(2);
% end

% delete testimpwkt.*; importwkt

% Find the coordinate output format
labelType = getappdata(figHandle,'LabelFormatType');
if (isempty(labelType))     % This is in fact an error test because labelType should never be empty.
    labelType = 'NotGeog';
end
% Make a copy of x & y to use in the distance case. Needed when input was transformed to a string below
x1 = x;     y1 = y;
switch labelType
    case {'DegDec', 'NotGeog'}       % This's the default. Just build the format string
        form_xy = ' %g, %g =';                
    case 'DegMin'
        out_x = degree2dms(x,'DDMM',0,'str');   x = [out_x.dd ':' out_x.mm];
        out_y = degree2dms(y,'DDMM',0,'str');   y = [out_y.dd ':' out_y.mm];
        form_xy = ' %s, %s =';
    case 'DegMinDec'
        out_x = degree2dms(x,'DDMM.x',2,'str');   x = [out_x.dd ':' out_x.mm];
        out_y = degree2dms(y,'DDMM.x',2,'str');   y = [out_y.dd ':' out_y.mm];
        form_xy = ' %s, %s =';
    case 'DegMinSec'
        out_x = degree2dms(x,'DDMMSS',0,'str');   x = [out_x.dd ':' out_x.mm ':' out_x.ss];
        out_y = degree2dms(y,'DDMMSS',0,'str');   y = [out_y.dd ':' out_y.mm ':' out_y.ss];
        form_xy = ' %s, %s =';
    case 'DegMinSecDec'
        out_x = degree2dms(x,'DDMMSS.x',2,'str');   x = [out_x.dd ':' out_x.mm ':' out_x.ss];
        out_y = degree2dms(y,'DDMMSS.x',2,'str');   y = [out_y.dd ':' out_y.mm ':' out_y.ss];
        form_xy = ' %s, %s =';
end
 
% figure out the new string
switch dbud.displayMode
case 'normal'   % Just display Z (or intensity) information
    if strcmp(imageType,'rgb') || strcmp(imageType,'indexed')
        if isa(img, 'uint8') &&  strcmp(imageType,'rgb')
            if (no_Zlim == 0 || haveGrid == 1)           % Hacked here
                pixval_str = sprintf([form_xy ' %6.3f'], x,y,pixel(1:end));
            else
                pixval_str = sprintf([form_xy ' %3d,%3d,%3d'], x,y,pixel(1:end));
            end
        elseif isa(img, 'uint16') && strcmp(imageType,'rgb')
            pixval_str = sprintf([form_xy ' %5d,%5d,%5d'], x,y,pixel(1:end));  
        elseif islogical(img) && strcmp(imageType,'rgb')
            pixval_str = sprintf([form_xy ' %1d,%1d,%1d'], x,y,pixel(1:end));  
        else  % all indexed images use double precision colormaps
            if (no_Zlim == 0 || haveGrid == 1)           % Hacked here
                pixval_str = sprintf([form_xy ' %6.3f'], x,y,pixel(1:end));
            else
                if (isequal(pixel(1),pixel(2),pixel(3)))
                    pixval_str = sprintf([form_xy ' %g'], x,y,pixel(1));
                else
                    pixval_str = sprintf([form_xy ' %g,%g,%g'], x,y,pixel(1:end));
                end
            end
        end
    else      % intensity
        if isa(img, 'uint8'),         pixval_str = sprintf([form_xy ' %3d'],x,y,pixel(1));
        elseif isa(img, 'uint16'),    pixval_str = sprintf([form_xy ' %5d'],x,y,pixel(1));
        elseif islogical(img),        pixval_str = sprintf([form_xy ' %1d'],x,y,pixel(1));
        else                          pixval_str = sprintf([form_xy ' %6.4f'],x,y,pixel(1));
        end
    end
   
case 'distance'
    handles = guidata(figHandle);
    delta_x = (x1 - dbud.x0);   delta_y = (y1 - dbud.y0);
    set(dbud.line, 'XData', [dbud.x0 x1], 'YData', [dbud.y0 y1]);
    if (handles.geog)
        switch handles.DefineMeasureUnit     % I have to do it here to allow midtime changes in preferences
            case 'n'        % Nautical miles
                scale = 1852;   str_dist = 'dist(NM)';
            case 'k'        % Kilometers
                scale = 1000;   str_dist = 'dist(km)';
            case 'm'        % Meters or user unites
                scale = 1;      str_dist = 'dist(m)';
            case 'u'        % Meters or user unites
                scale = 1;      str_dist = 'dist(usr)';
        end
        dist = vdist(dbud.y0,dbud.x0,y1,x1,handles.DefineEllipsoide) / scale;
    else
        dist = sqrt(delta_x^2 + delta_y^2);     str_dist = 'dist(usr)';
    end
    
   if strcmp(imageType,'rgb') || strcmp(imageType,'indexed')
      if (isa(img, 'uint8') &&  strcmp(imageType,'rgb') && ~haveGrid)
         pixval_str = sprintf([form_xy ' %3d,%3d,%3d  '  str_dist ' = %3.3f'], x,y,pixel(1:end),dist);
      elseif (isa(img, 'uint8') &&  strcmp(imageType,'rgb') && haveGrid)
         pixval_str = sprintf([form_xy ' %g '  str_dist ' = %4.4f'], x,y,pixel(1:end),dist);
      elseif isa(img, 'uint16') &&  strcmp(imageType,'rgb')
         pixval_str = sprintf([form_xy ' %5d,%5d,%5d  '  str_dist ' = %3.3f'], x,y,pixel(1:end),dist);
      elseif islogical(img) &&  strcmp(imageType,'rgb')
         pixval_str = sprintf([form_xy ' %1d,%1d,%1d  '  str_dist ' = %3.3f'], x,y,pixel(1:end),dist);
      else  % all indexed images use double precision colormaps
          if (no_Zlim == 0)
              pixval_str = sprintf([form_xy ' %6.3f  '  str_dist ' = %4.4f'], x,y,pixel(1:end),dist);
          else
              pixval_str = sprintf([form_xy ' %6.4f,%6.4f,%6.4f  '  str_dist ' = %3.3f'], x,y,pixel(1:end),dist);
          end
      end
   else % intensity
      if isa(img, 'uint8')
         pixval_str = sprintf([form_xy ' %3d  dist = %3.3f'], x,y,pixel(1),dist);
      elseif isa(img, 'uint16')
         pixval_str = sprintf([form_xy ' %5d  dist = %3.3f'], x,y,pixel(1),dist);
      elseif islogical(img)
         pixval_str = sprintf([form_xy ' %1d  dist = %3.3f'], x,y,pixel(1),dist);  
      else
         pixval_str = sprintf([form_xy ' %6.4f  dist = %3.3f'], x,y,pixel(1),dist);
      end
   end
end
set(displayBar, 'String', pixval_str, 'UserData', dbud);

%-----------------------------------------------------------------------------------------
function ButtonDownOnImage
imageHandle = gcbo;
figureHandle = get(get(imageHandle,'Parent'),'Parent');
displayBar = findobj(figureHandle, 'Tag', 'pixel value status bar');
dbud = get(displayBar, 'UserData');
axesHandle = get(imageHandle, 'Parent');
% Set the initial point (x0,y0)
pt = get(axesHandle, 'CurrentPoint');
dbud.x0 = pt(1,1);
dbud.y0 = pt(1,2);
dbud.line = line('Parent', axesHandle, 'erasemode', 'xor', 'color', [1 0 0], ...
   'Xdata', [dbud.x0 dbud.x0],'Ydata', [dbud.y0 dbud.y0]);
dbud.displayMode = 'distance';
dbud.buttonDownImage = imageHandle;
set(displayBar, 'UserData', dbud);
set(dbud.figureHandle, 'WindowButtonUpFcn', 'pixval_stsbar(''BackToNormalPixvalDisplay'')');
PixvalMotionFcn(displayBar);

%-----------------------------------------------------------------------------------------
function BackToNormalPixvalDisplay
displayBar = findobj(gcbo, 'Tag', 'pixel value status bar');
dbud = get(displayBar, 'UserData');
delete(dbud.line);
dbud.line = [];     dbud.x0 = []; dbud.y0 = [];
set(dbud.figureHandle, 'WindowButtonUpFcn', '');
dbud.displayMode = 'normal';
dbud.buttonDownImage = 0;
set(displayBar, 'UserData', dbud);
PixvalMotionFcn(displayBar);

%-----------------------------------------------------------------------------------------
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

if (max(size(dim)) ~= 1);   error('First argument must be a scalar.');  end
if (min(size(x)) > 1);      error('X must be a vector.');               end
xfirst = x(1);      xlast = x(max(size(x)));
if (dim == 1);      pixelx = axesx - xfirst + 1;    return;     end
xslope = (dim - 1) / (xlast - xfirst);
if ((xslope == 1) && (xfirst == 1));     pixelx = axesx;
else    pixelx = xslope * (axesx - xfirst) + 1;         end

%----------------------------------------------------------------------
function [A,state] = get_image_info(him)
if (isempty(him))                               % We didn't find an image.
    A = [];    state = 0;
elseif (strcmp(get(him, 'Type'), 'surface'))    % We found a texturemapped surface object.
    A = get(him, 'CData');    state = 2;
else                                            % We did find an image.  Find out about it.
    userdata = get(him, 'UserData');
    cdatamapping = get(him, 'CDataMapping');
    A = get(him, 'CData');
    if ((ndims(A) == 3) && (size(A,3) == 3))     % We have an RGB image
        state = 4;
    else                                        % Not an RGB image
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
            else                                % We have an indexed image.
                state = 1;
            end
        else                                    % CDataMapping is 'scaled'
            hax = get(him, 'Parent');
            clim = get(hax, 'CLim');
            if ((isa(A,'double') && isequal(clim,[0 1])) || ...
                  (isa(A,'uint8') && isequal(clim,[0 255])) || ...
                  (isa(A,'uint16') && isequal(clim,[0 65535])))
                % We have an intensity image.
                state = 2;
            else                                % We have a scaled image.
                state = 3;
            end
        end
    end
end

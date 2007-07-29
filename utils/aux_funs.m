%# mex
function  varargout = aux_funs(opt,varargin)
% This contains Mirone's auxiliay functions that are called by the several
% of the Mirone's callback functions. I puted them here to release somehow
% the burden of the non-stop groing length of the Mirone code.

switch opt
    case 'StoreZ'
        StoreZ(varargin{:})
    case 'msg_dlg'
        varargout = msg_dlg(varargin{:});
    case 'in_map_region'
        [varargout{1} varargout{2} varargout{3} varargout{4}] = in_map_region(varargin{:});
    case 'cleanGRDappdata'
        clean_GRDappdata(varargin{:})
    case 'guessGeog'
        varargout{1} = guessGeog(varargin{:});
    case 'colormap_bg'
        colormap_bg(varargin{:})
    case 'findFileType'
        varargout{1} = findFileType(varargin{:});
    case 'isProj'
        varargout{1} = isProj(varargin{:});
    case 'toProjPT'
        toProjPT(varargin{:});
    case 'polysplit'
        [varargout{1} varargout{2}] = localPolysplit(varargin{:});
    case 'adjust_lims'
        [varargout{1} varargout{2}] = adjust_lims(varargin{:});
    case 'axes2pix'
        varargout{1} = axes2pix(varargin{:});
    case 'insideRect'
        varargout{1} = insideRect(varargin{:});
    case 'strip_bg_color'
        varargout{1} = strip_bg_color(varargin{:});
    case 'help'
        str = [varargin{1}.home_dir filesep 'doc' filesep 'MironeMan.chm'];
        if (ispc),      dos(str);
        else            unix(str);
        end
    case 'min_max_single'
        [varargout{1} varargout{2}] = min_max_single(varargin{:});
end

% --------------------------------------------------------------------
function StoreZ(handles,X,Y,Z)
    % If grid size is not to big I'll store it
    if (numel(Z)*4 > handles.grdMaxSize),       return;     end
    if (~isa(Z,'single')),  setappdata(handles.figure1,'dem_z',single(Z));
    else                    setappdata(handles.figure1,'dem_z',Z);  end
    setappdata(handles.figure1,'dem_x',X);  setappdata(handles.figure1,'dem_y',Y);

% --------------------------------------------------------------------
function [x,y,indx,indy] = in_map_region(handles,x,y,tol,map_lims)
	%   Given X & Y vectors retain only the elements that are inside the current map region
	%   OPTIONS:
	%   TOL is used normally when ploting lines and serves to extend the map
	%       limits so that lines are allowed to be drawn until the image borders
	%   MAP_LIMS a 1x4 vector with [x_min x_max y_min y_max]. If not given it will be geted here
	%   NOTE THAT ALL ARGUMENTS MUST BE PROVIDED, EVEN IF THEY ARE EMPTY
	if (isempty(tol)),   tol = 0;    end
	if (isempty(map_lims))
        x_lim = get(handles.axes1,'Xlim');      y_lim = get(handles.axes1,'Ylim');
	else
        x_lim(1:2) = map_lims(1:2);    y_lim(1:2) = map_lims(3:4);
	end
	indx = find((x < x_lim(1)-tol) | (x > x_lim(2)+tol));
	x(indx) = [];           y(indx) = [];
	if (nargout == 2),      clear indx;     end         % Save memory
	indy = find((y < y_lim(1)-tol) | (y > y_lim(2)+tol));
	x(indy) = [];           y(indy) = [];
	axes(handles.axes1)     % This is for the GCP mode be able to plot on the Master Image

% --------------------------------------------------------------------
function colormap_bg(handles,Z,pal)
% Insert the background color in the palette for arrays that have NaNs
% [m,n] = size(Z);    dbl_size = m*n*8;
% if (dbl_size > handles.grdMaxSize)  % I'm almost shure that isnan uses doubles
%     n_stripes = round(dbl_size / handles.grdMaxSize) * 2;
%     d_stripe = fix(m / n_stripes);
%     is_it = 0;      tmp = 0;
%     for (i=1:n_stripes-1)
%         is = (i-1)*d_stripe + 1;    ie = is + d_stripe;
%         tmp = any(isnan(Z(is:ie,:)));
%         tmp = any(tmp);     % Make it a single value
%         is_it = is_it + tmp;
%     end
%     tmp = isnan(Z(ie+1:end,:));     tmp = any(tmp);
%     is_it = is_it + tmp;
%     if (is_it)              pal = [handles.bg_color; pal];   end
% else
%     if any(isnan(Z(:)))     pal = [handles.bg_color; pal];   end
% end
% tic;    is_it = any(Z(:)~=Z(:)); toc
% tic;    is_it = any(isnan(Z(:))); toc
% return

%if any(Z(:)~=Z(:))     pal = [handles.bg_color; pal];   end
if ( handles.have_nans & ~isequal(pal(1,:),handles.bg_color) )
    if (size(pal,1) == 256),    pal = [handles.bg_color; pal(2:end,:)];     % Jump firts color to not have
    else                        pal = [handles.bg_color; pal];              % a CMAP with more than 256 colors
    end
end
set(handles.figure1,'Colormap',pal)

% ----------------------------------------------------------------------------------
function out = findFileType(fname)
    % From the extension guess what function should be called to open this file
	out = [];	[PATH,FNAME,EXT] = fileparts(fname);
	if ( any(strcmpi(EXT,{'.grd' '.nc'})) )
        out = 'gmt';        return
	end
	if ( any(strcmpi(EXT,{'.jpg' '.png' '.bmp' '.gif' '.pcx' '.ras' '.ppm' '.pgm' '.xwd' '.shade' '.raw' '.bin'})) )
        out = 'generic';    return
	end
	if ( any(strcmpi(EXT,{'.tif' '.tiff' '.sid' '.ecw' '.jp2'})) )
        out = 'geotif';    return
	end
	if ( any(strcmpi(EXT,{'.n1' '.n14' '.n15' '.n16' '.n17'})) )
        out = 'multiband';    return
	end
	if ( any(strcmpi(EXT,'.img')) ),    out = 'envherd';    end

% --------------------------------------------------------------------
function out = msg_dlg(in,handles)
% Cast a window with a pre-deffined message
out = {0};    msg = [];     h = [];
if (nargin == 1),   handles.no_file = 0;    handles.image_type = 0;    end      % backawrd comp
if (handles.no_file == 1)
    msg = 'You didn''t even load a file. What are you expecting then?';    out = {1};
    h = msgbox(msg,'Error');      movegui(h,'center');      return;
end
switch in
    case 1,     msg = 'C''mon load a file first. It''s logic, isn''t it?';    out = {1};
    case 2,     msg = 'This operation is deffined only for images derived from DEM grids';  out = {1};
    case 3
        if (~handles.geog)
            msg = 'This operation is currently possible only for geographic type data';
            out = {1};
        end
    case 14
        if ~(handles.image_type == 1 || handles.image_type == 4 || handles.computed_grid)
            msg = 'This operation is deffined only for images derived from DEM grids';
            out = {1};
        end
    case 5
        if (~handles.is_projected && ~handles.geog)
            msg = 'This operation is only possible for geographic data OR when the Map Projection is known';
            out = {1};            
        end
end
if (~isempty(msg))
    h = msgbox(msg,'Error');      movegui(h,'center')
end
axes(handles.axes1)     % This is for the GCP mode be able to plot on the Master Image
if (~isempty(h))        % Bring the figure forward because it was hiden by the previous command
    figure(h)
end

% ----------------------------------------------------------------------------------
function handles = isProj(handles, opt)
    % Se if we have projected grid/images and act acordingly

    % Fish eventual proj strings
    projGMT = getappdata(handles.figure1,'ProjGMT');
    projWKT = getappdata(handles.figure1,'ProjWKT');
    proj4 = getappdata(handles.figure1,'Proj4');
    
    if (~handles.geog)
        if (~isempty(projWKT))              % We have a GDAL projected file
            handles.is_projected = 1;
            % Desable the Projections menu.
            if (ishandle(handles.Projections)),     set(handles.Projections,'Enable','off');   end
        elseif (~isempty(projGMT) || ~isempty(proj4))  % We have a GMT or Proj4 projection selection
            handles.is_projected = 1;
            if (ishandle(handles.Projections)),     set(handles.Projections,'Enable','on');   end
        else                                % We know nothing about these coords (e.g. an image)
            handles.is_projected = 0;
            if (ishandle(handles.Projections)),     set(handles.Projections,'Enable','on');   end
        end
    else
        handles.is_projected = 0;
        if (ishandle(handles.Projections)),     set(handles.Projections,'Enable','off');   end
    end
    if (handles.is_projected),      set(handles.hAxMenuLF, 'Vis', 'on', 'Separator','on')
    else                            set(handles.hAxMenuLF, 'Vis', 'off', 'Separator','off')
    end
    
    if (handles.is_projected)       % When read a gdal referenced file we need to set this right away
        toProjPT(handles)
        setappdata(handles.figure1,'DispInGeogs',0)     % We need to set it now (this is the default)
    end
    
    if (nargin == 2),       guidata(handles.figure1, handles);      end

% ----------------------------------------------------------------------------------
function toProjPT(handles)
    % The following is to deal with eventual coords display in geogs and is used in PIXVAL_STSBAR 

    displayBar = findobj(handles.figure1, 'Tag', 'pixValStsBar');
	dbud = get(displayBar, 'UserData');
    dbud.toProjPT = 0;
    if (handles.is_projected)
        % Fish eventual proj strings
        projGMT = getappdata(handles.figure1,'ProjGMT');
        projWKT = getappdata(handles.figure1,'ProjWKT');
        proj4 = getappdata(handles.figure1,'Proj4');
        if (~isempty(projWKT) || ~isempty(proj4))
            % In case we have both 'projWKT' takes precedence because it came from file metadata
            if (~isempty(projWKT)),     dbud.projStruc.SrcProjWKT = projWKT;
            else                        dbud.projStruc.SrcProjWKT = ogrproj(proj4);
            end
            dbud.toProjPT = 1;
        elseif (~isempty(projGMT))
            lims = [getappdata(handles.axes1,'ThisImageLims')];
            out = mapproject_m([lims(1) lims(3); lims(2) lims(4)],'-R-180/180/0/80','-I','-F',projGMT{:});    % Convert lims back to geogs
            x_min = min(out(:,1));        x_max = max(out(:,1));
            y_min = min(out(:,2));        y_max = max(out(:,2));
      	    dbud.opt_R = ['-R' sprintf('%f',x_min) '/' sprintf('%f',x_max) '/' sprintf('%f',y_min) '/' sprintf('%f',y_max)];
            dbud.projGMT = projGMT;
            dbud.toProjPT = 2;
        end
    	set(displayBar, 'UserData', dbud);
    end

% --------------------------------------------------------------------
function [z_min,z_max] = min_max_single(Z)
% Compute the min/max of single precision Z arrays. I need this due to (another) Matlab
% bug that gives wrong results when the Z (single) array has NaNs. Ouput are doubles.
z_min = double(min(Z(~isnan(Z(:)))));   z_max = double(max(Z(~isnan(Z(:)))));

% --------------------------------------------------------------------
function img = strip_bg_color(handles,img)
% Strip eventual row/columns with color equal to the Figure's background color
bg_color = uint8(get(handles.figure1,'color')*255);
c1 = bg_color(1);    c2 = bg_color(2);    c3 = bg_color(3);
center_row = round(size(img,1) / 2);
center_col = round(size(img,2) / 2);
h_Xlabel = get(handles.axes1,'Xlabel');
h_Ylabel = get(handles.axes1,'Ylabel');
% Strip north
i = 1;
while (img(i,center_col,1) == c1 && img(i,center_col,2) == c2 && img(i,center_col,3) == c3)
    i = i + 1;
end;    i = i - 1;
img(1:i,:,:) = [];      % Strip the eventual gray band
% Strip west
j = 0;
while (img(center_row,end-j,1) == c1 && img(center_row,end-j,2) == c2 && img(center_row,end-j,3) == c3)
    j = j + 1;
end;    j = j - 1;
img(:,end-j:end,:) = [];      % Strip the eventual gray band
% Strip south
i = 0;
while (img(end-i,center_col,1) == c1 && img(end-i,center_col,2) == c2 && img(end-i,center_col,3) == c3)
    i = i + 1;
end;    i = i - 1;
img(end-i:end,:,:) = [];      % Strip the eventual gray band
% Strip east
j = 1;
while (img(center_row,j,1) == c1 && img(center_row,j,2) == c2 && img(center_row,j,3) == c3)
    j = j + 1;
end;    j = j - 1;
img(:,1:j,:) = [];      % Strip the eventual gray band

% --------------------------------------------------------------------
function [X,Y] = adjust_lims(X,Y,m,n)
% Convert the image limits from pixel reg to grid reg
dx = (X(2) - X(1)) / n;         dy = (Y(2) - Y(1)) / m;
X(1) = X(1) + dx/2;             X(2) = X(2) - dx/2;
Y(1) = Y(1) + dy/2;             Y(2) = Y(2) - dy/2;

% --------------------------------------------------------------------
function geog = guessGeog(lims)
    % Make a good guess if LIMS are geographic
    geog = double( ( (lims(1) >= -180 && lims(2) <= 180) || (lims(1) >= 0 && lims(2) <= 360) )...
        && (lims(3) >= -90 && lims(4) <= 90) );

% --------------------------------------------------------------------
function res = insideRect(rect,pt)
    % Check which elements of the  [x y] (Mx2) PT array are inside the rectangle RECT
    % RECT = [x_min x_max y_min y_max]
    % RES is a logical column vector with length = size(PT,1)
    % NO ERROR TESTING
    res = ( pt(:,1) >= rect(1) & pt(:,1) <= rect(2) & pt(:,2) >= rect(3) & pt(:,2) <= rect(4) );

% --------------------------------------------------------------------
function clean_GRDappdata(handles)
try
    rmappdata(handles.figure1,'dem_x');     rmappdata(handles.figure1,'dem_y');     rmappdata(handles.figure1,'dem_z');
    rmappdata(handles.figure1,'GMThead');
end

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

% --------------------------------------------------------------------------------
function [latcells,loncells] = localPolysplit(lat,lon)
%POLYSPLIT Extract segments of NaN-delimited polygon vectors to cell arrays
%
%   [LATCELLS,LONCELLS] = POLYSPLIT(LAT,LON) returns the NaN-delimited
%   segments of the vectors LAT and LON as N-by-1 cell arrays with one
%   polygon segment per cell.  LAT and LON must be the same size and have
%   identically-placed NaNs.  The polygon segments are column vectors if
%   LAT and LON are column vectors, and row vectors otherwise.

% Copyright 1996-2006 The MathWorks, Inc.
% $Revision: 1.4.4.5 $    $Date: 2006/05/24 03:35:26 $

	[lat, lon] = localRemoveExtraNanSeps(lat, lon);
	indx = find(isnan(lat(:)));         % Find NaN locations.
	
	% Simulate the trailing NaN if it's missing.
	if ~isempty(lat) && ~isnan(lat(end))
        indx(end+1,1) = numel(lat) + 1;
	end
	
	%  Extract each segment into pre-allocated N-by-1 cell arrays, where N is
	%  the number of polygon segments.  (Add a leading zero to the indx array
	%  to make indexing work for the first segment.)
	N = numel(indx);
	latcells = cell(N,1);       loncells = cell(N,1);
	indx = [0; indx];
	for k = 1:N
        iStart = indx(k)   + 1;
        iEnd   = indx(k+1) - 1;
        latcells{k} = lat(iStart:iEnd);
        loncells{k} = lon(iStart:iEnd);
	end

% --------------------------------------------------------------------------------
function [xdata, ydata, zdata] = localRemoveExtraNanSeps(xdata, ydata, zdata)
    %removeExtraNanSeps  Clean up NaN separators in polygons and lines

	p = find(isnan(xdata(:)'));     % Determing the positions of each NaN.
	
	% Determine the position of each NaN that is not the final element in a sequence of contiguous NaNs.
	q = p(diff(p) == 1);
	
	% If there's a leading sequence of NaNs (a sequence starting with a NaN in
	% position 1), determine the position of each NaN in this sequence.
	if isempty(p),      r = [];
	else                r = find((p - (1:numel(p))) == 0);
	end
	
	% Determine the position of each excess NaN.
	if isempty(r),      s = q;
	else                s = [r q(q > r(end))];
	end
	
	% Remove the excess NaNs.
	xdata(s) = [];      ydata(s) = [];
	if (nargin >= 3),   zdata(s) = [];  end

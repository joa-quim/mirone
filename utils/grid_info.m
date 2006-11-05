function grid_info(handles,X,Y,hdr,Z)
%#function grdinfo_m
global  home_dir;

if (nargin == 3 && strcmp(Y,'gdal'))            % Just extact the relevant info from the attribute struct
    att2Hdr(handles,X);    return
elseif (nargin == 3 && strcmp(Y,'iminfo'))
    img2Hdr(handles,X);    return
elseif (nargin == 4 && strcmp(Y,'iminfo'))      % Used with internaly generated images
    img2Hdr(handles,X,hdr);    return
end

ud = get(handles.figure1,'UserData');           % Retrieve Image info 
if (handles.image_type == 1 && ~handles.computed_grid)          % Image derived from a grdfile
    if (handles.grdformat > 1), handles.grdname = [handles.grdname '=' num2str(handles.grdformat)];     end
    info1 = grdinfo_m(handles.grdname,'hdr_struct');    % info1 is a struct with the GMT grdinfo style
    info2 = grdutils(Z,'-H');                           % info2 is a vector with [z_min z_max i_zmin i_zmax n_nans mean std]
    w{1} = ['Title: ' info1.Title];
    w{2} = ['Command: ' info1.Command];
    w{3} = ['Remark: ' info1.Remark];
    w{4} = info1.Registration;
    w{5} = ['grdfile format: #' num2str(info1.Scale(3))];
    txt1 = num2str(hdr(1),'%.7f');      txt1 = wipe_zeros(txt1);    % x_min
    txt2 = num2str(hdr(2),'%.7f');      txt2 = wipe_zeros(txt2);    % x_max
    txt3 = num2str(hdr(8),'%.7f');      txt3 = wipe_zeros(txt3);    % x_inc
    w{6} = ['x_min: ' txt1 '  x_max: ' txt2 '  x_inc: ' txt3 '  nx: ' num2str(info1.X_info(4))];
    txt1 = num2str(hdr(3),'%.7f');      txt1 = wipe_zeros(txt1);    % y_min
    txt2 = num2str(hdr(4),'%.7f');      txt2 = wipe_zeros(txt2);    % y_max
    txt3 = num2str(hdr(8),'%.7f');      txt3 = wipe_zeros(txt3);    % y_inc
    w{7} = ['y_min: ' txt1 '  y_max: ' txt2 '  y_inc: ' txt3 '  ny: ' num2str(info1.Y_info(4))];
    txt1 = num2str(hdr(5),'%.3f');      txt1 = wipe_zeros(txt1);    % z_min
    txt2 = num2str(hdr(6),'%.3f');      txt2 = wipe_zeros(txt2);    % z_max
    
    if (hdr(6)),    half = 0.5;
    else             half = 0;       end
    x_min = hdr(1) + (rem(info2(3), length(X)) + half) * hdr(8);    % x of z_min
    x_max = hdr(1) + (rem(info2(4), length(X)) + half) * hdr(8);    % x of z_max
    y_min = hdr(4) - (rem(info2(3), length(Y)) + half) * hdr(9);    % y of z_min
    y_max = hdr(4) - (rem(info2(4), length(Y)) + half) * hdr(9);    % y of z_max
    txt_x1 = num2str(x_min,'%.7f');     txt_x1 = wipe_zeros(txt_x1);
    txt_x2 = num2str(x_max,'%.7f');     txt_x2 = wipe_zeros(txt_x2);
    txt_y1 = num2str(y_min,'%.7f');     txt_y1 = wipe_zeros(txt_y1);
    txt_y2 = num2str(y_max,'%.7f');     txt_y2 = wipe_zeros(txt_y2);
    
    w{8} = ['z_min: ' txt1 ' at x = ' txt_x1 ' y = ' txt_y1 '  z_max: ' txt2 ' at x = ' ...
            txt_x2 ' y = ' txt_y2];

    w{9} = ['scale factor: ' num2str(info1.Scale(1)) ' add_offset: ' num2str(info1.Scale(2))];
    txt1 = num2str(info2(6),'%.3f');    txt1 = wipe_zeros(txt1);    % mean
    txt2 = num2str(info2(7),'%.3f');    txt2 = wipe_zeros(txt2);    % stdev
    w{10} = ['mean: ' txt1 '  stdev: ' txt2];
    if (info2(5))       % We have NaNs, report them also
        w{11} = ['nodes set to NaN: ' num2str(info2(5))];
    end
    msgbox(w,'Grid Info');
elseif (handles.computed_grid)  % Computed array
    w{1} = '    INTERNALY COMPUTED GRID';   w{2} = ' ';
    w{3} = ['   Xmin:  ' num2str(handles.head(1)) '    Xmax: ' num2str(handles.head(2))];
    w{4} = ['   Ymin:  ' num2str(handles.head(3)) '    Ymax: ' num2str(handles.head(4))];
    w{5} = ['   Zmin:  ' num2str(handles.head(5)) '    Zmax: ' num2str(handles.head(6))];
    w{6} = ['   Xinc:  ' num2str(handles.head(8)) '    Yinc: ' num2str(handles.head(9))];
    one_or_zero = ~(handles.head(7) == 1);      % To give correct nx,ny with either grid or pixel registration
    nx = round((handles.head(2) - handles.head(1))/handles.head(8) + one_or_zero);
    ny = round((handles.head(4) - handles.head(3))/abs(handles.head(9)) + one_or_zero);
    w{7} = ['   nx:  ' num2str(nx) '    ny: ' num2str(ny)];
    msgbox(w,'Grid Info');
else
    InfoMsg = getappdata(handles.axes1,'InfoMsg');
    if (~isempty(InfoMsg))
        msgbox(InfoMsg,'Image Info');
    else
        msgbox('Info missing or nothing to info about?','???')
    end
end

% -----------------------------------------------------------
function txt = wipe_zeros(txt)
	% Wipe zeros at the end of the TXT string
	while (txt(end) == '0')
        txt(end) =[];
	end

% --------------------------------------------------------------------
function Hdr = att2Hdr(handles,att)
	% Fill a header with the info from the att struct issued by gdalread

    w{1} = ['Driver : ' att.DriverShortName];
    w{2} = att.ProjectionRef;
    w{3} = [];
    w{4} = ['Width:  ' num2str(att.RasterXSize) '    Height:  ' num2str(att.RasterYSize)];
    w{5} = ['Pizel Size:  (' num2str(att.GMT_hdr(8)) ',' num2str(att.GMT_hdr(9)) ')'];
    w{6} = 'Projected corner coordinates';
    w{7}  = ['   Xmin:  ' num2str(att.Corners.LL(1)) '    Xmax: ' num2str(att.Corners.UR(1))];
    w{8}  = ['   Ymin:  ' num2str(att.Corners.LL(2)) '    Ymax: ' num2str(att.Corners.UR(2))];
    if (~isempty(att.GEOGCorners))
        w{9} = 'Geographical corner coordinates';
        w{10} = ['   Lon min:  ' att.GEOGCorners{1,1} '    Lon max: '  att.GEOGCorners{4,1}];
        w{11} = ['   Lat min:  ' att.GEOGCorners{1,2} '    Lat max: '  att.GEOGCorners{2,2}];
    end
    w{end+1} = ['   Zmin:  ' num2str(att.GMT_hdr(5)) '   Zmax: ' num2str(att.GMT_hdr(6))];
    w{end+1} = ['Color Type:  ' att.ColorInterp];
    
    setappdata(handles.axes1,'InfoMsg',w)
    
    if (~isempty(att.ProjectionRef))            % Save Proj WKT for eventual later use
        setappdata(handles.axes1,'ProjWKT',att.ProjectionRef)
        out = decodeProjectionRef(att.ProjectionRef);
        setappdata(handles.axes1,'DatumProjInfo',out)
    end

% --------------------------------------------------------------------
function out = decodeProjectionRef(strProj)
    ind = findstr(strProj,char(10));
    out.datum = [];     out.ellipsoid = [];     out.projection = [];
    if (numel(ind) <= 1),   return;    end
    
    ind = [0 ind length(strProj)-1];
    for (i=1:numel(ind)-1)
        str = strProj(ind(i)+1:ind(i+1)-1);
        if ~isempty(findstr(str,'GEOGCS'))      % Get datum
            xx = findstr(str,'"');
            if (numel(xx) < 2),     continue;   end
            out.datum = str(xx(1)+1:xx(2)-1);
        end
        if ~isempty(findstr(str,'SPHEROID'))      % Get ellipsoid
            if (numel(xx) < 2),     continue;   end
            xx = findstr(str,'"');
            out.ellipsoid = str(xx(1)+1:xx(2)-1);
        end
        if ~isempty(findstr(str,'PROJECTION'))      % Get ellipsoid
            if (numel(xx) < 2),     continue;   end
            xx = findstr(str,'"');
            out.projection = str(xx(1)+1:xx(2)-1);
        end
    end

% --------------------------------------------------------------------
function img2Hdr(handles,imgName,img)
    w = [];
    if (nargin == 2)
        try
            info_img = imfinfo(imgName);
            w{1} = ['File Name:    ' info_img.Filename];
            w{2} = ['Image Size:    ' num2str(info_img.FileSize) '  Bytes'];
            w{3} = ['Width:  ' num2str(info_img.Width) '    Height:  ' num2str(info_img.Height)];
            w{4} = ['Bit Depth:  ' num2str(info_img.BitDepth)];
            w{5} = ['Color Type:  ' info_img.ColorType];
        end
    else
        [m n k] = size(img);
        w{1} = 'File Name:  none (imported array)';
        w{2} = ['Image Size:    ' num2str(m*n*k) '  Bytes'];
        w{3} = ['Width:  ' num2str(n) '    Height:  ' num2str(m)];
        if (k == 1)
            w{4} = 'Bit Depth:  8 bits';
            w{5} = 'Color Type:  Indexed';
        else
            w{4} = 'Bit Depth:  24 bits';
            w{5} = 'Color Type:  True Color';
        end
    end
    setappdata(handles.axes1,'InfoMsg',w)

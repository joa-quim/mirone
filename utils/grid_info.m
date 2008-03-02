function grid_info(handles,X,Y,hdr)
%#function grdinfo_m

if (nargin == 3 && strcmp(Y,'gdal'))            % Just extract the relevant info from the attribute struct
    att2Hdr(handles,X);    return
elseif (nargin == 3 && strcmp(Y,'iminfo'))
    img2Hdr(handles,X);    return
elseif (nargin == 4 && strcmp(Y,'iminfo'))      % Used with internaly generated images
    img2Hdr(handles,X,hdr);    return
end

if (handles.no_file),     return;      end
%ud = get(handles.figure1,'UserData');           % Retrieve Image info 
if (handles.image_type == 1 && ~handles.computed_grid)          % Image derived from a GMT grdfile
	info1 = grdinfo_m(handles.grdname,'hdr_struct');    % info1 is a struct with the GMT grdinfo style
	X = getappdata(handles.figure1,'dem_x');    Y = getappdata(handles.figure1,'dem_y');
	Z = getappdata(handles.figure1,'dem_z');    % We want the Z for statistics (might not be in argin)
	if (~isempty(Z))
		info2 = grdutils(Z,'-H');               % info2 is a vector with [z_min z_max i_zmin i_zmax n_nans mean std]
		info2(3:4) = info2(3:4) + 1;            % Info from grdutils is zero based
	else
		info2 = zeros(7,1);     X = [0 1];      Y = [0 1];  % Just to no error hang bellow
	end
	w{1} = ['Title: ' info1.Title];
	w{2} = ['Command: ' info1.Command];
	w{3} = ['Remark: ' info1.Remark];
	w{4} = info1.Registration;
	w{5} = ['grdfile format: #' num2str(info1.Scale(3))];
	txt1 = sprintf('%.6f',handles.head(1));      txt1 = wipe_zeros(txt1);    % x_min
	txt2 = sprintf('%.6f',handles.head(2));      txt2 = wipe_zeros(txt2);    % x_max
	txt3 = sprintf('%.7f',handles.head(8));      txt3 = wipe_zeros(txt3);    % x_inc
	w{6} = ['x_min: ' txt1 '  x_max: ' txt2 '  x_inc: ' txt3 '  nx: ' num2str(info1.X_info(4))];
	txt1 = sprintf('%.6f',handles.head(3));      txt1 = wipe_zeros(txt1);    % y_min
	txt2 = sprintf('%.6f',handles.head(4));      txt2 = wipe_zeros(txt2);    % y_max
	txt3 = sprintf('%.7f',handles.head(9));      txt3 = wipe_zeros(txt3);    % y_inc
	w{7} = ['y_min: ' txt1 '  y_max: ' txt2 '  y_inc: ' txt3 '  ny: ' num2str(info1.Y_info(4))];
	txt1 = sprintf('%.4f',info2(1));            txt1 = wipe_zeros(txt1);    % z_min
	txt2 = sprintf('%.4f',info2(2));            txt2 = wipe_zeros(txt2);    % z_max

	if (handles.head(7)),   half = 0.5;
	else                    half = 0;       end
	x_min = handles.head(1) + (fix(info2(3) / length(X)) + half) * handles.head(8);    % x of z_min
	x_max = handles.head(1) + (fix(info2(4) / length(X)) + half) * handles.head(8);    % x of z_max
	y_min = handles.head(3) + (rem(info2(3)-1, length(Y)) + half) * handles.head(9);   % y of z_min
	y_max = handles.head(3) + (rem(info2(4)-1, length(Y)) + half) * handles.head(9);   % y of z_max
	txt_x1 = sprintf('%.6f',x_min);     txt_x1 = wipe_zeros(txt_x1);
	txt_x2 = sprintf('%.6f',x_max);     txt_x2 = wipe_zeros(txt_x2);
	txt_y1 = sprintf('%.6f',y_min);     txt_y1 = wipe_zeros(txt_y1);
	txt_y2 = sprintf('%.6f',y_max);     txt_y2 = wipe_zeros(txt_y2);

	w{8} = ['z_min: ' txt1 ' at x = ' txt_x1 ' y = ' txt_y1]; 
	w{9} = ['z_max: ' txt2 ' at x = ' txt_x2 ' y = ' txt_y2];

	w{10} = ['scale factor: ' num2str(info1.Scale(1)) ' add_offset: ' num2str(info1.Scale(2))];
	if (~isequal(info2,0))
		txt1 = num2str(info2(6),'%.3f');    txt1 = wipe_zeros(txt1);    % mean
		txt2 = num2str(info2(7),'%.3f');    txt2 = wipe_zeros(txt2);    % stdev
		w{11} = ['mean: ' txt1 '  stdev: ' txt2];
	else
		w{11} = 'WARNING: GRID WAS NOT IN MEMORY SO SOME INFO MIGHT NO BE ENTIRELY CORRECT.';
	end
	if (info2(5))       % We have NaNs, report them also
		w{12} = ['nodes set to NaN: ' num2str(info2(5))];
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
		meta = getappdata(handles.hImg,'meta');
		if (~isempty(meta))						% If we have metadata use the window message to display everything
			InfoMsg = [InfoMsg'; {' '; ' '}; meta];
			message_win('create',InfoMsg,'position','east');
		else
	        msgbox(InfoMsg,'Image Info');
		end
    else
        msgbox('Info missing or nothing to info about?','???')
    end
end

% -----------------------------------------------------------
function txt = wipe_zeros(txt)
	% Wipe zeros at the end of the TXT string
	while (txt(end) == '0')
        txt(end) = [];
	end

% --------------------------------------------------------------------
function att2Hdr(handles,att)
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
	if (~isempty(att.Metadata)),	setappdata(handles.hImg,'meta', att.Metadata),	end
	aux_funs('appP', handles, att.ProjectionRef)		% If we have a WKT proj store it, otherwise clean eventual predecessors

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

    % Maybe not the most apropriate place to do this but ...
    if (isappdata(handles.figure1,'ProjWKT')),    rmappdata(handles.figure1,'ProjWKT'); end
    if (isappdata(handles.figure1,'ProjGMT')),    rmappdata(handles.figure1,'ProjGMT'); end
    if (isappdata(handles.figure1,'Proj4')),      rmappdata(handles.figure1,'Proj4'); end

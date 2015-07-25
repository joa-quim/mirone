function grid_info(handles,X,Y,hdr)
% Report info of grid/image currently displayed in Mirone

%	Copyright (c) 2004-2015 by J. Luis
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

% $Id$

	if (nargin == 3 && strcmp(Y,'gdal'))			% Just extract the relevant info from the attribute struct
		att2Hdr(handles,X);			return
	elseif (nargin == 3 && strcmp(Y,'iminfo'))
		img2Hdr(handles,X);			return
	elseif (nargin == 4 && strcmp(Y,'iminfo'))      % Used when an non-referenced image was sent in input
		img2Hdr(handles,X,hdr);		return
	elseif (nargin == 4 && strcmp(Y,'referenced'))	% Used when referenced img/grid was sent in input
		ref2Hdr(handles,X,hdr);		return
	end

	if (handles.no_file),	return,		end

	if (handles.image_type == 1)					% Image derived from a GMT grdfile
		Z = getappdata(handles.figure1,'dem_z');	% We want the Z for statistics (might not be in argin)
		if (isempty(Z))
			warndlg('Grid too big to be in memory so I cannot report info','Warning')
			return
		end
		nx = size(Z,2);		ny = size(Z,1);
		try
			if (~handles.computed_grid)
				info1 = c_grdinfo(handles.grdname,'hdr_struct');    % info1 is a struct with the GMT grdinfo style
				info_from_GMT_limits = [info1.X_info(1:2) info1.Y_info(1:2)];
			end
			erro = false;
		catch
			erro = true;
		end

		if (handles.computed_grid || erro)
			info1.Title = 'Unknown';	info1.Command = ' ';	info1.Remark = 'INTERNALY COMPUTED GRID';
			info1.X_info(4) = nx;		info1.Y_info(4) = ny;
			info1.Scale = [1 0 -1];
			info1.Registration = 'Registration:  grid (guessed)';
			info_from_GMT_limits = [];
		end

		info2 = grdutils(Z,'-H');				% info2 is a vector with [z_min z_max i_zmin i_zmax n_nans mean std]
		info2(3:4) = info2(3:4) + 1;			% Info from grdutils is zero based

		% In case grid is Pix reg the handles.head no longer holds that info, so check against info1
		head = handles.head;
		if (~isempty(info_from_GMT_limits) && (abs(head(1) - info_from_GMT_limits(1)) > head(8)/4) && strncmp(info1.Registration,'Pixel',5) )
			head(1:4) = info_from_GMT_limits;
			head(7) = 1;
		end

		w{1} = ['Title: ' info1.Title];
		w{2} = ['Command: ' info1.Command];
		w{3} = ['Remark: ' info1.Remark];
		w{4} = info1.Registration;
		w{5} = ['grdfile format: #' sprintf('%d',info1.Scale(3))];
		txt1 = sprintf('%.8g',head(1));			% x_min
		txt2 = sprintf('%.8g',head(2));			% x_max
		txt3 = sprintf('%.8g',head(8));			% x_inc
		w{6} = ['x_min: ' txt1 '  x_max: ' txt2 '  x_inc: ' txt3 '  nx: ' sprintf('%d',info1.X_info(4))];
		txt1 = sprintf('%.8g',head(3));			% y_min
		txt2 = sprintf('%.8g',head(4));			% y_max
		txt3 = sprintf('%.8g',head(9));			% y_inc
		w{7} = ['y_min: ' txt1 '  y_max: ' txt2 '  y_inc: ' txt3 '  ny: ' sprintf('%d',info1.Y_info(4))];
		txt1 = sprintf('%.8g',info2(1));		% z_min
		txt2 = sprintf('%.8g',info2(2));		% z_max

		if (head(7)),	half = 0.5;
		else			half = 0;
		end
		xx_min = head(1) + (fix((info2(3)-1) / ny) + half) * head(8);	% x of z_min
		xx_max = head(1) + (fix((info2(4)-1) / ny) + half) * head(8);	% x of z_max
		yy_min = head(3) + (rem(info2(3)-1, ny) + half) * head(9);		% y of z_min
		yy_max = head(3) + (rem(info2(4)-1, ny) + half) * head(9);		% y of z_max
		txt_x1 = sprintf('%.8g',xx_min);
		txt_x2 = sprintf('%.8g',xx_max);
		txt_y1 = sprintf('%.8g',yy_min);
		txt_y2 = sprintf('%.8g',yy_max);

		w{8} = ['z_min: ' txt1 '   at x = ' txt_x1 '   y = ' txt_y1]; 
		w{9} = ['z_max: ' txt2 '   at x = ' txt_x2 '   y = ' txt_y2];

		w{10} = ['scale factor: ' num2str(info1.Scale(1)) '        add_offset: ' num2str(info1.Scale(2))];
		if (~isequal(info2,0))
			w{11} = sprintf('mean: %.8g    stdev: %.8g',info2(6), info2(7));
		else
			w{11} = 'WARNING: GRID WAS NOT IN MEMORY SO SOME INFO MIGHT NO BE ENTIRELY CORRECT.';
		end
		if (info2(5))       % We have NaNs, report them also
			w{12} = ['nodes set to NaN: ' sprintf('%d',info2(5))];
		end
		hMsg = message_win('create',w, 'figname','Grid Info', 'button','yes');

		% Give message_win the info needed to plot the grid's min/Max
		ud.hAxe = handles.axes1;	ud.p1 = [xx_min yy_min];	ud.p2 = [xx_max yy_max];
		set(hMsg, 'UserData', ud)
	else
		InfoMsg = getappdata(handles.axes1,'InfoMsg');
		if (~isempty(InfoMsg))
			meta = getappdata(handles.hImg,'meta');
			if (~isempty(meta))
				InfoMsg = [InfoMsg; {' '; ' '}; meta];
			end
			hMsg = message_win('create',InfoMsg, 'figname', 'Grid Info');
		else
			hMsg = msgbox('Info missing or nothing to info about?','???');
		end
	end

	% Add the Message figure handle to the carraças list
	plugedWin = getappdata(handles.figure1,'dependentFigs');
	plugedWin = [plugedWin hMsg];
	setappdata(handles.figure1,'dependentFigs',plugedWin);

	if (handles.version7 < 8.4 && strncmp(computer,'PC',2))
		WindowAPI(hMsg, 'TopMost')
	end

% --------------------------------------------------------------------
function att2Hdr(handles,att)
	% Fill a header with the info from the att struct issued by gdalread

	w = cell(13,1);
    w{1} = ['Driver : ' att.DriverShortName];
    w{2} = att.ProjectionRef;
    w{3} = [];
    w{4} = 'Proj4 string:';
	try
		if (~isempty(att.ProjectionRef))
			w{5} = ogrproj(att.ProjectionRef);		% Get equivalent in Proj4 format
		else
			w{5} = 'Ghrrr!! Don''t know (tried to find it but failed)';
		end
	catch
		w{5} = 'Fiu Fiu Fiu! Translating from rubish ... obviously failed';
	end
    w{6} = [];
    w{7} = ['Width:  ' num2str(att.RasterXSize) '    Height:  ' num2str(att.RasterYSize)];
    w{8} = ['Pizel Size:  (' num2str(att.GMT_hdr(8)) ',' num2str(att.GMT_hdr(9)) ')'];
    w{9} = 'Projected corner coordinates';
    w{10}  = ['   Xmin:  ' num2str(att.Corners.LL(1)) '    Xmax: ' num2str(att.Corners.UR(1))];
    w{11}  = ['   Ymin:  ' num2str(att.Corners.LL(2)) '    Ymax: ' num2str(att.Corners.UR(2))];
    if (~isempty(att.GEOGCorners))
		w{12} = 'Geographical corner coordinates';
		w{13} = ['   Lon min:  ' att.GEOGCorners{1,1} '    Lon max: '  att.GEOGCorners{4,1}];
		w{14} = ['   Lat min:  ' att.GEOGCorners{1,2} '    Lat max: '  att.GEOGCorners{2,2}];
    end
    w{end+1} = ['   Zmin:  ' num2str(att.GMT_hdr(5)) '   Zmax: ' num2str(att.GMT_hdr(6))];
    w{end+1} = ['Color Type:  ' att.ColorInterp];
    
    setappdata(handles.axes1,'InfoMsg',w)
	if (~isempty(att.Metadata)),	setappdata(handles.hImg,'meta', att.Metadata),	end
	aux_funs('appP', handles, att.ProjectionRef)		% If we have a WKT proj store it, otherwise clean eventual predecessors

% --------------------------------------------------------------------
function ref2Hdr(handles, srs, img)
% Deal with the cases when a grid or an image with a WKT proj string was sent in input to Mirone
% We than build a fake minimalist 'att' like it outputs from gdalread and call ATT2HDR to do the job

	att.DriverShortName = 'Input array';
	att.ProjectionRef = srs;
	att.RasterXSize = size(img,2);
	att.RasterYSize = size(img,1);
	att.Corners.LL = [handles.head(1) handles.head(3)];
	att.Corners.UR = [handles.head(2) handles.head(4)];
	att.GEOGCorners = [];
	att.GMT_hdr = handles.head;
	att.Metadata = '';
	att.ColorInterp = 'Don''t know';
	att2Hdr(handles,att)				% That's where the info string is built and stored

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

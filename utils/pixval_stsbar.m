function pixval_stsbar(arg1)
% Based on the defunct pixval, but heavily hacked in several ways.

% Coffeeright J. Luis 2004-2017

% $Id$

	if (nargin == 0),	arg1 = [];		end

	if ischar(arg1) % execute a callback
		switch lower(arg1)
			case 'exit'
				displayBar = findobj(gcf, 'Tag', 'pixValStsBar');
				userData = get(displayBar, 'UserData');
				if (~isempty(userData)),	uirestore_j(userData.oldFigureUIState),		end
				if ishandle(displayBar),    delete(displayBar);     end
				return
		end 
	else
		if (nargin == 1 && ishandle(arg1))
			hFig = arg1;
		else
			hFig = gcf;
		end
		userData.hFig = hFig;
		userData.displayMode = 'normal';
		% Save the interactive state of the figure. 
		userData.oldFigureUIState = uisuspend_j(hFig, true);
		% Create position vectors for Text
		handsBar = getappdata(hFig,'CoordsStBar');
		tmp1 = get(handsBar(3),'YData');
		tmp2 = get(handsBar(2),'XData');
		p = get(handsBar(1),'Pos');			% Assumes units are pixels.
		pos = [p(1) tmp1(1) tmp2(2) abs(tmp1(2)-tmp1(1))] + [4 2 -5 -4];

		% Create the display bar
		userData.buttonDownImage = 0;		% Image 0 never exists      
		DisplayBar = uicontrol('Parent', hFig, ...
							'Style','text', ...
							'Units','pixels', ...
							'Position',pos, ...
							'Foreground', [0 0 0], ...
							'Horiz','left', ...
							'Tag', 'pixValStsBar', ...
							'fontname', 'FixedWidth', ...
							'BusyAction', 'queue', ...
							'enable', 'inactive', ...
							'Interruptible', 'off');

		% Register the motion function and button up function
		set(hFig, 'WindowButtonMotionFcn', {@PixvalMotionFcn,DisplayBar})
		% Reset the original 'KeyPressFcn' since we dont need a private one here
		set(hFig, 'KeyPressFcn', userData.oldFigureUIState.KeyPressFcn)
		handles = guidata(hFig);
		userData.is_projected = handles.is_projected;
		userData.head  = handles.head;
		userData.axes1 = handles.axes1;
		userData.geog  = handles.geog;
		if isempty(getappdata(hFig,'dem_x')),   userData.haveGrid = 0;
		else                                    userData.haveGrid = 1;
		end
		userData.toProjPT = 0;

		if (isempty(findobj(hFig, '-depth',2, 'type','image')))	% To know if being called from Mirone or Ecran
			userData.IAmImage = false;
		else
			userData.IAmImage = true;
		end

		set(DisplayBar, 'UserData', userData);
		PixvalMotionFcn([], [], DisplayBar);
	end

%-----------------------------------------------------------------------------------------
function PixvalMotionFcn(obj, evt, displayBar)

	userData = get(displayBar, 'UserData');

	if (strcmp(userData.displayMode, 'normal'))			% See if we're above the displayBar
		  [hImg, imageType, img, x, y] = OverImage(userData.hFig);
		  if (hImg ~= 0)								% Update the Pixel Value display
			  UpdatePixelValues(userData.hFig, hImg, imageType, displayBar, img, x, y);
		  end
	elseif strcmp(userData.displayMode, 'distance')		% If we're in distance mode and in another image, clean up a bit.
	   [hImg, imageType, img, x, y] = OverImage(userData.hFig);
	   if (hImg ~= 0)
		  if (hImg == userData.buttonDownImage)
			  UpdatePixelValues(userData.hFig, hImg, imageType, displayBar, img, x, y);
		  end
	   end
	end
	set(displayBar, 'UserData', userData);

%-----------------------------------------------------------------------------------------
function [hImg,imageType,img,x,y] = OverImage(hFig)
% Return the index of which image we are over, and return a 0 if we aren't above an image.
	images = findobj(hFig, '-depth',2, 'type','image');
	if isempty(images)
		[hImg, imageType, x, y] = OverAxes(hFig);		% Than see if we are over an axes (the 'ecran' case)
		img = [];
		return
	end
	% Make sure that the Image's Button Down & Up functions will queue
	set(images, 'ButtonDownFcn', {@ButtonDownOnImage,hFig}, 'Interruptible','off', 'BusyAction','Queue');
	hAxes = get(images, {'Parent'});
	axCurPt = get([hAxes{:}], {'CurrentPoint'});

	% Loop over the axes, see if we are above any of them
	hImg = 0;
	for (k = 1:numel(hAxes))
		XLim = get(hAxes{k}, 'XLim');	YLim = get(hAxes{k}, 'YLim');
		pt = axCurPt{k};
		x = pt(1,1);	y = pt(1,2);
		if (x >= XLim(1) && x <= XLim(2) && y >= YLim(1) && y <= YLim(2))
			hImg = images(k);
			break
		end
	end

	% Figure out image type
	if (hImg ~= 0)
		[img, flag] = get_image_info(hImg);
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
				error('Images:pixval_stsbar:invalidImage', '%s', ['Invalid image, GETIMAGE returned flag = ' flag '.']);
		end
	else
		hImg = 0;	imageType = '';		img = [];		x = 0;	y = 0;
	end

%-----------------------------------------------------------------------------------------
function [hAx, imageType, x, y] = OverAxes(hFig)
% Return the index of which axes we are over, and return a 0 if we aren't above an axes.
	hAxes = findobj(hFig, '-depth',1, 'type','axes');
	if isempty(hAxes)
		hAx = 0;	imageType = '';		x = 0;	y = 0;
		return
	end
	% Make sure that the Axes's Button Down & Up functions will queue
	set(hAxes, 'ButtonDownFcn', {@ButtonDownOnImage,hFig}, 'Interruptible','off', 'BusyAction','Queue');

	% Loop over the axes, see if we are above any of them
	axCurPt = get(hAxes, {'CurrentPoint'});

	% Loop over the axes, see if we are above any of them
	hAx = 0;
	for (k = 1:numel(hAxes))
		XLim = get(hAxes(k), 'XLim');	YLim = get(hAxes(k), 'YLim');
		pt = axCurPt{k};
		x = pt(1,1);	y = pt(1,2);
		if (x >= XLim(1) && x <= XLim(2) && y >= YLim(1) && y <= YLim(2))
			hAx = hAxes(k);
			break
		end
	end

	if (hAx ~= 0)
		imageType = 'intensity';		% To true either but we use this code to find our way in the algo.
	else
		hAx = 0;	imageType = '';		x = 0;	y = 0;
	end

%-----------------------------------------------------------------------------------------
function UpdatePixelValues(hFig, hImg, imageType, displayBar, img, x, y)
%   This is the motion callback for when we are displaying pixels.
%   Either we are in automatic display mode and the mouse pointer is
%   moving around or we are in normal mode and there has been a button-down
%   but not yet a button up. I get the current point and update the string.

	userData = get(displayBar, 'UserData');
	axHandle = userData.axes1;
	rcMode = false;

	if (userData.haveGrid)
		pMode = getappdata(hFig,'PixelMode');		rcMode = getappdata(hFig,'RCMode');
		if (isempty(pMode)),		pMode = false;		end
		if (isempty(rcMode)),		rcMode = false;		end
 		if (pMode || rcMode)
			% Inside each grid cell, which is a pixel on the screen, display only the grid node value
			rows = size(img, 1);		cols = size(img, 2);
			rp = getPixel_coords(rows, get(hImg,'YData'),y);
			cp = getPixel_coords(cols, get(hImg,'XData'),x);
			r = min(rows, max(1, round(rp)));   c = min(cols, max(1, round(cp)));
			if (pMode)
				Z = getappdata(hFig,'dem_z');
				pixel = double(Z(r,c));
			else
				pixel = [r c];
			end
		else
			pixel = bi_linear(getappdata(hFig,'dem_x'),getappdata(hFig,'dem_y'),getappdata(hFig,'dem_z'),x,y);
		end

	else						% work on a image type
		rows = size(img, 1);		cols = size(img, 2);
		if (userData.IAmImage)
			rp = getPixel_coords(rows, get(hImg,'YData'),y);
			cp = getPixel_coords(cols, get(hImg,'XData'),x);
		else					% Assume hImg is in fact an Axes handle
			rp = 1;		cp = 1;		rows = 1;	cols = 1;	% None of these are really used but need to exist 
			pixel = 0;
		end
		r = min(rows, max(1, round(rp)));   c = min(cols, max(1, round(cp)));
		if strcmp(imageType,'indexed')
			map = get(userData.hFig, 'Colormap');
			idx = img(r,c);		ind = double(idx);
			if (~isa(idx,'double')), 	idx = double(idx)+1;   end
			idx = round(idx);
			if (idx <= size(map,1)),	pixel = [ind round(map(idx,:)*255)];
			else						pixel = [ind round(map(end,:)*255)];
			end
		elseif (userData.IAmImage)		% Otherwise 'pixel' has no real meaning
			pixel = double(img(r,c,:));
		end
	end

	if (userData.toProjPT > 0 && getappdata(hFig,'DispInGeogs'))
        if (userData.toProjPT == 1)
			xy_p = ogrproj([x y],userData.projStruc);
        else
			xy_p = c_mapproject([x y], userData.opt_R, '-F', '-I', userData.projGMT{:});
        end
		x = xy_p(1);    y = xy_p(2);
	end

	% Find the coordinate output format
	labelType = getappdata(axHandle,'LabelFormatType');
	if (isempty(labelType))     % Actualy this is an error test because labelType should never be empty.
		labelType = 'NotGeog';
	end
	% Make a copy of x & y to use in the distance case. Needed when input was transformed to a string below
	x1 = x;		y1 = y;
	switch labelType
		case {'DegDec', 'NotGeog'}       % This's the default. Just build the format string
			form_xy = ' %8.3f,%7.3f =';
		case 'DegMin'
			out_x = degree2dms(x,'DDMM',0,'str');		x = [out_x.dd ':' out_x.mm];
			out_y = degree2dms(y,'DDMM',0,'str');		y = [out_y.dd ':' out_y.mm];
			form_xy = ' %s, %s =';
		case 'DegMinDec'
			out_x = degree2dms(x,'DDMM.x',2,'str');		x = [out_x.dd ':' out_x.mm];
			out_y = degree2dms(y,'DDMM.x',2,'str');		y = [out_y.dd ':' out_y.mm];
			form_xy = ' %s, %s =';
		case 'DegMinSec'
			out_x = degree2dms(x,'DDMMSS',0,'str');		x = [out_x.dd ':' out_x.mm ':' out_x.ss];
			out_y = degree2dms(y,'DDMMSS',0,'str');		y = [out_y.dd ':' out_y.mm ':' out_y.ss];
			form_xy = ' %s, %s =';
		case 'DegMinSecDec'
			out_x = degree2dms(x,'DDMMSS.x',2,'str');	x = [out_x.dd ':' out_x.mm ':' out_x.ss];
			out_y = degree2dms(y,'DDMMSS.x',2,'str');	y = [out_y.dd ':' out_y.mm ':' out_y.ss];
			form_xy = ' %s, %s =';
		case 'Date'
			x = datestr(x,'dd-mmm-yyyy HH:MM:SS');
			form_xy = ' %s,%7.3f =';
	end
 
	% figure out the new string
	switch userData.displayMode
		case 'normal'			% Just display Z (or intensity) information
			if strcmp(imageType,'rgb') || strcmp(imageType,'indexed')
				if isa(img, 'uint8') && strcmp(imageType,'rgb')
					if (userData.haveGrid)			% Hacked here
						pixval_str = sprintf([form_xy ' %6.3f'], x,y,pixel(1:end));
					else
						pixval_str = sprintf([form_xy ' %3d,%3d,%3d'], x,y,pixel(1:end));
					end
				elseif isa(img, 'uint16') && strcmp(imageType,'rgb')
					pixval_str = sprintf([form_xy ' %5d,%5d,%5d'], x,y,pixel(1:end));  
				elseif islogical(img) && strcmp(imageType,'rgb')
					pixval_str = sprintf([form_xy ' %1d,%1d,%1d'], x,y,pixel(1:end));  
				else	% all indexed images use double precision colormaps
					if (userData.haveGrid)				% Hacked here
						pixval_str = sprintf([form_xy ' %6.3f'], x,y,pixel(1:end));
						if (rcMode)
							pixval_str = sprintf([form_xy(1:end-1) '\tRow = %d  Col = %d'], x,y,pixel(1:end));
						end
					else
						if (numel(pixel) == 4)		% Display indexed color plus matrix value
							pixval_str = sprintf([form_xy ' %3d,%3d,%3d [%d]'], x,y,pixel(2:4),pixel(1));
						else
							pixval_str = sprintf([form_xy ' %.0f,%.0f,%.0f'], x,y,pixel(1:end));
						end
					end
				end
			else		% intensity
				if (~userData.IAmImage),		pixval_str = sprintf(form_xy(1:end-1),x,y);		% A call from Ecran
				elseif isa(img, 'uint8'),		pixval_str = sprintf([form_xy ' %g'],x,y,pixel(1));
				elseif isa(img, 'uint16'),		pixval_str = sprintf([form_xy ' %g'],x,y,pixel(1));
				elseif islogical(img),			pixval_str = sprintf([form_xy ' %g'],x,y,pixel(1));
				else							pixval_str = sprintf([form_xy ' %6.4f'],x,y,pixel(1));
				end
			end

		case 'distance'
			handles = guidata(hFig);
			delta_x = (x1 - userData.x0);		delta_y = (y1 - userData.y0);
			set(userData.line, 'XData', [userData.x0 x1], 'YData', [userData.y0 y1]);
			
			% The difference is that handles.geog can be changed in mid-time but the one in userData is frozen to its initial value.
			geog = handles.geog;
			if (~userData.IAmImage),		geog = userData.geog;	end		% Serving 'ecran()'
			if (geog)
				switch handles.DefineMeasureUnit(1)     % I have to do it here to allow midtime changes in preferences
					case 'n',		scale = 1852;   str_dist = 'dist(NM)';		% Nautical miles
					case 'k',		scale = 1000;   str_dist = 'dist(km)';		% Kilometers
					case 'm',		scale = 1;      str_dist = 'dist(m)';		% Meters or user unites
					case 'u',		scale = 1;      str_dist = 'dist(usr)';		% Meters or user unites
				end
				if (delta_x == 0 && delta_y == 0)		% Its true on every first click 
					dist = 0;	az = 0;
				else
					dist = vdist(userData.y0,userData.x0,y1,x1,handles.DefineEllipsoide([1 3])) / scale;
					D2R = pi/180;
					lat1 = userData.y0*D2R;     lon1 = userData.x0*D2R;
					lat2 = y1*D2R;          lon2 = x1*D2R;
					f2 = cos(lat1) * sin(lat2);
					f3 = sin(lat1) * cos(lat2) * cos(lon2-lon1);
					az = 90 - atan2(cos(lat2) * sin(lon2-lon1), f2-f3) / D2R;
				end
			else
				dist = sqrt(delta_x^2 + delta_y^2);     str_dist = 'dist(usr)';
				az = atan2(delta_y, delta_x) * 180/pi;
				if(strcmp(get(handles.axes1,'YDir'),'reverse')),    az = -az;   end
			end

			if strcmp(imageType,'rgb') || strcmp(imageType,'indexed')
				if (isa(img, 'uint8') &&  strcmp(imageType,'rgb') && ~userData.haveGrid)
					pixval_str = sprintf([form_xy ' %3d,%3d,%3d  '  str_dist ' = %3.3f ang = %.1f'], x,y,pixel(1:end),dist,az);
				elseif (isa(img, 'uint8') &&  strcmp(imageType,'rgb') && userData.haveGrid)
					pixval_str = sprintf([form_xy ' %g '  str_dist ' = %4.4f ang = %.1f'], x,y,pixel(1:end),dist,az);
				elseif isa(img, 'uint16') &&  strcmp(imageType,'rgb')
					pixval_str = sprintf([form_xy ' %5d,%5d,%5d  '  str_dist ' = %3.3f ang = %.1f'], x,y,pixel(1:end),dist,az);
				elseif islogical(img) &&  strcmp(imageType,'rgb')
					pixval_str = sprintf([form_xy ' %1d,%1d,%1d  '  str_dist ' = %3.3f ang = %.1f'], x,y,pixel(1:end),dist,az);
				else		% all indexed images use double precision colormaps
					pixval_str = sprintf([form_xy ' %g  '  str_dist ' = %4.4f ang = %.1f'], x,y,pixel(1),dist,az);
				end
			else		% intensity
				if (~userData.IAmImage)			% A call from Ecran
					pixval_str = sprintf([form_xy '  dist_X = %g dist_Y = %g'], x, y, delta_x, delta_y);
				elseif isa(img, 'uint8')
					pixval_str = sprintf([form_xy ' %3d  dist = %3.3f ang = %.1f'], x, y, pixel(1), dist,az);
				elseif isa(img, 'uint16')
					pixval_str = sprintf([form_xy ' %5d  dist = %3.3f ang = %.1f'], x, y, pixel(1), dist,az);
				elseif islogical(img)
					pixval_str = sprintf([form_xy ' %1d  dist = %3.3f ang = %.1f'], x, y, pixel(1), dist,az);  
				else
					pixval_str = sprintf([form_xy ' %6.4f  dist = %3.3f ang = %.1f'], x, y, pixel(1), dist,az);
				end
			end
	end
	set(displayBar, 'String', pixval_str, 'UserData', userData);

%-----------------------------------------------------------------------------------------
function ButtonDownOnImage(hImg, evt, hFig)
	stype = get(hFig,'selectiontype');
	IAmImage = strcmp(get(hImg, 'type'),'image');		% To know if called from Mirone or Ecran
	if (strcmp(stype,'alt'))		% A right-click, either return or pass control to ...
		if (IAmImage && ~isempty(getappdata(hFig,'LinkedTo')))
			linkDisplays(hFig)
		else						% uictx to options to fill gaps (NaN holes) 
			Z = getappdata(hFig,'dem_z');
			if (~isempty(Z))
				pt = get(get(hImg, 'Parent'), 'CurrentPoint');
				[rows, cols, rp] = size(get(hImg,'CData'));
				rp = getPixel_coords(rows, get(hImg,'YData'),pt(1,2));
				cp = getPixel_coords(cols, get(hImg,'XData'),pt(1,1));
				r = min(rows, max(1, round(rp)));	c = min(cols, max(1, round(cp)));
				z = Z(r,c);
				if (isnan(z))		% Create the UIContextMenu
					cmenu = uicontextmenu('Parent', hFig, 'Tag','clickedHole');
					set(hImg, 'UIContextMenu', cmenu);
					uimenu(cmenu, 'Label', 'Digitize body', 'Call', 'mirone(''ImageEdgeDetect_CB'',guidata(gcbo),''apalpa_body'')');
					item0 = uimenu(cmenu, 'Label', 'Digitize body with pad');
					uimenu(item0, 'Label', 'pad = 1 cell', 'Call', 'mirone(''ImageEdgeDetect_CB'',guidata(gcbo),''apalpa_body_1'')');
					uimenu(item0, 'Label', 'pad = 2 cell', 'Call', 'mirone(''ImageEdgeDetect_CB'',guidata(gcbo),''apalpa_body_2'')');
					uimenu(item0, 'Label', 'other...',  'Call', @other_pad);
					uimenu(cmenu, 'Label', 'Digitize this hole', 'Call', 'mirone(''ImageEdgeDetect_CB'',guidata(gcbo),''apalpa_um'')', 'Sep','on');
					uimenu(cmenu, 'Label', 'Inpaint this hole', 'Call', 'inpaint_nans(guidata(gcbo), ''single'')');
					item1 = uimenu(cmenu, 'Label', 'Fill this hole ...', 'Sep','on');
					uimenu(item1, 'Label', 'with 2nd Grid (sharp edges)', 'Call', 'transplants([], ''one_sharp'', true, guidata(gcbo))');
					uimenu(item1, 'Label', 'with 2nd Grid (smooth edges)','Call', 'transplants([], ''one_smooth'', true, guidata(gcbo))');
					item2 = uimenu(cmenu, 'Label', 'Fill all holes ...');
					uimenu(item2, 'Label', 'with 2nd Grid (sharp edges)', 'Call', '');
					uimenu(item2, 'Label', 'with 2nd Grid (smooth edges)', 'Call', '');
				else
					delete(findobj('Tag','clickedHole'))	% Delete uctx so that it does show up on right-clicks on ANY point
				end
			end
		end
		return
	elseif (strcmp(stype,'extend'))
		if (IAmImage),	magnify(get(hImg, 'Parent')),	end
		return
	end
	displayBar = findobj(hFig, 'Tag', 'pixValStsBar');
	userData = get(displayBar, 'UserData');
	if (IAmImage)
		hAxes = get(hImg, 'Parent');
	else		% When called from 'ecran' we have no image, just axes
		hAxes = hImg;
	end
	% Set the initial point (x0,y0)
	pt = get(hAxes, 'CurrentPoint');
	userData.x0 = pt(1,1);
	userData.y0 = pt(1,2);
	userData.line = line('Parent', hAxes, 'color', [1 0 0], ...
		'Xdata', [userData.x0 userData.x0],'Ydata', [userData.y0 userData.y0],'Tag','polyline');
	userData.displayMode = 'distance';
	userData.buttonDownImage = hImg;
	set(displayBar, 'UserData', userData);
	set(userData.hFig, 'WindowButtonUpFcn', {@BackToNormalPixvalDisplay, displayBar});
	PixvalMotionFcn([], [], displayBar);

%-----------------------------------------------------------------------------------------
function other_pad(obj,evt)
% Set the padding for a digitize request
	resp  = str2double(inputdlg({'Enter new pad width in number of cells'}, 'Pad width', [1 30], {'3'}));
	if (isnan(resp) || resp < 0),	return,		end
	mirone('ImageEdgeDetect_CB', guidata(obj), sprintf('apalpa_body_%d', round(resp)))

%-----------------------------------------------------------------------------------------
function BackToNormalPixvalDisplay(obj, evt, displayBar)
	userData = get(displayBar, 'UserData');
	stype = get(userData.hFig,'selectiontype');
	if (strcmp(stype,'extend'))			% A mouse midle-button click, otherwise let the line survive
		if (userData.IAmImage)
			draw_funs(userData.line,'line_uicontext')		% Set lines's uicontextmenu
		else
			draw_funs([], 'set_line_uicontext_XY', userData.line)	% Set lines's uicontextmenu
		end
		handles = guidata(userData.hFig);
		% Use 'try' because the fields below may not exist in ecran()
		try		set(userData.line, 'Color', handles.DefLineColor),	end
		try		set(userData.line, 'LineWidth', handles.DefLineThick),	end

		% But we also need to remove several uicontexts when the line is on a ecran() fig
	else
		delete(userData.line)
	end
	userData.line = [];     userData.x0 = []; userData.y0 = [];
	set(userData.hFig, 'WindowButtonUpFcn', '');
	userData.displayMode = 'normal';
	userData.buttonDownImage = 0;
	set(displayBar, 'UserData', userData);
	PixvalMotionFcn([], [], displayBar);

% -------------------------------------------------------------------------------------
function pix_coords = getPixel_coords(img_length, XData, axes_coord)
% Convert coordinates from axes (real coords) to image (pixel) coordinates.
% IMG_LENGTH is the image width (n_columns)
% XDATA is the image's [x_min x_max] in axes coordinates
% AXES_COORD is the (x,y) coordinate of the point(s) to be converted

	slope = (img_length - 1) / (XData(end) - XData(1));
	if ((XData(1) == 1) && (slope == 1))
		pix_coords = axes_coord;
	else
		pix_coords = slope * (axes_coord - XData(1)) + 1;
	end

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
		if ((ndims(A) == 3) && (size(A,3) == 3))	% We have an RGB image
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

% ---------------------------------------------------------------------------------
function linkDisplays(hFig)
%
	linkedFig = getappdata(hFig,'LinkedTo');
	if (~ishandle(linkedFig)),		return,		end			% Linked Fig was killed
	handMir_linked  = guidata(linkedFig);					% Fish the pair of handles
	handMir_clicked = guidata(hFig);

	% This chunk is a trick to patch a f... ML bug that when it pleases stop calling the unset part
	linked_old = getappdata(0,'linkedStruct');
	if (~isempty(linked_old) && linked_old.toggled)
		set(handMir_clicked.hImg,'CData', linked_old.imgClicked, ...
			'XData',linked_old.xdataClicked, 'YData',linked_old.ydataClicked)
		set(hFig, 'colormap', linked_old.cmapClicked)
		linked_old.toggled = false;
		setappdata(0, 'linkedStruct', linked_old)
		return
	end

	imgLinked  = get(handMir_linked.hImg,'CData');
	linked.imgClicked = get(handMir_clicked.hImg,'CData');			% We need to backup this before it is replaced by the linked
	linked.cmapClicked = get(hFig,'colormap');
	linked.xdataClicked = get(handMir_clicked.hImg, 'XData');	linked.ydataClicked = get(handMir_clicked.hImg, 'YData');

	set(handMir_clicked.hImg,'CData', get(handMir_linked.hImg,'CData'), ...
		'XData',get(handMir_linked.hImg, 'XData'), 'YData',get(handMir_linked.hImg, 'YData'))
	if (ndims(imgLinked) == 2)
		set(hFig, 'colormap',get(handMir_linked.figure1,'colormap'))
	end

	linked.toggled = true;
	setappdata(0, 'linkedStruct', linked)

	set(hFig, 'WindowButtonUpFcn', {@link_bu, hFig, handMir_clicked.hImg});

% ---------------------------------------------------------------------------------
function link_bu(obj, event, hFig, hImg)
% Undo the image swapping and reset the 'WindowButtonUpFcn' to ''

	handMir_clicked = guidata(hFig);
	linked = getappdata(0,'linkedStruct');
	set(handMir_clicked.hImg,'CData', linked.imgClicked)
	set(hFig,'colormap', linked.cmapClicked, 'WindowButtonUpFcn', '')
	set(hImg, 'XData',linked.xdataClicked, 'YData',linked.ydataClicked)
	linked.toggled = false;
	setappdata(0, 'linkedStruct', linked)

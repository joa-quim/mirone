function varargout = image_enhance(varargin)
%   IMCONTRAST(H) creates an Adjust Contrast tool associated with the image
%   specified by the handle H. H may be an image, axes, or figure handle.
%
%   The Adjust Contrast tool presents a scaled histogram of pixel values
%   (overly represented pixel values are truncated for clarity). Dragging on the
%   left red bar in the histogram display changes the minimum value. The minimum
%   value (and any value less than the minimum) displays as black. Dragging on the
%   right red bar in the histogram changes the maximum value. The maximum value
%   (and any value greater than the maximum) displays as white. Values in between
%   the red bars display as intermediate shades of gray.
%
%   Together the minimum and maximum values create a "window". Stretching the window
%   reduces contrast. Shrinking the window increases contrast. Changing the center of
%   the window changes the brightness of the image. It is possible to manually enter
%   the minimum, maximum, width, and center values for the window. Changing one value
%   automatically updates the other values and the image.
%

%	Copyright (c) 2004-2012 by J. Luis
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

	if (isempty(varargin))	return,		end

	hObject = figure('Vis','on');
	image_enhance_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hObject,'right')

	handles.new_pointer = [NaN	NaN	NaN	NaN	NaN	NaN	NaN	2	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN; ...
			NaN	NaN	NaN	NaN	NaN	NaN	2	1	2	NaN	NaN	NaN	NaN	NaN	NaN	NaN; ...
			NaN	NaN	NaN	NaN	NaN	2	1	1	1	2	NaN	NaN	NaN	NaN	NaN	NaN; ...
			NaN	NaN	NaN	NaN	2	1	1	1	1	1	2	NaN	NaN	NaN	NaN	NaN; ...
			NaN	NaN	NaN	2	1	1	1	1	1	1	1	2	NaN	NaN	NaN	NaN; ...
			2	2	2	2	2	2	1	1	1	2	2	2	2	2	2	NaN; ...
			2	1	1	1	1	1	1	1	1	1	1	1	1	1	2	NaN; ...
			2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	NaN; ...
			2	1	1	1	1	1	1	1	1	1	1	1	1	1	2	NaN; ...
			2	2	2	2	2	2	1	1	1	2	2	2	2	2	2	NaN; ...
			NaN	NaN	NaN	2	1	1	1	1	1	1	1	2	NaN	NaN	NaN	NaN; ...
			NaN	NaN	NaN	NaN	2	1	1	1	1	1	2	NaN	NaN	NaN	NaN	NaN; ...
			NaN	NaN	NaN	NaN	NaN	2	1	1	1	2	NaN	NaN	NaN	NaN	NaN	NaN; ...
			NaN	NaN	NaN	NaN	NaN	NaN	2	1	2	NaN	NaN	NaN	NaN	NaN	NaN	NaN; ...
			NaN	NaN	NaN	NaN	NaN	NaN	NaN	2	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN; ...
			NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN]';

	handMir = guidata(varargin{1});
	handles.hMirAxes = handMir.axes1;
	handles.hMirImg = handMir.hImg;
	img = (get(handles.hMirImg,'CData'));

	if (isempty(img))
		errordlg('C''mon load a file first. It''s logic, isn''t it?','Error')
		delete(hObject);
		return
	elseif (ndims(img) == 3)
		handles.isRGB = 1;
		set(handles.edit_percentOutliersR,'BackgroundColor','r')
		handles.percentOutliers(1:3) = [2 2 2];
	elseif (ndims(img) == 2)
		set(handles.hMirImg,'CDataMapping','scaled')   % We will change CLim in intensity images
		handles.isRGB = 0;
		handles.percentOutliers = 2;
		% OK, here we must reshape the GUI, but we can delete some uis first
		delete(handles.edit_percentOutliersG);    delete(handles.edit_percentOutliersB)
		delete(handles.radio_R);    delete(handles.radio_G);    delete(handles.radio_B);
		delete(handles.push_scaterPlot);    delete(handles.push_decorrStrectch);
		pos_fig = get(hObject,'pos');
		pos_a1 = get(handles.axes1,'pos');
		new_height = pos_fig(4) - pos_a1(2) + 40;
		all_childs = get(hObject,'Children');
		set(hObject,'pos',[pos_fig(1:3) new_height])
		for i = 1:length(all_childs)
			set(all_childs(i),'pos',get(all_childs(i),'pos')+[0 -(pos_fig(4) - new_height) 0 0])
		end
		set(handles.axes2,'Visible','off')
		pos_a1 = get(handles.axes1,'pos');      % Get new axe1 position to be used as reference
		posf = get(handles.frame_axes,'pos');
		set(handles.frame_axes,'pos', [posf(1) pos_a1(2)-23 posf(3) pos_a1(4)+42])
		post = get(handles.text_tip,'pos');
		set(handles.text_tip,'pos',[post(1) pos_a1(2)-40 post(3) post(4)])
		post = get(handles.text_stat,'pos');
		set(handles.text_stat,'pos',[post(1) pos_a1(2)-40 post(3) post(4)])
		handles.currAxes = zeros(1,3);
	else
		errordlg('This tool only works with Intensity or RGB images.','Error')
		delete(hObject);    return
	end

	%------------ Give a Pro look (3D) to the frame boxes  -------------------------------
	new_frame3D(hObject, [handles.text_dataRange handles.text_window handles.text_scaleRange])
	%------------- END Pro look (3D) -----------------------------------------------------

	handles.imgClass = class(img);
	handles.satistic(1) = {' '};

	if (~handles.isRGB)
		handles.currAxes = 1;	set(handles.figure1,'CurrentAxes',handles.axes1);
		localImhist(handles,img);
		handles = guidata(handles.figure1);			% The handles has been saved in plot_results
		handles.minCData = double(min(min(img)));	handles.maxCData = double(max(max(img)));
	else
		handles.currAxes = 3;	set(handles.figure1,'CurrentAxes',handles.axes3);
		localImhist(handles,img(:,:,3));
		handles = guidata(handles.figure1);			% The handles has been saved in plot_results
		handles.currAxes = 2;	set(handles.figure1,'CurrentAxes',handles.axes2);
		localImhist(handles,img(:,:,2));
		handles = guidata(handles.figure1);			% The handles has been saved in plot_results
		handles.currAxes = 1;	set(handles.figure1,'CurrentAxes',handles.axes1);
		localImhist(handles,img(:,:,1));
		handles = guidata(handles.figure1);			% The handles has been saved in plot_results
		handles.minCData(1) = double(min(min(img(:,:,1))));
		handles.maxCData(1) = double(max(max(img(:,:,1))));
		handles.minCData(2) = double(min(min(img(:,:,2))));
		handles.maxCData(2) = double(max(max(img(:,:,2))));
		handles.minCData(3) = double(min(min(img(:,:,3))));
		handles.maxCData(3) = double(max(max(img(:,:,3))));
		handles.satistic(2:3) = {' '};
	end
	handles.orig_img = img;             % Store a copy of the original image

	set(handles.edit_minRange,'String',num2str(double(handles.minCData(1))))
	set(handles.edit_maxRange,'String',num2str(double(handles.maxCData(1))))

	if (nargout),   varargout{1} = hObject;     end
	guidata(hObject, handles);

	set(hObject,'Vis','on');

% --------------------------------------------------------------------------
function edit_minWindow_CB(hObject, handles)
	x_min = str2double(get(hObject,'String'));
	x_max = str2double(get(handles.edit_maxWindow,'String'));
	if (isnan(x_min) || x_min < 0 || x_min > x_max),   set(hObject,'String','');    return;  end
	updateAll(handles,[x_min x_max],'final')

% --------------------------------------------------------------------------
function edit_maxWindow_CB(hObject, handles)
	x_max = str2double(get(hObject,'String'));
	x_min = str2double(get(handles.edit_minWindow,'String'));
	if (isnan(x_max) || x_min > x_max),   set(hObject,'String',''); return;  end
	updateAll(handles,[x_min x_max],'final')

% --------------------------------------------------------------------------
function edit_widthWindow_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx) || xx < 1),   set(hObject,'String',''); return;  end
	updateAll(handles,[xx xx],'width')

% --------------------------------------------------------------------------
function edit_centerWindow_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx) || xx < 1),   set(hObject,'String',''); return;  end
	updateAll(handles,[xx xx],'center')

% --------------------------------------------------------------------------
function edit_percentOutliersR_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx) || xx < 0 || xx > 100),   set(hObject,'String','2');  end

% --------------------------------------------------------------------------
function edit_percentOutliersG_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx) || xx < 0 || xx > 100),   set(hObject,'String','2');  end

% --------------------------------------------------------------------------
function edit_percentOutliersB_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx) || xx < 0 || xx > 100),   set(hObject,'String','2');  end

% --------------------------------------------------------------------------
function radio_matchRange_CB(hObject, handles)
	if (get(hObject,'Value'))
		set(handles.radio_delOutliers,'Value',0)
	else
		set(hObject,'Value',1)
	end

% --------------------------------------------------------------------------
function radio_delOutliers_CB(hObject, handles)
	if (get(hObject,'Value'))
		set(handles.radio_matchRange,'Value',0)
	else
		set(hObject,'Value',1)
	end

% --------------------------------------------------------------------------
function radio_R_CB(hObject, handles)
	if (get(hObject,'Value')),	swap_radios(handles,1)
	else						set(hObject, 'Value',1)
	end

% --------------------------------------------------------------------------
function radio_G_CB(hObject, handles)
	if (get(hObject,'Value')),	swap_radios(handles,2)
	else						set(hObject, 'Value',1)
	end

% --------------------------------------------------------------------------
function radio_B_CB(hObject, handles)
	if (get(hObject,'Value')),	swap_radios(handles,3)
	else						set(hObject, 'Value',1)
	end

% --------------------------------------------------------------------------
function push_applyRange_CB(hObject, handles)
if (get(handles.radio_delOutliers,'Value'))
    if (~handles.isRGB)
    	outlierPct = str2double(get(handles.edit_percentOutliersR,'String')) / 100;
    	newClim = localStretchlim(handles.orig_img, outlierPct / 2);
    else
        img = handles.orig_img;
        switch handles.currAxes
            case 1
                outlierPct = str2double(get(handles.edit_percentOutliersR,'String')) / 100;
                img(:,:,2:3) = [];
            case 2
                outlierPct = str2double(get(handles.edit_percentOutliersG,'String')) / 100;
                img(:,:,[1 3]) = [];
            case 3
                outlierPct = str2double(get(handles.edit_percentOutliersB,'String')) / 100;
                img(:,:,1:2) = [];
        end
    	newClim = localStretchlim(img, outlierPct / 2);
    end
    if isequal(newClim, [0;1])
        if outlierPct > 0.02
            errordlg('The specified percentage is too great. Coose a smaller one','Percentage Too Large')
        elseif outlierPct ~= 0
            errordlg('This image contains too few grayscale values to eliminate outliers.','Cannot Eliminate Outliers')
        end
        return
    end
    updateAll(handles,newClim,'convert')
else                            % Use pre-stored data range
    newClim = [handles.minCData(handles.currAxes) handles.maxCData(handles.currAxes)];
    updateAll(handles,newClim,'final')
end

% --------------------------------------------------------------------------
function push_contStrectch_CB(hObject, handles)
	limsR = get(handles.patch(1),'XData');
	if (handles.isRGB)
		limsG = get(handles.patch(2),'XData');
		limsB = get(handles.patch(3),'XData');
		low_high = [limsR(1) limsG(1) limsB(1); limsR(4) limsG(4) limsB(4)] / 255;
	else
		low_high = [limsR(1) limsR(4)] / 255;
	end
	set(handles.figure1,'pointer','watch')
	img = img_fun('imadjust_j',handles.orig_img,low_high,[]);
	set(handles.hMirImg, 'Cdata', img)
	set(handles.figure1,'pointer','arrow')

% --------------------------------------------------------------------------
function push_decorrStrectch_CB(hObject, handles)
	set(handles.figure1,'pointer','watch')
	tol = str2double(get(handles.edit_percentOutliersR,'String')) / 100 / 2;
	img = img_fun('decorrstretch',handles.orig_img, 'Tol', tol);
	set(handles.hMirImg, 'Cdata', img)
	set(handles.figure1,'pointer','arrow')

% --------------------------------------------------------------------------
function push_scaterPlot_CB(hObject, handles)
[m,n,k] = size(handles.orig_img);
if (m*n > 256*256)      % For scater plots we don't need all points (If they are many)
    fac = min(256/m,256/n);
    img = localImresize(get(handles.hMirImg,'CData'),fac);
    r = img(:,:,1);     g = img(:,:,2);     b = img(:,:,3);
else
    r = handles.orig_img(:,:,1);     g = handles.orig_img(:,:,2);     b = handles.orig_img(:,:,3);
end
h = figure;
plot3(r(:),g(:),b(:),'.')
grid('on');     xlabel('Red Band');     ylabel('Green Band');   zlabel('Blue Band')

% --------------------------------------------------------------------------
function img = localImresize(img,m)
% Specialized version of IMRESIZE to reduce the size of IMG by factor M
% Define old and new image sizes, and actual scaling
	[so(1),so(2),thirdD] = size(img);   % old image size
	sn = max(floor(m*so(1:2)),1);       % new image size=(integer>0)
	sc = so./sn + eps;                  % actual scaling between images, eps is needed
	img = img(floor(sc(1)/2:sc(1):end)+1,floor(sc(2)/2:sc(2):end)+1,:);

% --------------------------------------------------------------------------
function [y,x] = localImhist(handles,varargin)

	[a, n, isScaled, top] = parse_inputs(varargin{:});

	if islogical(a)
		y(2) = sum(a(:));
		y(1) = numel(a) - y(2);
		y = y';
	else
		y = imhistc(a, n, isScaled, top);       % Call MEX file to do work.
	end

	classin = class(a);
	switch classin
		case 'uint8',   x = linspace(0,255,n)';
		case 'uint16',  x = linspace(0,65535,n)';
		case 'double',  x = linspace(0,1,n)';
		case 'logical', x = [0,1]';
	otherwise
		errordlg('The input image must be uint8, uint16, double, or logical.','ERROR');
		return
	end

	if (nargout == 0),    plot_result(x, y, handles);     end

% --------------------------------------------------------------------------
function [a, n, isScaled, top] = parse_inputs(varargin)

	a = varargin{1};        n = 256;

	if (isa(a,'double'))
		isScaled = 1;   top = 1;
	elseif (isa(a,'uint8'))
		isScaled = 1;   top = 255;
	elseif (islogical(a))
		n = 2;      isScaled = 1;    top = 1;
	else % a must be uint16
		isScaled = 1;   top = 65535;
	end

	if (nargin ==2)
		if (numel(varargin{2}) == 1)        % IMHIST(I, N)
			n = varargin{2};
		elseif (size(varargin{2},2) == 3)   % IMHIST(X,MAP) or invalid second argument
			n = size(varargin{2},1);
			isScaled = 0;
			top = n;
		else
			messageId = sprintf('Images:%s:invalidSecondArgument', mfilename);
			error(messageId, '%s', 'Second argument must be a colormap or a positive integer.'); 
		end  
	end

% --------------------------------------------------------------------------
function plot_result(x, y, handles)

	hAxes = get(handles.figure1,'CurrentAxes');
	handles.histo{handles.currAxes} = [x y];      handles.n_tot(handles.currAxes) = sum(y);     % And the no-value pixels?

	% Create a patch that encloses (x,y) limits
	x_min = min(x);     x_max = max(x);     y_max = max(y);
	if (handles.isRGB)
		switch handles.currAxes
			case 1,     cor = [255 124 117] / 255;
			case 2,     cor = [120 255 114] / 255;
			case 3,     cor = [119 119 255] / 255;
		end
	else
		cor = [120 255 114] / 255;
	end
	handles.patch(handles.currAxes) = patch([x_min x_min x_max x_max],[0 y_max y_max 0],cor);

	localStem(x,y)
	limits = [get(hAxes,'XLim') get(hAxes,'YLim')];
	limits(1) = -10;    limits(2) = 266;
	set(hAxes,'XLim',limits(1:2),'YLim',limits(3:4),'XLimMode','manual','YLimMode','manual');

	if (handles.isRGB && handles.currAxes < 3)  % We don't want XTicks on R & G axes of a RGB display
		set(hAxes, 'Xtick',[])
		label = get(hAxes,'YTickLabel');
		label(1) = ' ';     % We also don't want the '0' because is overlaps the next axes YTickLabel(end,:)
		set(hAxes,'YTickLabel',label)
	end
	set(handles.patch(handles.currAxes),'YData',[0 limits(4) limits(4) 0])    % Needed to update y_max

	% Create three vertical lines
	h_min = line('XData',[x_min x_min],'YData',[0 limits(4)],'color','r','LineWidth',2);
	set(h_min,'UserData',1)     % Flag to indicate left line
	h_med = line('XData',[(x_min+x_max)/2 (x_min+x_max)/2],'YData',[0 limits(4)],'color','r','LineWidth',2,'LineStyle','--');
	set(h_med,'UserData',2)     % Flag to indicate middle line
	h_max = line('XData',[x_max x_max],'YData',[0 limits(4)],'color','r','LineWidth',2);
	set(h_max,'UserData',3)     % Flag to indicate right line
	handles.h_vert_lines{handles.currAxes} = [h_min h_med h_max];
	handles.min_max{handles.currAxes} = [x_min x_max];        % In fact this is always [0 255]
	set(handles.figure1,'WindowButtonMotionFcn',{@wbm_vertLine,handles})

	% OK, now fill the edit boxes that are still empty
	set(handles.edit_minWindow,'String',sprintf('%d',round(x_min)))
	set(handles.edit_maxWindow,'String',sprintf('%d',round(x_max)))
	set(handles.edit_widthWindow,'String',sprintf('%d',round(x_max-x_min)))
	set(handles.edit_centerWindow,'String',sprintf('%d',round((x_min+x_max)/2)))

	handles.minWindow(handles.currAxes) = round(x_min);
	handles.maxWindow(handles.currAxes) = round(x_max);
	handles.widthWindow(handles.currAxes) = round(x_max-x_min);
	handles.centerWindow(handles.currAxes) = round((x_min+x_max)/2);

	guidata(handles.figure1,handles)

% --------------------------------------------------------------------------
function localStem(x,y)
%   This is a hacked code from the stem function
%   Only strictly necessary code was kept

	% Set up data using fancing indexing
	[m,n] = size(x);
	xx = zeros(3*m,n);      yy = xx;
	xx(1:3:3*m,:) = x;      yy(2:3:3*m,:) = y;
	xx(2:3:3*m,:) = x;      yy(3:3:3*m,:) = NaN;
	xx(3:3:3*m,:) = NaN;

	line('XData',xx,'YData',yy,'color','k');

% --------------------------------------------------------------------------
function lowhigh = localStretchlim(img,tol)
%STRETCHLIM Find limits to contrast stretch an image.
%   LOW_HIGH = STRETCHLIM(I,TOL) returns a pair of gray values that can be
%   used by IMADJUST to increase the contrast of an image.

	if (nargin == 2)
	if (numel(tol) == 1),        tol(2) = 1 - tol;    end
	if (tol(1) >= tol(2))
		errordlg('ERROR in strechlim: TOL(1) must be less than TOL(2).','ERROR');
	end
	if ( any(tol < 0) || any(tol > 1) || any(isnan(tol)) )
		errordlg('ERROR in strechlim: TOL must be in the range [0 1].','ERROR');
	end
	else
	tol = [.01 .99];
	end

	if isa(img,'uint8'),    nbins = 256;
	else                    nbins = 65536;
	end

	tol_low = tol(1);       tol_high = tol(2); 
	p = size(img,3);

	if (tol_low < tol_high)
	ilowhigh = zeros(2,p);
	for i = 1:p                         % Find limits, one plane at a time
		N = localImhist([],img(:,:,i),nbins);   % The [] is instead of handles (not used here)
		cdf = cumsum(N)/sum(N);         %cumulative distribution function
		ilow = min(find(cdf > tol_low));
		ihigh = min(find(cdf >= tol_high));
		ilowhigh(:,i) = [ilow;ihigh];
		if (ilow == ihigh)
			ilowhigh(:,i) = [1; nbins]; % limits are same, use default   
		end
	end
	lowhigh = (ilowhigh - 1)/(nbins-1);  % convert to range [0 1]
	else
	%   tol_low >= tol_high, this tolerance does not make sense. For example, if
	%   the tolerance is .5 then no pixels would be left after saturating
	%   low and high pixel values. In all of these cases, STRETCHLIM returns [0; 1].
	lowhigh = repmat([0;1],1,p);
	end

% --------------------------------------------------------------------------
function updateAll(handles,newClim, opt)
% This function is called by the edit boxes and the Apply button
	if (strcmp(opt,'convert'))
		newClim = double(loc_intmax(handles.imgClass)) * newClim;
	elseif (strcmp(opt,'width'))
		x_center = get(handles.h_vert_lines{handles.currAxes}(2),'XData');
		newClim = [x_center-newClim(1)/2 x_center+newClim(1)/2];
		if (newClim(1) < 0),    newClim(1) = 0;     end
		if (newClim(2) > handles.min_max{handles.currAxes}(2)),   newClim(2) = handles.min_max{handles.currAxes}(2);   end
	elseif (strcmp(opt,'center'))
		xx = get(handles.patch,'XData');    halfDx = (xx(4) - xx(1)) / 2;
		newClim = [newClim(1)-halfDx newClim(1)+halfDx];
		if (newClim(1) < 0),    newClim(1) = 0;     end
		if (newClim(2) > handles.min_max{handles.currAxes}(2)),   newClim(2) = handles.min_max{handles.currAxes}(2);   end    
	elseif (strcmp(opt,'final'))
		% Nothing to do. Everything is recomputed below
	end

	set(handles.edit_minRange,'String',num2str(double(handles.minCData(handles.currAxes))))
	set(handles.edit_maxRange,'String',num2str(double(handles.maxCData(handles.currAxes))))
	set(handles.edit_minWindow,'String',round(newClim(1)))
	set(handles.edit_maxWindow,'String',round(newClim(2)))
	set(handles.edit_widthWindow,'String',round(newClim(2)-newClim(1)))
	set(handles.edit_centerWindow,'String',round((newClim(1)+newClim(2))/2))
	handles.minWindow(handles.currAxes) = round(newClim(1));
	handles.maxWindow(handles.currAxes) = round(newClim(2));
	handles.widthWindow(handles.currAxes) = round(newClim(2)-newClim(1));
	handles.centerWindow(handles.currAxes) = round((newClim(1)+newClim(2))/2);

	if (~handles.isRGB),    set(handles.hMirAxes,'CLim',newClim);  end
	set(handles.patch(handles.currAxes),'XData',[newClim(1) newClim(1) newClim(2) newClim(2)]);
	set(handles.h_vert_lines{handles.currAxes}(1),'XData',[newClim(1) newClim(1)])
	set(handles.h_vert_lines{handles.currAxes}(2),'XData',[(newClim(1)+newClim(2))/2 (newClim(1)+newClim(2))/2])
	set(handles.h_vert_lines{handles.currAxes}(3),'XData',[newClim(2) newClim(2)])
	%n_pixels_in_patch = sum(handles.histo{handles.currAxes}(newClim(1):newClim(2),2));
	patch_lims = get(handles.patch(handles.currAxes),'XData');
	idx = (handles.histo{handles.currAxes}(:,1) >= patch_lims(1)) & (handles.histo{handles.currAxes}(:,1) <= patch_lims(4));
	n_pixels_in_patch = sum(handles.histo{handles.currAxes}(idx,2));
	str = sprintf('N = %d\t(%.2f%%)',[n_pixels_in_patch n_pixels_in_patch/handles.n_tot(handles.currAxes)*100]);
	set(handles.text_stat,'String',str)
	handles.satistic{handles.currAxes} = str;
	guidata(handles.figure1,handles)

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
			error('image_enhance:loc_intmax:invalidClassName','Invalid class name.')
	end

% --------------------------------------------------------------------------
function wbm_vertLine(obj,eventdata,handles)
	handles = guidata(handles.figure1);

	hAxes = get(handles.figure1,'CurrentAxes');
	pt = get(hAxes, 'CurrentPoint');        lim = [get(hAxes,'XLim') get(hAxes,'YLim')];
	if (pt(1,1) < lim(1)) || (pt(1,1) > lim(2)) || (pt(1,2) < lim(3)) || (pt(1,2) > lim(4));
		set(handles.figure1,'Pointer','arrow','WindowButtonDownFcn',{@wbd_strayClick,handles})
		return
	end
	xx = get(handles.h_vert_lines{handles.currAxes},'XData');
	if ( abs(pt(1,1) - xx{1}(1)) < 2)
		set(handles.figure1,'Pointer','custom','PointerShapeCData',handles.new_pointer,'PointerShapeHotSpot',[8 8], ...
		'WindowButtonDownFcn',{@wbd_vertLine,handles.h_vert_lines{handles.currAxes}(1),handles},'WindowButtonUpFcn',{@wbu_vertLine,handles})
	elseif ( abs(pt(1,1) - xx{2}(1)) < 2 )
		set(handles.figure1,'Pointer','custom','PointerShapeCData',handles.new_pointer,'PointerShapeHotSpot',[8 8], ...
		'WindowButtonDownFcn',{@wbd_vertLine,handles.h_vert_lines{handles.currAxes}(2),handles}, ...
		'WindowButtonUpFcn',{@wbu_vertLine,handles})
	elseif ( abs(pt(1,1) - xx{3}(1)) < 2 )
		set(handles.figure1,'Pointer','custom','PointerShapeCData',handles.new_pointer,'PointerShapeHotSpot',[8 8], ...
		'WindowButtonDownFcn',{@wbd_vertLine,handles.h_vert_lines{handles.currAxes}(3),handles}, ...
		'WindowButtonUpFcn',{@wbu_vertLine,handles})
	else
		set(handles.figure1,'Pointer','arrow','WindowButtonDownFcn','','WindowButtonUpFcn',{@wbd_strayClick,handles})
	end

% -----------
function wbd_strayClick(obj,eventdata,handles)
% When the user clicks on a axes it automatically will become the active one. Register that.
	handles = guidata(handles.figure1);
	tag = get(get(handles.figure1,'CurrentAxes'),'Tag');
	what_axes = str2double(tag(end));
	if (handles.isRGB)      % Update Band radiotuttons and handles.currAxes
		swap_radios(handles,what_axes)
	end

% -----------
function wbd_vertLine(obj,eventdata,h,handles)
	n = get(h,'UserData');
	set(handles.figure1,'WindowButtonMotionFcn',{@drag_vertLine,h,handles,n});

% -----------
function drag_vertLine(obj,eventdata,h,handles,n)
%pt = get(handles.(['axes' sprintf('%d',handles.currAxes)]), 'CurrentPoint');    x = pt(1,1);
	pt = get(get(handles.figure1,'CurrentAxes'), 'CurrentPoint');    x = pt(1,1);
	xp = get(handles.patch(handles.currAxes),'XData');
	if (x < handles.min_max{handles.currAxes}(1)-0.4 || x > handles.min_max{handles.currAxes}(2)),  return;     end
	set(h,'XData',[x x]);
	% OK, now adapt the patch accordingly
	switch n
		case 1				% Update the left side
			x = (xp(1) + xp(4)) / 2;
			new_patch_lims = [round(pt(1,1)) round(pt(1,1)) xp(3) xp(4)];
			set(handles.h_vert_lines{handles.currAxes}(2),'XData',[x x])    % Update central line
			set(handles.edit_minWindow,'String',round(pt(1,1)))
			set(handles.edit_widthWindow,'String',round(xp(4)-pt(1,1)))
			handles.minWindow(handles.currAxes) = round(pt(1,1));
			handles.widthWindow(handles.currAxes) = round(xp(4)-pt(1,1));
		case 2				% Update both left and right sides
			dx = xp(4) - xp(1);
			xp(1) = max(handles.min_max{handles.currAxes}(1),pt(1,1)-dx/2);
			xp(4) = min(handles.min_max{handles.currAxes}(2),pt(1,1)+dx/2);
			new_patch_lims = [xp(1) xp(1) xp(4) xp(4)];
			set(handles.h_vert_lines{handles.currAxes}(1),'XData',[xp(1) xp(1)])    % Update left line
			set(handles.h_vert_lines{handles.currAxes}(3),'XData',[xp(4) xp(4)])    % Update right line
			set(handles.edit_minWindow,'String',round(xp(1)))
			set(handles.edit_maxWindow,'String',round(xp(4)))
			set(handles.edit_widthWindow,'String',round(xp(4)-xp(1)))
			handles.minWindow(handles.currAxes) = round(xp(1));
			handles.maxWindow(handles.currAxes) = round(xp(4));
			handles.widthWindow(handles.currAxes) = round(xp(4)-xp(1));
		case 3				% Update the right side
			x = (xp(1) + xp(4)) / 2;
			new_patch_lims = [xp(1) xp(2) round(pt(1,1)) round(pt(1,1))];
			set(handles.h_vert_lines{handles.currAxes}(2),'XData',[x x])    % Update central line
			set(handles.edit_maxWindow,'String',round(pt(1,1)))
			set(handles.edit_widthWindow,'String',round(pt(1,1)-xp(1)))
			handles.maxWindow(handles.currAxes) = round(pt(1,1));
			handles.widthWindow(handles.currAxes) = round(pt(1,1)-xp(1));
	end

	handles.centerWindow(handles.currAxes) = round(x);
	set(handles.hMirAxes,'CLim',[new_patch_lims(1) new_patch_lims(4)])
	set(handles.patch(handles.currAxes),'XData',new_patch_lims)
	set(handles.edit_centerWindow,'String',round(x))
	idx = (handles.histo{handles.currAxes}(:,1) >= new_patch_lims(1)) & (handles.histo{handles.currAxes}(:,1) <= new_patch_lims(4));
	n_pixels_in_patch = sum(handles.histo{handles.currAxes}(idx,2));
	str = sprintf('N = %d\t(%.2f%%)',[n_pixels_in_patch n_pixels_in_patch/handles.n_tot(handles.currAxes)*100]);
	set(handles.text_stat,'String',str)
	guidata(handles.figure1,handles)

% -----------
function wbu_vertLine(obj,eventdata,handles)
	set(handles.figure1,'Pointer','arrow', 'WindowButtonUpFcn','', 'WindowButtonDownFcn','', ...
	'WindowButtonMotionFcn',{@wbm_vertLine,handles});

	xp = get(handles.patch(handles.currAxes),'XData');
	set(handles.hMirAxes,'CLim',[xp(1) xp(4)])

% --------------------------------------------------------------------------
function swap_radios(handles,n)
	switch n
		case 1
			set([handles.radio_G handles.radio_B], 'Value',0)
			set(handles.radio_R,'Value',1)     % Needed when called from wbd_strayClick
			axes(handles.axes1)
		case 2
			set([handles.radio_R handles.radio_B], 'Value',0)
			set(handles.radio_G,'Value',1)     % Needed when called from wbd_strayClick
			axes(handles.axes2)
		case 3
			set([handles.radio_R handles.radio_G], 'Value',0)
			set(handles.radio_B,'Value',1)     % Needed when called from wbd_strayClick
			axes(handles.axes3)
	end
	handles.currAxes = n;

	if (handles.isRGB)      % Update also the edit boxes
		set(handles.edit_minWindow,'String',handles.minWindow(n))
		set(handles.edit_maxWindow,'String',handles.maxWindow(n))
		set(handles.edit_widthWindow,'String',handles.widthWindow(n))
		set(handles.edit_centerWindow,'String',handles.centerWindow(n))
		set(handles.edit_minRange,'String',handles.minCData(n))
		set(handles.edit_maxRange,'String',handles.maxCData(n))
		set(handles.text_stat,'String',handles.satistic{n})
	end

	guidata(handles.figure1,handles)


% --- Creates and returns a handle to the GUI figure. 
function image_enhance_LayoutFcn(h1)

set(h1, 'Color',get(0,'factoryUicontrolBackgroundColor'),...
'DoubleBuffer','on',...
'MenuBar','none',...
'Name','Adjust Contrast',...
'NumberTitle','off',...
'Position',[520 293 536 540],...
'Resize','off',...
'Tag','figure1');

uicontrol('Parent',h1, 'Position',[5 29 526 391], 'Style','frame', 'Tag','frame_axes');
uicontrol('Parent',h1, 'Position',[5 422 526 116], 'Style','frame');
uicontrol('Parent',h1, 'Position',[360 432 161 85], 'Style','frame');
uicontrol('Parent',h1, 'Position',[129 462 222 55], 'Style','frame');
uicontrol('Parent',h1, 'Position',[10 462 111 55], 'Style','frame');

uicontrol('Parent',h1, 'BackgroundColor',[0.753 0.753 0.753],...
'Enable','inactive', 'HorizontalAlignment','center',...
'Position',[66 489 46 20], 'Style','edit', 'Tag','edit_minRange');

uicontrol('Parent',h1, 'HorizontalAlignment','right',...
'Position',[5 491 60 15], 'String','Minimum:', 'Style','text');

uicontrol('Parent',h1, 'BackgroundColor',[0.753 0.753 0.753],...
'Enable','inactive', 'HorizontalAlignment','center',...
'Position',[66 467 46 20], 'Style','edit', 'Tag','edit_maxRange');

uicontrol('Parent',h1, 'HorizontalAlignment','right',...
'Position',[5 472 60 15], 'String','Maximum:', 'Style','text');

uicontrol('Parent',h1, 'BackgroundColor',[1 1 1],...
'Call',@main_uiCB,...
'HorizontalAlignment','center', 'Position',[186 489 47 20],...
'Style','edit', 'Tag','edit_minWindow');

uicontrol('Parent',h1, 'HorizontalAlignment','right',...
'Position',[135 491 50 15], 'String','Minimum:', 'Style','text');

uicontrol('Parent',h1, 'BackgroundColor',[1 1 1],...
'Call',@main_uiCB,...
'HorizontalAlignment','center','Position',[186 467 47 20],...
'Style','edit', 'Tag','edit_maxWindow');

uicontrol('Parent',h1, 'HorizontalAlignment','right',...
'Position',[130 472 55 15], 'String','Maximum:', 'Style','text');

uicontrol('Parent',h1, 'BackgroundColor',[1 1 1],...
'Call',@main_uiCB,...
'HorizontalAlignment','center','Position',[294 488 47 20],...
'Style','edit', 'Tag','edit_widthWindow');

uicontrol('Parent',h1, 'HorizontalAlignment','right',...
'Position',[258 490 35 15], 'String','Width:', 'Style','text');

uicontrol('Parent',h1, 'BackgroundColor',[1 1 1],...
'Call',@main_uiCB,...
'HorizontalAlignment','center','Position',[294 466 47 20],...
'Style','edit', 'Tag','edit_centerWindow');

uicontrol('Parent',h1, 'HorizontalAlignment','right',...
'Position',[252 469 40 15], 'String','Center:', 'Style','text');

uicontrol('Parent',h1, 'BackgroundColor',[1 1 1],...
'Call',@main_uiCB,...
'Position',[474 470 25 18],'String','2', 'Style','edit',...
'Tooltip','Specifies the fraction of the image (or R band) to saturate at low and high intensities.',...
'Tag','edit_percentOutliersR');

uicontrol('Parent',h1,...
'Call',@main_uiCB,...
'Position',[367 486 130 17],'String','Match Data Range',...
'Style','radiobutton',...
'Tooltip','Try this it see what it does',...
'Value',1, 'Tag','radio_matchRange');

uicontrol('Parent',h1,...
'Call',@main_uiCB,...
'Position',[367 465 107 17],'String','Eliminate outliers',...
'Style','radiobutton',...
'Tooltip','Eleminates % outliers as selected in the right side boxe(s)',...
'Tag','radio_delOutliers');

uicontrol('Parent',h1, 'HorizontalAlignment','left',...
'Position',[502 470 10 15],'String','%', 'Style','text');

uicontrol('Parent',h1, 'Position',[24 509 70 15],...
'String','Data Range', 'Style','text','Tag','text_dataRange');

uicontrol('Parent',h1, 'Position',[153 509 50 15],...
'String','Window', 'Style','text', 'Tag','text_window');

uicontrol('Parent',h1,...
'Call',@main_uiCB,...
'FontName','Helvetica', 'FontSize',9,...
'Position',[369 438 80 21],...
'String','Apply', 'Tag','push_applyRange');

uicontrol('Parent',h1, 'Position',[370 509 115 15],...
'String','Scale Display Range', 'Style','text', 'Tag','text_scaleRange');

axes('Parent',h1, 'Units','pixels',...
'Color',[0.9 1 1], 'Position',[50 291 431 115], 'Tag','axes1');

uicontrol('Parent',h1, 'HorizontalAlignment','left','Position',[6 7 221 15],...
'String','Adjust the histogram above, or click and drag',...
'Style','text', 'Tag','text_tip');

uicontrol('Parent',h1, 'HorizontalAlignment','left',...
'Position',[398 7 132 15], 'Style','text',...
'Tooltip','Statistics of the selected histogram portion',...
'Tag','text_stat');

axes('Parent',h1, 'Units','pixels',...
'Color',[0.9 1 1], 'Position',[50 173 431 115], 'Tag','axes2');

axes('Parent',h1, 'Units','pixels',...
'Color',[0.9 1 1], 'Position',[50 54 431 115], 'Tag','axes3');

uicontrol('Parent',h1,...
'Call',@main_uiCB,...
'FontName','Helvetica', 'FontSize',9,'Position',[16 432 110 20],...
'String','Contrast Stretch',...
'Tooltip','Increase the image contrast  with parameters deduced from the "colored" histograms portion',...
'Tag','push_contStrectch');

uicontrol('Parent',h1,...
'Call',@main_uiCB,...
'FontName','Helvetica', 'FontSize',9,'Position',[132 432 130 20],...
'String','Decorrelation Stretch',...
'Tooltip','Do a decorrelation stretch using as tolerance the value of the R band outliers percentage',...
'Tag','push_decorrStrectch');

uicontrol('Parent',h1, 'BackgroundColor',[0 1 0],...
'Call',@main_uiCB,...
'Position',[474 453 25 18], 'String','2', 'Style','edit',...
'Tooltip','Specifies the fraction of the Green band to saturate at low and high intensities.',...
'Tag','edit_percentOutliersG');

uicontrol('Parent',h1, 'BackgroundColor',[0 0 1],...
'Call',@main_uiCB,...
'Position',[474 436 25 18], 'String','2', 'Style','edit',...
'Tooltip','Specifies the fraction of the Blue band to saturate at low and high intensities.',...
'Tag','edit_percentOutliersB');

uicontrol('Parent',h1,...
'Call',@main_uiCB,...
'FontName','Helvetica','Position',[489 336 35 15],...
'String','R','Style','radiobutton','Value',1,'Tag','radio_R');

uicontrol('Parent',h1,...
'Call',@main_uiCB,...
'FontName','Helvetica','Position',[489 223 35 15],...
'String','G','Style','radiobutton','Tag','radio_G');

uicontrol('Parent',h1,...
'Call',@main_uiCB,...
'FontName','Helvetica','Position',[489 105 35 15],...
'String','B','Style','radiobutton','Tag','radio_B');

uicontrol('Parent',h1,...
'Call',@main_uiCB,...
'FontName','Helvetica', 'FontSize',9,...
'Position',[268 432 80 20],...
'String','ScaterPlot',...
'Tooltip','Sow a scater plot of the 3 bands',...
'Tag','push_scaterPlot');

function main_uiCB(hObject, eventdata)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));

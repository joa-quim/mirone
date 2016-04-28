function varargout = thresholdit(varargin)
% Helper window to binarize an image by threshold computation.

%	Copyright (c) 2004-2014 by J. Luis
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

% $Id: thresholdit.m 4285 2014-01-17 04:04:39Z j $

	if (isempty(varargin)),		return,		end

	hObject = figure('Vis','off');
	thresholdit_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hObject,'right')

	handMir = guidata(varargin{1});
	handles.hMirFig = varargin{1};
	handles.hImgMir = handMir.hImg;

	if (handMir.no_file)
		errordlg('You didn''t even load a file. What are you expecting then?','Error')
		delete(hObject);    return
	end

	handles.new_pointer = [ ...
			NaN	NaN	NaN	NaN	NaN	NaN	NaN	2	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN; ...
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

	handles.bg_color = uint8(handMir.bg_color * 255);
	handles.I_never_changed_imgMir = true;
	handles.img = get(handMir.hImg,'CData');

	% Must find out if Mirone's backup img is the same as the displayed (to eventual use is "Apply to original"
	handles.fine_to_use_imgMir_orig = true;			% Start by assuming a yes
	handles.fine_to_use_imgMir_orig = (numel(handles.img) == numel(handMir.origFig));	% simple test on images size
	if (handles.fine_to_use_imgMir_orig)			% make a further test of equality on 10 randomly picked points
		rands = fix( rand(1,10) * (size(handles.img,1) * size(handles.img,2) - 1) ) + 1;
		handles.fine_to_use_imgMir_orig = isequal(handles.img(rands), handMir.origFig(rands));
	end

	if (ndims(handles.img) == 3)
		handles.img = cvlib_mex('color',handles.img,'rgb2gray');
	end

	[bw,level] = img_fun('im2bw',handles.img);
	handles.hImg = image(bw,'Parent',handles.axes1);
	set(handles.hImg,'CDataMapping','scaled');
	set(handles.figure1,'ColorMap',gray(2));

	[m,n] = size(handles.img);
	set(handles.axes1,'PlotBoxAspectRatio',[1 m/n 1],'Visible','off')
	if (strcmp(get(handMir.axes1,'YDir'),'normal'))
		set(handles.axes1,'YDir','normal')
	end

	% ---------------- Construct the histogram
	localImhist(handles.img, handles);

	% ---------------- display level as vertical line
	x_lim = get(handles.axes2,'xlim');		y_lim = get(handles.axes2,'ylim');
	handles.hVline = line('XData',[level level], 'YData',y_lim, 'Parent',handles.axes2,...
		'LineWidth',2,'color',0.5*[1 1 1],'handlevisibility','off');

	handles.hText = text(level,-1,sprintf('%d',level),'HorizontalAlignment','Center',...
		'VerticalAlignment','top','FontWeight','Bold','color',.4*[1 1 1]);
	movex_text(handles,level)

	% ---------------- Create one box plus three vertical lines
	handles.hBox = patch('XData',[0 0 x_lim(2) x_lim(2)],'YData',[0 y_lim(2) y_lim(2) 0], ...
		'FaceColor','none', 'Parent',handles.axes2, 'Vis','off');
	handles.hVertLines(1) = line('XData',[0 0],'YData',[0 y_lim(2)],'color','r','LineWidth',2, 'UserData',1, 'Vis','off');
	handles.hVertLines(2) = line('XData',[x_lim(2) x_lim(2)]/2,'YData',[0 y_lim(2)],'color','r', ...
	'LineWidth',2,'LineStyle','--', 'UserData',2, 'Vis','off');
	handles.hVertLines(3) = line('XData',[x_lim(2) x_lim(2)],'YData',[0 y_lim(2)],'color','r', ...
		'LineWidth',2, 'UserData',3, 'Vis','off');

	handles = move_vline(handles);			% attach draggable behavior for user to change level

	% Add this figure handle to the carraças list
	plugedWin = getappdata(handles.hMirFig,'dependentFigs');
	plugedWin = [plugedWin hObject];
	setappdata(handles.hMirFig,'dependentFigs',plugedWin);

	set(hObject,'Visible','on');
	guidata(hObject, handles);
	if (nargout),	varargout{1} = hObject;		end

% --------------------------------------------------------------------------
function [y,x] = localImhist(img, handles)
	[n, isScaled, top] = parse_inputs(img);
	y = imhistc(img, n, isScaled, top);       % Call MEX file to do work.

	switch class(img)
		case 'uint8',   x = linspace(0,255,n)';
		case 'uint16',  x = linspace(0,65535,n)';
		case 'double',  x = linspace(0,1,n)';
		case 'logical', x = [0,1]';
	otherwise
		errordlg('The input image must be uint8, uint16, double, or logical.','ERROR');
		return
	end

	if (nargout == 0)
		[xx,yy] = localStem(x,y);
		plot(xx,yy, 'Parent',handles.axes2, 'LineWidth',2.5);
		set(handles.axes2,'XLim',[min(x) max(x)],'YLim',[min(y) max(y)], 'Xtick',[min(x) max(x)], 'Ytick',[]);
	end

% --------------------------------------------------------------------------
function [n, isScaled, top] = parse_inputs(img)

	n = 256;
	if (isa(img,'double'))
		isScaled = 1;   top = 1;
	elseif (isa(img,'uint8'))
		isScaled = 1;   top = 255;
	elseif (islogical(img))
		n = 2;      isScaled = 1;    top = 1;
	else		% img must be uint16
		isScaled = 1;   top = 65535;
	end

% --------------------------------------------------------------------------
function [xx,yy] = localStem(x,y)
%   This is a hacked code from the stem function

	% Set up data using fancing indexing
	[m,n] = size(x);
	xx = zeros(3*m,n);      xx(1:3:3*m,:) = x;
	xx(2:3:3*m,:) = x;      xx(3:3:3*m,:) = NaN;

	[m,n] = size(y);        yy = zeros(3*m,n);
	yy(2:3:3*m,:) = y;      yy(3:3:3*m,:) = NaN;

%----------------------------------------------------------------------
function movex_text(handles,x_pos)
% Update the threshold level text position to the same of vertical line
    msg = sprintf('%d',round(x_pos));
    pos = get(handles.hText,'position');
    pos(1) = x_pos;
    set(handles.hText,'Position',pos,'String',msg)

%----------------------------------------------------------------------
function handles = move_vline(handles)
% implements horizontal movement of line(s).

	% This seems to lock the axes position
	set(handles.figure1,'Nextplot','Replace')
	set(handles.hVline,'ButtonDownFcn',{@VLineDownFcn,handles})
	handles.hFunWBM_VL  = @VLineDownFcn;
	handles.hFunWBM_Box = @BoxDownFcn;

% -------------------------------------------------------------
function VLineDownFcn(hObject,evt,handles)
% Set of functions to deal with the manual threshold binarization
	set(handles.figure1,'WindowButtonMotionFcn',{@VLMoveFcn,handles},'WindowButtonUpFcn',{@VLUpFcn,handles})

function VLMoveFcn(hObject,evt,handles)
	cp = get(handles.axes2,'CurrentPoint');                 %
	xpos = cp(1);                                           %
	x_range = get(handles.axes2,'xlim');                    %
	if (xpos < x_range(1)),		xpos = x_range(1);	end     %
	if (xpos > x_range(2)),		xpos = x_range(2);	end     %
	XData = get(handles.hVline,'XData');                    %
	XData(:) = xpos;
	set(handles.hVline,'xdata',XData)
	movex_text(handles,xpos)        % update text

function VLUpFcn(hObject,evt,handles)
	set(handles.figure1,'WindowButtonMotionFcn',[])
	level = mean(get(handles.hVline,'xdata'));
	if (~get(handles.check_revert,'Value'))
		bw = (handles.img > level);
	else
		bw = (handles.img <= level);
	end
	set(handles.hImg,'cdata',bw)
% -------------------------------------------------------------

% -------------------------------------------------------------
function BoxDownFcn(hObject,evt,handles)
% Set of functions to deal with the Kaba binarization
	hLine = gco;
	n = get(hLine,'UserData');
	set(handles.figure1, 'Pointer','custom','PointerShapeCData',handles.new_pointer,'PointerShapeHotSpot',[8 8], ...
		'WindowButtonMotionFcn',{@BoxMoveFcn,handles,hLine,n},'WindowButtonUpFcn',{@BoxUpFcn,handles})

function BoxMoveFcn(obj,evt,handles, hLine, n)
	pt = get(handles.axes2, 'CurrentPoint');
	xp = get(handles.hBox,'XData');
	set(hLine,'XData',[pt(1) pt(1)]);
	% OK, now adapt the patch accordingly
	switch n
		case 1				% Update the left side
			x = (xp(1) + xp(4)) / 2;
			new_patch_lims = [round(pt(1,1)) round(pt(1,1)) xp(3) xp(4)];
			set(handles.hVertLines(2),'XData',[x x])		% Update central line
		case 2				% Update both left and right sides
			dx = xp(4) - xp(1);
			xp(1) = max(0, pt(1,1)-dx/2);
			xp(4) = min(255,pt(1,1)+dx/2);
			new_patch_lims = [xp(1) xp(1) xp(4) xp(4)];
			set(handles.hVertLines(1),'XData',[xp(1) xp(1)])	% Update left line
			set(handles.hVertLines(3),'XData',[xp(4) xp(4)])	% Update right line
		case 3
			x = (xp(1) + xp(4)) / 2;
			new_patch_lims = [xp(1) xp(2) round(pt(1,1)) round(pt(1,1))];
			set(handles.hVertLines(2),'XData',[x x])		% Update central line
	end
	set(handles.hBox,'XData',new_patch_lims)

function BoxUpFcn(obj,evt,handles)
	set(handles.figure1,'Pointer','arrow', 'WindowButtonMotionFcn','');
	x_min = get(handles.hVertLines(1),'xdata');
	x_max = get(handles.hVertLines(3),'xdata');
	if (~get(handles.check_revert,'Value'))
		bw = (handles.img >= x_min(1) & handles.img <= x_max(1));
	else
		bw = (handles.img < x_min(1) | handles.img > x_max(1));
	end
	set(handles.hImg,'cdata',bw)
% -------------------------------------------------------------

%----------------------------------------------------------------------
function push_Otsu_CB(hObject, handles)
% ...
	[bw, threshold] = img_fun('im2bw',handles.img);
	if (get(handles.check_revert,'Value')),		bw = ~bw;	end
	set(handles.hImg,'cdata',bw)
	set(handles.hVline,'xdata',[threshold threshold])
	movex_text(handles, threshold)

%----------------------------------------------------------------------
function push_maxEntropy_CB(hObject, handles)
% ...
	threshold = maxentropy(handles);
	bw = binarize_it(handles.img, threshold, get(handles.check_revert,'Value'));
	set(handles.hImg,'cdata',bw)
	set(handles.hVline,'xdata',[threshold threshold])
	movex_text(handles, threshold)

%----------------------------------------------------------------------
function push_minCrossEntropy_CB(hObject, handles)
% ...
	threshold = minCE(handles);
	bw = binarize_it(handles.img, threshold, get(handles.check_revert,'Value'));
	set(handles.hImg,'cdata',bw)
	set(handles.hVline,'xdata',[threshold threshold])
	movex_text(handles, threshold)

%----------------------------------------------------------------------
function push_isodata_CB(hObject, handles)
	threshold = isodata(handles);
	bw = binarize_it(handles.img, threshold, get(handles.check_revert,'Value'));
	set(handles.hImg,'cdata',bw)
	set(handles.hVline,'xdata',[threshold threshold])
	movex_text(handles, threshold)

%----------------------------------------------------------------------
function push_triang_CB(hObject, handles)
	threshold = triangle_th(handles);
	bw = binarize_it(handles.img, threshold, get(handles.check_revert,'Value'));
	set(handles.hImg,'cdata',bw)
	set(handles.hVline,'xdata',[threshold threshold])
	movex_text(handles, threshold)

%----------------------------------------------------------------------
function push_cleanDust_CB(hObject, handles)
% Clean the isolated small group of black pixels inside white zones or vice-versa
	bw = cleanDust(get(handles.hImg,'cdata'),2);	% Clean white dust over black background
	set(handles.hImg,'cdata',bw)

%----------------------------------------------------------------------
function push_fillHoles_CB(hObject, handles)
% Name says it all
	bw = img_fun('imfill', get(handles.hImg,'cdata'),'holes');
	set(handles.hImg,'cdata',bw)

% --------------------------------------------------------------------------
function bw = cleanDust(bw, opt)
% This function comes from autofaults and I don't remember anymore if it was me who wrote it
	bw2 = img_fun('bwareaopen', bw, 15);	% Remove all connected components that have fewer than 15 pixels
	removed = xor(bw, bw2);
	if (opt == 1)		% Not used here (what does ir do?) so comment for compiling safety
% 		bw3 = imdilate(bw2, strel('disk', 5));
% 		overlaps = bw3 & removed;
% 		bw = (bw2 | overlaps);
	else
		%D = img_fun('bwdist', bw2);
		D = cvlib_mex('distance',~bw2, 2, 4);	% Negate because OpenCV computes distance to zero and ML to non-zero
		within_hailing_distance = (D <= 5);		clear D
		put_back_pixels = removed & within_hailing_distance;	clear removed
		bw = bw2 | put_back_pixels;
	end

%----------------------------------------------------------------------
function radio_singleLine_CB(hObject, handles)
% Make the binarize by threshold active
	if (~get(hObject,'Val')),	set(hObject,'Val',1),	return,		end
	set(handles.radio_window,'Val',0)
	set([handles.hVertLines handles.hBox], 'Vis', 'off')
	set([handles.hVline handles.hText], 'Vis', 'on')
	set(handles.hVline,'ButtonDownFcn',{handles.hFunWBM_VL,handles})

%----------------------------------------------------------------------
function radio_window_CB(hObject, handles)
% Make the binarize by interval (window) active
	if (~get(hObject,'Val')),	set(hObject,'Val',1),	return,		end
	set(handles.radio_singleLine,'Val',0)
	set([handles.hVertLines handles.hBox], 'Vis', 'on')
	set([handles.hVline handles.hText], 'Vis', 'off')
	set(handles.hVertLines,'ButtonDownFcn',{handles.hFunWBM_Box,handles})

%----------------------------------------------------------------------
function check_revert_CB(hObject, handles)
% Compute the image complement and signal to do that while button is checked
	set( handles.hImg,'cdata', ~get(handles.hImg,'cdata') )

%----------------------------------------------------------------------
function push_OK_CB(hObject, handles)
% ...
	handMir = guidata(handles.hMirFig);			% Fish the Mirone handles
	if (get(handles.check_applyOrig,'Val'))
		bw  = get(handles.hImg,'cdata');
		% Here we have to be carefull because the contents of the Mirone fig may have been changed
		% BEFORE calling this tool, or BY this tool and the procedure is different for each case.
		if (handles.I_never_changed_imgMir)
			img = get(handMir.hImg,'cdata');		% Current image displayed in Mirone fig
		elseif (~handles.I_never_changed_imgMir && handles.fine_to_use_imgMir_orig) % Uf, Mirone image backup is fine to use.
			img = handMir.origFig;
		else
			warndlg('Sorry but I have no more copies of the original image that I can change once again. Bye','Warning')
			return
		end
		if (ndims(img) == 3)
			x = img(:,:,1);		x(bw) = handles.bg_color(1);	img(:,:,1) = x;
			x = img(:,:,2);		x(bw) = handles.bg_color(2);	img(:,:,2) = x;
			x = img(:,:,3);		x(bw) = handles.bg_color(3);	img(:,:,3) = x;
		else
			img(bw) = handles.bg_color(1);		% Likely not what the user expects
		end
		set(handles.hImgMir,'cdata', img);
		handles.I_never_changed_imgMir = false;
		guidata(handles.figure1, handles)
	else
		if (handMir.image_type == 2)			% If we are processing a simple image
			h = mirone(get(handles.hImg,'cdata'));
			set(h,'Name','Threshold mask')
			return
		end
		% Ok, if we get here it's because the image has coordinates. Fish it all
		imgLims = getappdata(handMir.axes1,'ThisImageLims');
		tmp.X = imgLims(1:2);		tmp.Y = imgLims(3:4);		tmp.geog = handMir.geog;
		tmp.head = handMir.head;	tmp.head(5) = 0;			tmp.head(5) = 1;
		ProjWKT = getappdata(handMir.figure1,'ProjWKT');
		if (~isempty(ProjWKT))
			tmp.srsWKT = ProjWKT;
		else
			proj4 = getappdata(handMir.figure1,'Proj4');
			if (~isempty(proj4)),	tmp.srsWKT = proj4;		end
		end
		tmp.name = 'Threshold mask';
		mirone(get(handles.hImg,'cdata'), tmp)
	end

% ---------------------------------------------------------------------
function I = binarize_it(img, threshold, negative)
% Compute the binary image from the threshold level and eventual imcomplement it
	I = false(size(img));
	if (~negative)
		I(img < threshold) = 0;
		I(img > threshold) = 1;
	else
		I(img < threshold) = 1;
		I(img > threshold) = 0;
	end

% ---------------------------------------------------------------------
function level = triangle_th(handles, num_bins)
%     Triangle algorithm
%     This technique is due to Zack (Zack GW, Rogers WE, Latt SA (1977), 
%     "Automatic measurement of sister chromatid exchange frequency", 
%     J. Histochem. Cytochem. 25 (7): 741–53, )
%     A line is constructed between the maximum of the histogram at 
%     (b) and the lowest (or highest depending on context) value (a) in the 
%     histogram. The distance L normal to the line and between the line and 
%     the histogram h[b] is computed for all values from a to b. The level
%     where the distance between the histogram and the line is maximal is the 
%     threshold value (level). This technique is particularly effective 
%     when the object pixels produce a weak peak in the histogram.

%     Use Triangle approach to compute threshold (level) based on a 1D histogram. 
%
%     INPUTS
%         handles :   handles.img - The image
%         num_bins:   number of bins (e.g. gray levels)
%     OUTPUT
%         level   :   threshold value in the range [0 255];
% 
%     Dr B. Panneton, June, 2010
%     Agriculture and Agri-Food Canada
%     St-Jean-sur-Richelieu, Qc, Canad
%     bernard.panneton@agr.gc.ca
% http://www.mathworks.com/matlabcentral/fileexchange/28047

	if (nargin == 1),		num_bins = 256;		end
%   Find maximum of histogram and its location along the x axis
	lehisto = localImhist(handles.img);
    [h,xmax]= max(lehisto);
    xmax = round(mean(xmax));   %can have more than a single value!
    h=lehisto(xmax);
    
%   Find location of first and last non-zero values.
%   Values<h/10000 are considered zeros.
    indi = find(lehisto > h/10000);
    fnz=indi(1);
    lnz=indi(end);

%   Pick side as side with longer tail. Assume one tail is longer.
    lspan=xmax-fnz;
    rspan=lnz-xmax;
    if rspan>lspan  % then flip lehisto
        lehisto=fliplr(lehisto');
        a=num_bins-lnz+1;
        b=num_bins-xmax+1;
        isflip=1;
    else
        lehisto=lehisto';
        isflip=0;
        a=fnz;
        b=xmax;
    end
    
%   Compute parameters of the straight line from first non-zero to peak
%   To simplify, shift x axis by a (bin number axis)
    m=h/(b-a);
    
%   Compute distances
    x1=0:(b-a);
    y1=lehisto(x1+a);
    beta=y1+x1/m;
    x2=beta/(m+1/m);
    y2=m*x2;
    L=((y2-y1).^2+(x2-x1).^2).^0.5;

%   Obtain threshold as the location of maximum L.    
    level=find(max(L)==L);
    level=a+mean(level);
    
%   Flip back if necessary
    if isflip
        level = num_bins - level + 1;
    end
    
    level = (level / num_bins) * 255;
 
% ---------------------------------------------------------------------
function threshold = maxentropy(handles)
% maxentropie is a function for thresholding using Maximum Entropy
% 
% output = threshold ==> the threshold choosen by maxentropy
%  
% F.Gargouri
% Originaly from http://www.mathworks.com/matlabcentral/fileexchange/35158
%
% Lots of cleanings (no pre-allocations) but still room for vectorizations
% Joaquim Luis

	[n,m] = size(handles.img);
	h = localImhist(handles.img);
	%normalize the histogram ==>  hn(k)=h(k)/(n*m) ==> k  in [1 256]
	hn = h / (n*m);

	%Cumulative distribution function
	c = cumsum(hn);

	hl = zeros(1,256);
	hh = zeros(1,256);
	for t = 1:256
		%low range entropy
		cl = c(t);
		if cl > 0
			for (i = 1:t)
				if hn(i) > 0
					hl(t) = hl(t) - (hn(i)/cl)*log(hn(i)/cl);                      
				end
			end
		end

		%high range entropy
		ch = 1.0 - cl;  %constraint cl+ch=1
		if ch > 0
			for i = t+1:256
				if hn(i) > 0
					hh(t) = hh(t) - (hn(i)/ch)*log(hn(i)/ch);
				end
			end
		end
	end

	% choose best threshold

	h_max = hl(1) + hh(1);
	threshold = 0;
	entropie = hl + hh;
	entropie(1) = h_max;
	for t = 2:256
		if (entropie(t) > h_max)
			h_max = entropie(t);
			threshold = t-1;
		end
	end

% ---------------------------------------------------------------------
function threshold = minCE(handles)
% MinCE is a function for thresholding using Minimum Cross Entropy
% 
% output = threshold ==> the threshold choosen by MinCE
%  
% F.Gargouri
% Originaly from http://www.mathworks.com/matlabcentral/fileexchange/35157
%
% Lots of cleanings (no pre-allocations) but still room for vectorizations
% Joaquim Luis

	h = localImhist(handles.img);
	[n,m]=size(handles.img);
	%normalize the histogram ==>  hn(k)=h(k)/(n*m) ==> k  in [1 256]
	hn=h/(n*m);

	% entropy of gray level image
	imEntropy = 0;
	for i=1:256
		imEntropy = imEntropy + (i*hn(i)*log(i));
	end

	% MCE
	CE(1,256) = 0;
	for (t = 1:256)
		% moyenne de Low range image
		lowValue = 0;
		lowSum   = 0;
		for (i = 1:t)
			lowValue = lowValue + i*hn(i);
			lowSum = lowSum + hn(i);
		end
		if (lowSum > 0)
			lowValue = lowValue/lowSum;
		else
			lowValue = 1;
		end

		% average of High range image
		highValue = 0;
		highSum   = 0;
		for (i = t+1:256)
			highValue = highValue + i*hn(i);
			highSum   = highSum + hn(i);
		end
		if (highSum > 0)
			highValue = highValue/highSum;
		else
			highValue = 1;
		end

		% Entropy of low range 
		lowEntropy = 0;		logLowVal = log(lowValue);
		for (i = 1:t)
			lowEntropy = lowEntropy + i*hn(i)*logLowVal;
		end      

		% Entropy of high range 
		highEntropy = 0;	logHighVal = log(highValue);
		for (i = t+1:256)
			highEntropy = highEntropy + i*hn(i)*logHighVal;
		end

		% Cross Entropy 
		CE(t) = imEntropy - lowEntropy - highEntropy; 
	end

	% choose the best threshold

	D_min = CE(1);
	entropie = CE;
	entropie(1) = D_min;
	threshold = 0;
	for t = 2:256
		if (entropie(t) < D_min)
			D_min = entropie(t);
			threshold = t - 1;
		end
	end

% ---------------------------------------------------------------------
function level = isodata(handles)
%   ISODATA Compute global image threshold using iterative isodata method.
%   LEVEL = ISODATA(I) computes a global threshold (LEVEL) that can be
%   used to convert an intensity image to a binary image with IM2BW. LEVEL
%   is a normalized intensity value that lies in the range [0 255].
%   This iterative technique for choosing a threshold was developed by Ridler and Calvard .
%   The histogram is initially segmented into two parts using a starting threshold value such as 0 = 2B-1, 
%   half the maximum dynamic range. 
%   The sample mean (mf,0) of the gray values associated with the foreground pixels and the sample mean (mb,0) 
%   of the gray values associated with the background pixels are computed. A new threshold value 1 is now computed 
%   as the average of these two sample means. The process is repeated, based upon the new threshold, 
%   until the threshold value does not change any more. 
  
%
%   Example
%   -------
%       I = imread('blood1.tif');
%       level = graythresh(I);
%       BW = im2bw(I,level);
%       imshow(BW)
%
%
% Reference :T.W. Ridler, S. Calvard, Picture thresholding using an iterative selection method, 
%            IEEE Trans. System, Man and Cybernetics, SMC-8 (1978) 630-632.
%
% Originaly from http://www.mathworks.com/matlabcentral/fileexchange/3195
% Licence BSD

	% STEP 1: Compute mean intensity of image from histogram, set T=mean(I)
	[counts,N] = localImhist(handles.img);
	i    = 1;
	mu   = cumsum(counts);
	T(i) = (sum(N.*counts))/mu(end);
	T(i) = round(T(i));

	% STEP 2: compute Mean above T (MAT) and Mean below T (MBT) using T from step 1
	mu2 = cumsum(counts(1:T(i)));
	MBT = sum(N(1:T(i)).*counts(1:T(i)))/mu2(end);

	mu3 = cumsum(counts(T(i):end));
	MAT = sum(N(T(i):end).*counts(T(i):end))/mu3(end);
	i = i + 1;
	T(i) = round((MAT+MBT)/2);

	% STEP 3 to n: repeat step 2 if T(i)~=T(i-1)
	while abs(T(i)-T(i-1))>=1
		mu2 = cumsum(counts(1:T(i)));
		MBT = sum(N(1:T(i)).*counts(1:T(i)))/mu2(end);

		mu3 = cumsum(counts(T(i):end));
		MAT = sum(N(T(i):end).*counts(T(i):end))/mu3(end);

		i = i + 1;
		T(i) = round((MAT+MBT)/2); 
		Threshold = T(i);
	end

	% Normalize the threshold to the range [0 255].
	level = (Threshold - 1) / (N(end) - 1) * 255;


% --- Creates and returns a handle to the GUI figure. 
function thresholdit_LayoutFcn(h1)

set(h1, 'Position',[520 301 680 499],...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'MenuBar','none',...
'Name','Thresholdit',...
'NumberTitle','off',...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1');

axes('Parent',h1,...
'Units','pixels',...
'Position',[10 88 511 411],...
'CameraPosition',[0.5 0.5 9.16025403784439],...
'Tag','axes1');

axes('Parent',h1,...
'Units','pixels',...
'Position',[10 18 661 71],...
'CameraPosition',[0.5 0.5 9.16025403784439],...
'Tag','axes2');

uicontrol('Parent',h1, 'Position',[530 466 141 23],...
'Call',@thresholdit_uiCB,...
'String','Otsu''s method',...
'Tag','push_Otsu');

uicontrol('Parent',h1, 'Position',[530 436 141 23],...
'Call',@thresholdit_uiCB,...
'String','Maximum Entropy',...
'Tag','push_maxEntropy');

uicontrol('Parent',h1, 'Position',[530 406 141 23],...
'Call',@thresholdit_uiCB,...
'String','Minimum Cross-Entropy',...
'Tag','push_minCrossEntropy');

uicontrol('Parent',h1, 'Position',[530 376 141 23],...
'Call',@thresholdit_uiCB,...
'String','Isodata method',...
'Tag','push_isodata');

uicontrol('Parent',h1, 'Position',[530 310 73 19],...
'Call',@thresholdit_uiCB,...
'String','Single line',...
'Style','radiobutton',...
'TooltipString','Drag one line only to separate black and white',...
'Value',1,...
'Tag','radio_singleLine');

uicontrol('Parent',h1, 'Position',[530 346 141 23],...
'Call',@thresholdit_uiCB,...
'String','Triangle method',...
'Tag','push_triang');

uicontrol('Parent',h1, 'Position',[531 268 120 19],...
'Call',@thresholdit_uiCB,...
'String','Revert (negative)',...
'Style','checkbox',...
'TooltipString','Invert (negative) thresholded image.',...
'Tag','check_revert');

uicontrol('Parent',h1, 'Position',[530 228 141 23],...
'Call',@thresholdit_uiCB,...
'String','Clean white dust',...
'TooltipString','Remove the isolated small group of white pixels inside black zones.',...
'Tag','push_cleanDust');

uicontrol('Parent',h1, 'Position',[530 200 141 23],...
'Call',@thresholdit_uiCB,...
'String','Fill holes',...
'TooltipString','Fill holes.',...
'Tag','push_fillHoles');

uicontrol('Parent',h1, 'Position',[608 310 67 19],...
'Call',@thresholdit_uiCB,...
'String','Window',...
'Style','radiobutton',...
'TooltipString','Iside window is 1, outside is 0. You can drag the widow and/or resize it.',...
'Tag','radio_window');

uicontrol('Parent',h1, 'Position',[530 142 120 19],...
'String','Apply to original',...
'Style','checkbox',...
'TooltipString','Mask the original image when "Goog, I like it"',...
'Tag','check_applyOrig');

uicontrol('Parent',h1, 'Position',[530 118 141 23],...
'Call',@thresholdit_uiCB,...
'FontSize',10,...
'FontWeight','bold',...
'String',' Good, I like it',...
'TooltipString','Create final threshold image',...
'Tag','push_OK');

function thresholdit_uiCB(hObject, eventdata)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));
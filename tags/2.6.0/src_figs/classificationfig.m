function varargout = classificationfig(varargin)
%   Perform image classification, supervised an unsupervised, using k-means

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
%
% The dcKMeans function is Copyright (C) David Corney 2000		D.Corney@cs.ucl.ac.uk

% $Id: $

	hObject = figure('Vis','off');
	classificationfig_LayoutFcn(hObject);
	handles = guihandles(hObject);

	if (~isempty(varargin))
		handles.hCallingFig = varargin{1};
		handMir = guidata(handles.hCallingFig);
		if (handMir.no_file)
			errordlg('You didn''t even load a file. What are you expecting then?','ERROR')
			delete(hObject);    return
		end
	else
		delete(hObject)
		return
	end

	% Try to position this figure glued to the right of calling figure
	posThis = get(hObject,'Pos');
	posParent = get(handles.hCallingFig,'Pos');
	ecran = get(0,'ScreenSize');
	xLL = posParent(1) + posParent(3) + 6;
	xLR = xLL + posThis(3);
	if (xLR > ecran(3))         % If figure is partially out, bring totally into screen
		xLL = ecran(3) - posThis(3);
	end
	yLL = (posParent(2) + posParent(4)/2) - posThis(4) / 2;
	set(hObject,'Pos',[xLL yLL posThis(3:4)])

	handlesMir = guidata(handles.hCallingFig);
	handles.head = handlesMir.head;
	handles.image_type = handlesMir.image_type;
	handles.hCallingAxes = handlesMir.axes1;
	handles.hImg = handlesMir.hImg;
	handles.image_type = handlesMir.image_type;
	handles.imgSize = size(get(handles.hImg,'CData'));
	handles.supervised = 1;
	handles.nClasses = 3;
	handles.selectedClasses = 0;
	handles.nNeighbors = 3;
	handles.bg_color = uint8(round(handlesMir.bg_color * 255));

	% Add this figure handle to the carra?as list
	plugedWin = getappdata(handles.hCallingFig,'dependentFigs');
	plugedWin = [plugedWin hObject];
	setappdata(handles.hCallingFig,'dependentFigs',plugedWin);

	%------------ Give a Pro look (3D) to the frame boxes  -------------------------------
	new_frame3D(hObject, handles.text_isolate)
	%------------- END Pro look (3D) -----------------------------------------------------

	str = sprintf(['With the supervised classification\n'...
		'seed point color is computed by\n'...
		'averaging pixels inside a square\n'...
		'window with this number of points side']);
	set(handles.edit_nNeighbors,'Tooltip',str)
	str = sprintf(['Background color used in color isolation.\n'...
		'Push to select a different color.']);
	set(handles.push_bgColor,'BackgroundColor', handlesMir.bg_color, 'ForegroundColor', 1-handlesMir.bg_color, 'Tooltip',str)
    
	guidata(hObject, handles);
	set(hObject,'Visible','on');
	if (nargout),   varargout{1} = handles.output;  end

% ----------------------------------------------------------------------------
function radio_supervised_CB(hObject, handles)
	if (get(hObject,'Val'))
		set(handles.radio_unsupervised,'Val',0)
		set(handles.edit_nClasses,'Enable','off')
		handles.supervised = 1;
		guidata(handles.figure1,handles)
	else
		set(hObject,'Val',1)
	end

% ----------------------------------------------------------------------------
function radio_unsupervised_CB(hObject, handles)
    if (get(hObject,'Val'))
        set(handles.radio_supervised,'Val',0)
        set(handles.edit_nClasses,'Enable','on')
        handles.supervised = 0;
        guidata(handles.figure1,handles)
    else
        set(hObject,'Val',1)
    end

% ----------------------------------------------------------------------------
function edit_nClasses_CB(hObject, handles)
    xx = str2double(get(hObject,'String'));
    if (isnan(xx)),     set(hObject,'String',3);    return;     end
    handles.nClasses = round(xx);
    if (handles.nClasses < 2)       % Minimum allowed is 2
        set(hObject,'String',2)
        handles.nClasses = 2;
    end
    guidata(handles.figure1,handles)

% ----------------------------------------------------------------------------
function toggle_clickDefine_CB(hObject, handles)
    figure(handles.hCallingFig)         % Bring the figure containing image forward

    h = findobj(handles.hCallingAxes,'Type','line','Tag','ClassifyPoly');
    if (~isempty(h)),   delete(h);      end     % If it exists, kill it. There can be only one.
    
    [xp,yp] = getline_j(handles.hCallingFig);
    set(hObject,'Val',0)
    nClusters = length(xp);
    if (nClusters < 2),     return;     end
    h = line('XData', xp, 'YData', yp,'Parent',handles.hCallingAxes,'Color','k','LineStyle','none',...
        'Marker','o','MarkerFaceColor','y','MarkerSize',5,'Tag','ClassifyPoly');
  
    cmenuHand = uicontextmenu('Parent',handles.hCallingFig);
    set(h, 'UIContextMenu', cmenuHand);
    uimenu(cmenuHand, 'Label', 'Delete', 'Call', 'delete(gco)');
    ui_edit_polygon(h)

% -------------------------------------------------------------------------------------
function push_compute_CB(hObject, handles)
%
	if (handles.supervised)
		h = findobj(handles.hCallingAxes,'Type','line','Tag','ClassifyPoly');
		if (isempty(h)),		return,		end
		xp = get(h,'XData');	yp = get(h,'YData');
		nClusters = length(xp);
		[c,r] = getpixcoords(handles,xp,yp);
		r = round(r);    c = round(c);
	else
		nClusters = handles.nClasses;
	end

	if (~ishandle(handles.hImg))
		delete(handles.figure1)
		warndlg('Apparently you loaded another file. So you need to call the "k-means classification" option again','Warning')
		return
	end

	img = double(get(handles.hImg,'CData')) / 255;      % Get the image
	nx = size(img,2);    ny = size(img,1);
	W = fix(handles.nNeighbors / 2);                    % If W = 3 it means a 3x3 window centered on point i

	if (ndims(img) == 2)
		segcolors = img(:);
		to_resize = size(img);
		if (handles.supervised)
			colors = zeros(nClusters,1);
			for (i=1:nClusters)
				colors(i) = mean(mean( img( max(r(i)-W,1):min(r(i)+W,ny), max(c(i)-W,1):min(c(i)+W,nx) ) ));
				%colors(i,:) = img(r(i),c(i));
			end
		end
	else
		%Get individual RGB colors included in segmentation
		rc = img(:,:,1);        gc = img(:,:,2);        bc = img(:,:,3);
		to_resize = size(rc);
		segcolors = [rc(:) gc(:) bc(:)];
		if (handles.supervised)
			colors = zeros(nClusters,3);
			for (i=1:nClusters)
				tmp1 = mean(mean( img( max(r(i)-W,1):min(r(i)+W,ny), max(c(i)-W,1):min(c(i)+W,nx),1 ) ));
				tmp2 = mean(mean( img( max(r(i)-W,1):min(r(i)+W,ny), max(c(i)-W,1):min(c(i)+W,nx),2 ) ));
				tmp3 = mean(mean( img( max(r(i)-W,1):min(r(i)+W,ny), max(c(i)-W,1):min(c(i)+W,nx),3 ) ));
				%colors(i,:) = img(r(i),c(i),:);
				colors(i,:) = [tmp1 tmp2 tmp3];
			end
		end
	end

	if (handles.supervised)
		%Idx = kmeans_j(segcolors,nClusters,'start',colors,'emptyaction','singleton','display','win');
		Idx = dcKMeans(segcolors, nClusters, colors, 50);
	else
		%[Idx,colors] = kmeans_j(segcolors,nClusters,'emptyaction','singleton','display','win');
		[Idx,colors] = dcKMeans(segcolors, nClusters, [], 50);
	end

	if (isempty(Idx)),		return,		end		% Aborted

	if (size(colors,2) == 1)    % Unsupervised on a BW image
		colors = [colors(:) colors(:) colors(:)];
	end

	Idx = uint8(Idx-1);
	Idx = reshape(Idx,to_resize);

	if (handles.image_type == 2 || handles.image_type == 20)
		h = mirone(Idx);
		set(h,'ColorMap',colors,'Name','Image classification')
	else
		tmp.X = handles.head(1:2);  tmp.Y = handles.head(3:4);  tmp.head = handles.head;
		tmp.cmap = colors;
		tmp.name = 'Image classification';
		h = mirone(Idx,tmp);
	end

	handMir = guidata(h);
	set(handMir.hImg, 'CDataMapping', 'scaled')

	set(handles.listbox_classes,'String',0:nClusters-1, 'Max', nClusters-1)		% To allow multiple class selection
	set(handles.push_getClass,'Enable','on')
	handles.hMirClass.hImg = handMir.hImg;
	guidata(handles.figure1, handles)
    
% -------------------------------------------------------------------------
function radio_asColor_CB(hObject, handles)
	if (~get(hObject,'Value')),		set(hObject,'Value',1),		return,		end
	set(handles.radio_asMask,'Val',0)

% -------------------------------------------------------------------------
function radio_asMask_CB(hObject, handles)
	if (~get(hObject,'Value')),		set(hObject,'Value',1),		return,		end
	set(handles.radio_asColor,'Val',0)

% -------------------------------------------------------------------------
function listbox_classes_CB(hObject, handles)
	handles.selectedClasses = get(hObject,'Value') - 1;		% -1 because first class is 0
	guidata(handles.figure1, handles)

% -------------------------------------------------------------------------
function push_getClass_CB(hObject, handles)
	if (~ishandle(handles.hImg))
		errordlg('You said bye bye to the original image and now ask to get someting out of it! Funny isn''t it?','Error'), return
	end
	img_orig  = get(handles.hImg, 'CData');
	img_class = get(handles.hMirClass.hImg,'CData');
	threeD = 3;
	if (ndims(img_orig) == 2),		threeD = 1;	end
	if (get(handles.radio_asColor, 'Val'))			% Color extraction
		img = img_orig;
		for (k = 1:numel(handles.selectedClasses))
			mask = (img_class == handles.selectedClasses(k));
			if (k == 1)
				for (n = 1:threeD)
					tmp = img(:,:,n);		tmp(~mask) = handles.bg_color(n);
					img(:,:,n) = tmp;
				end
				clear tmp
			else				% Only comes here if more than one class was selected
				img(repmat(mask, [1 1 threeD])) = img_orig(repmat(mask, [1 1 threeD]));
			end
		end
		clear mask
	else											% Mask extraction
		img = false(size(img_class));
		for (k = 1:numel(handles.selectedClasses))
			mask = (img_class == handles.selectedClasses(k));
			img  = ( img | mask );
		end		
	end

	if (handles.image_type == 2 || handles.image_type == 20)
		h = mirone(img);
		set(h,'Name','Class Mask')
	else
		tmp.X = handles.head(1:2);  tmp.Y = handles.head(3:4);  tmp.head = handles.head;
		if (islogical(img)),        tmp.name = 'Class Mask';
		else						tmp.name = 'Class Image';
		end
		mirone(img,tmp)
	end

% -------------------------------------------------------------------------------------
function push_bgColor_CB(hObject, handles)
    c = uisetcolor;
    if (isempty(c)),	return,		end
	handles.bg_color = uint8(round(c*255));
	set(hObject,'BackgroundColor', c, 'ForegroundColor', 1-c)
    guidata(handles.figure1,handles)

% -------------------------------------------------------------------------------------
function edit_nNeighbors_CB(hObject, handles)
    xx = str2double(get(hObject,'String'));
    if (isnan(xx)),		set(hObject,'String',3),	return,		end
    handles.nNeighbors = round(abs(xx));
    guidata(handles.figure1,handles)

% -------------------------------------------------------------------------------------
function [x,y] = getpixcoords(handles,x,y)
    % Convert x,y to pixel coordinates (they are not when the image has other coordinates)
    if (handles.head(7))                % Image is pixel registered
        X = [handles.head(1) handles.head(2)] + [handles.head(8) -handles.head(8)]/2;
        Y = [handles.head(3) handles.head(4)] + [handles.head(9) -handles.head(9)]/2;
    else                                % Image is grid registered
        X = [handles.head(1) handles.head(2)];
        Y = [handles.head(3) handles.head(4)];
    end
    x = getPixel_coords(handles.imgSize(2),X,x);
    y = getPixel_coords(handles.imgSize(1),Y,y);

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

% -------------------------------------------------------------------------------------
function [Classes,Centres,FinalDistance] = dcKMeans(Data,k,InitCentres,MaxIters)
%DCKMEANS Performs k-means clustering
% Classes=DCKMEANS(Data,k,InitCentres) where Data is a matrix with
% variables in columns and records in rows, k is the number of 
% clusters to search for and InitCentres is a list of initial centres.
% Classes=DCKMEANS(Data,k) randomly chooses initial centres from data set.
% [Classes,Centres,FinalDistance]=DCKMEANS(Data,k,InitCentres) also returns
% the centres and distances of each point to its nearest centre.

% Copyright (C) David Corney 2000		D.Corney@cs.ucl.ac.uk

	if (nargin < 3)
		InitCentres = ChooseInitialCentres(Data,k);		%randomly choose starting point (where needed)
	end
	if (isempty(InitCentres)),  InitCentres = ChooseInitialCentres(Data,k);     end
	Centres = InitCentres;
	OldCentres = Centres;
	if (nargin < 4),   MaxIters = 100;  end
	boing = 0;

	R = size(Data,1);

	DataSq = repmat(sum(Data.^2,2),1,k);	%sum squared data - save re-calculating repeatedly later
	%Do we need DataSq? It's constant, and we're minimsing things...

	hWait = aguentabar(0,'title','Please wait ...','CreateCancelBtn');
	for i = 1:MaxIters
		Dist = DataSq + repmat(sum((Centres.^2)',1),R,1) - 2.*(Data*(Centres'));   %i.e. d^2 = (x-c)^2 = x^2 + c^2 -2xc
		[D,Centre] = min(Dist,[],2);		%label of nearest centre for each point

		for (j = 1:k)
			idx = find(Centre == j);
			if (~isempty(idx))
				Centres(j,:) = mean(Data(idx,:));
			end
		end

		Change = sum(sum(abs(OldCentres-Centres)));
		if (Change < 1e-8),		break,		end		%Have we converged yet?
		OldCentres = Centres;
		hWait = aguentabar(i/MaxIters);
		if (isnan(hWait))		% The Cancel button was hit and the aguentabar destroyed
			boing = 1;
			break
		end
	end
	clear D Centre DataSq
	if (ishandle(hWait)),	delete(hWait),	end		% When convergence finish before MaxIters

	if (~boing)					% Normal termination
		[FinalDistance,Classes] = min(Dist,[],2);		%label points one last time
	else						% Aborted
		Classes = [];   Centres = [];   FinalDistance = [];
	end

% -------------------------------------------------------------------------
function [C,perm] = ChooseInitialCentres(set,k)
%CHOOSEINITIALCENTRES Randomly picks sample points
% Centres=CHOOSEINITIALCENTRES(Data,k) where Data is
% the data matrix, and k is the number of points required.
% Used to initialise various clustering routines.
%
%(C) David Corney (2000)   		D.Corney@cs.ucl.ac.uk

	perm = randperm(size(set,1));
	perm = perm(1:k);
	C = set(perm,:);	

% --- Creates and returns a handle to the GUI figure. 
function classificationfig_LayoutFcn(h1)

set(h1,...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'MenuBar','none',...
'Name','kmeans classification',...
'NumberTitle','off',...
'Pos',[520 710 210 210],...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@classification_uiCB,...
'Enable','off',...
'Pos',[164 158 31 21],...
'String','3',...
'Style','edit',...
'Tooltip','Number of classes to identify',...
'Tag','edit_nClasses');

uicontrol('Parent',h1,...
'Call',@classification_uiCB,...
'FontName','Helvetica',...
'Pos',[11 159 76 21],...
'String','click-define',...
'Style','togglebutton',...
'Tooltip','Left click on the target colors/gray levels. Stop with a right click',...
'Tag','toggle_clickDefine');

uicontrol('Parent',h1, 'Pos',[100 127 3 85], 'Style','frame', 'Tag','frame1');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'Pos',[111 162 51 15],...
'String','N classes',...
'Style','text');

uicontrol('Parent',h1,...
'Call',@classification_uiCB,...
'FontName','Helvetica',...
'FontSize',10,...
'Pos',[113 126 80 21],...
'String','Compute',...
'Tag','push_compute');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@classification_uiCB,...
'Pos',[62 128 31 21],...
'String','3',...
'Style','edit',...
'Tag','edit_nNeighbors');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'Pos',[9 132 52 15],...
'String','Neigbours',...
'Style','text');

uicontrol('Parent',h1,'Pos',[102 154 109 3],'Style','frame','Tag','frame2');

uicontrol('Parent',h1,...
'Call',@classification_uiCB,...
'FontName','Helvetica',...
'FontSize',9,...
'Pos',[10 186 85 16],...
'String','Supervised',...
'Style','radiobutton',...
'Tooltip','Do supervised image classification',...
'Value',1,...
'Tag','radio_supervised');

uicontrol('Parent',h1,...
'Call',@classification_uiCB,...
'FontName','Helvetica',...
'FontSize',9,...
'Pos',[109 186 98 16],...
'String','Unsupervised',...
'Style','radiobutton',...
'Tooltip','Do unsupervised image classification',...
'Tag','radio_unsupervised');

uicontrol('Parent',h1,'Pos',[0 107 211 3],'Style','frame','Tag','frame_isolate');

uicontrol('Parent',h1,...
'Call',@classification_uiCB,...
'FontName','Helvetica',...
'FontSize',9,...
'Pos',[66 2 140 21],...
'String','Isolate selected class',...
'Tooltip','Do it.',...
'Enable','off',...
'Tag','push_getClass');

uicontrol('Parent',h1,...
'Call',@classification_uiCB,...
'FontName','Helvetica',...
'FontSize',8,...
'Pos',[156 69 50 21],...
'String','BG color',...
'Tag','push_bgColor');

uicontrol('Parent',h1,...
'Call',@classification_uiCB,...
'Pos',[80 47 110 15],...
'String','Isolate as color',...
'Style','radiobutton',...
'Tooltip','Extract the original colors corresponding to the selected class(es)',...
'Value',1,...
'Tag','radio_asColor');

uicontrol('Parent',h1,...
'Call',@classification_uiCB,...
'Pos',[80 28 110 15],...
'String','Isolate as mask',...
'Style','radiobutton',...
'Tooltip','Extract the selected class(es) as binary image(s)',...
'Tag','radio_asMask');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'HorizontalAlignment','left',...
'Pos',[55 74 51 15],...
'String','Classes',...
'Style','text');

uicontrol('Parent',h1,'FontName','Helvetica',...
'FontSize',8,...
'Pos',[19 100 150 16],...
'String','Isolate classes (if wished)',...
'Style','text',...
'Tag','text_isolate');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@classification_uiCB,...
'Pos',[4 3 50 85],...
'Style','listbox',...
'Tooltip','Select one or more classes to islolate.',...
'Value',1,...
'Tag','listbox_classes');

function classification_uiCB(hObject, eventdata)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));

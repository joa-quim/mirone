function varargout = classificationfig(varargin)
%   Perform image classification, supervised an unsupervised, using kmeans
%
% M-File changed by desGUIDE 
 
	hObject = figure('Tag','figure1','Visible','off');
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
	bgcolor = get(0,'DefaultUicontrolBackgroundColor');
	framecolor = max(min(0.65*bgcolor,[1 1 1]),[0 0 0]);
	h_f = [handles.frame1 handles.frame2 handles.frame_isolate];
	for i=1:numel(h_f)
		frame_size = get(h_f(i),'Position');
		f_bgc = get(h_f(i),'BackgroundColor');
		usr_d = get(h_f(i),'UserData');
		if abs(f_bgc(1)-bgcolor(1)) > 0.01           % When the frame's background color is not the default's
			frame3D(hObject,frame_size,framecolor,f_bgc,usr_d)
		else
			frame3D(hObject,frame_size,framecolor,'',usr_d)
			delete(h_f(i))
		end
	end
	delete(handles.text_isolate)		% Recreate this text. Better than call the stupid uistack function
	handles.text_isolate = uicontrol('Parent',hObject,'FontName','Helvetica','FontSize',9,'Position',[19 100 160 16],...
		'String','Isolate classes (if wished)','Style','text');
	%------------- END Pro look (3D) -------------------------------------------------------

	str = sprintf(['With the supervised classification\n'...
		'seed point color is computed by\n'...
		'averaging pixels inside a square\n'...
		'window with this number of points side']);
	set(handles.edit_nNeighbors,'TooltipString',str)
	str = sprintf(['Background color used in color isolation.\n'...
		'Push to select a different color.']);
	set(handles.push_bgColor,'BackgroundColor', handlesMir.bg_color, 'ForegroundColor', 1-handlesMir.bg_color, 'TooltipString',str)
    
	guidata(hObject, handles);
	set(hObject,'Visible','on');
	if (nargout),   varargout{1} = handles.output;  end

% ----------------------------------------------------------------------------
function radio_supervised_Callback(hObject, eventdata, handles)
	if (get(hObject,'Val'))
		set(handles.radio_unsupervised,'Val',0)
		set(handles.edit_nClasses,'Enable','off')
		handles.supervised = 1;
		guidata(handles.figure1,handles)
	else
		set(hObject,'Val',1)
	end

% ----------------------------------------------------------------------------
function radio_unsupervised_Callback(hObject, eventdata, handles)
    if (get(hObject,'Val'))
        set(handles.radio_supervised,'Val',0)
        set(handles.edit_nClasses,'Enable','on')
        handles.supervised = 0;
        guidata(handles.figure1,handles)
    else
        set(hObject,'Val',1)
    end

% ----------------------------------------------------------------------------
function edit_nClasses_Callback(hObject, eventdata, handles)
    xx = str2double(get(hObject,'String'));
    if (isnan(xx)),     set(hObject,'String',3);    return;     end
    handles.nClasses = round(xx);
    if (handles.nClasses < 2)       % Minimum allowed is 2
        set(hObject,'String',2)
        handles.nClasses = 2;
    end
    guidata(handles.figure1,handles)

% ----------------------------------------------------------------------------
function toggle_clickDefine_Callback(hObject, eventdata, handles)
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
    uimenu(cmenuHand, 'Label', 'Delete', 'Callback', 'delete(gco)');
    ui_edit_polygon(h)

% -------------------------------------------------------------------------------------
function push_compute_Callback(hObject, eventdata, handles)
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
function radio_asColor_Callback(hObject, eventdata, handles)
	if (~get(hObject,'Value')),		set(hObject,'Value',1),		return,		end
	set(handles.radio_asMask,'Val',0)

% -------------------------------------------------------------------------
function radio_asMask_Callback(hObject, eventdata, handles)
	if (~get(hObject,'Value')),		set(hObject,'Value',1),		return,		end
	set(handles.radio_asColor,'Val',0)

% -------------------------------------------------------------------------
function listbox_classes_Callback(hObject, eventdata, handles)
	handles.selectedClasses = get(hObject,'Value') - 1;		% -1 because first class is 0
	guidata(handles.figure1, handles)

% -------------------------------------------------------------------------
function push_getClass_Callback(hObject, eventdata, handles)
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
        h = mirone(img,tmp);
    end

% -------------------------------------------------------------------------------------
function push_bgColor_Callback(hObject, eventdata, handles)
    c = uisetcolor;
    if (isempty(c)),	return,		end
	handles.bg_color = uint8(round(c*255));
	set(hObject,'BackgroundColor', c, 'ForegroundColor', 1-c)
    guidata(handles.figure1,handles)

% -------------------------------------------------------------------------------------
function edit_nNeighbors_Callback(hObject, eventdata, handles)
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
    x = localAxes2pix(handles.imgSize(2),X,x);
    y = localAxes2pix(handles.imgSize(1),Y,y);

% -------------------------------------------------------------------------------------
function pixelx = localAxes2pix(dim, x, axesx)
%   Convert axes coordinates to pixel coordinates.
%   PIXELX = AXES2PIX(DIM, X, AXESX) converts axes coordinates
%   (as returned by get(gca, 'CurrentPoint'), for example) into
%   pixel coordinates.  X should be the vector returned by
%   X = get(image_handle, 'XData') (or 'YData').  DIM is the
%   number of image columns for the x coordinate, or the number
%   of image rows for the y coordinate.

	xfirst = x(1);      xlast = x(max(size(x)));	
	if (dim == 1)
        pixelx = axesx - xfirst + 1;        return;
	end
	xslope = (dim - 1) / (xlast - xfirst);
	if ((xslope == 1) & (xfirst == 1))
        pixelx = axesx;
	else
        pixelx = xslope * (axesx - xfirst) + 1;
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

[R,C]=size(Data);

DataSq = repmat(sum(Data.^2,2),1,k);	%sum squared data - save re-calculating repeatedly later
%Do we need DataSq? It's constant, and we're minimsing things...

hWait = waitbar(0,'Please wait ...','CreateCancelBtn','delete(gcf)');
for i = 1:MaxIters
	Dist = DataSq + repmat(sum((Centres.^2)',1),R,1) - 2.*(Data*(Centres'));   %i.e. d^2 = (x-c)^2 = x^2 + c^2 -2xc
	[D,Centre] = min(Dist,[],2);		%label of nearest centre for each point
	
	for j=1:k
		idx = find(Centre == j);
		if (~isempty(idx))
			Centres(j,:) = mean(Data(idx,:));
		end
	end
	
	Change = sum(sum(abs(OldCentres-Centres)));
	if (Change < 1e-8),		break,		end		%Have we converged yet?
	OldCentres=Centres;
	waitbar(i/MaxIters)
	if (~ishandle(hWait)),        % The Cancel button was hit and the waitbar destroyed
		boing = 1;  break
	end
end
clear D Centre DataSq

if (~boing)             % Normal termination
    delete(hWait)
    [FinalDistance,Classes] = min(Dist,[],2);		%label points one last time
else                    % Aborted
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
'PaperUnits',get(0,'defaultfigurePaperUnits'),...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'MenuBar','none',...
'Name','kmeans classification',...
'NumberTitle','off',...
'Position',[520 710 210 210],...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@classificationfig_uicallback,h1,'edit_nClasses_Callback'},...
'Enable','off',...
'Position',[164 158 31 21],...
'String','3',...
'Style','edit',...
'TooltipString','Number of classes to identify',...
'Tag','edit_nClasses');

uicontrol('Parent',h1,...
'Callback',{@classificationfig_uicallback,h1,'toggle_clickDefine_Callback'},...
'FontName','Helvetica',...
'Position',[11 159 76 21],...
'String','click-define',...
'Style','togglebutton',...
'TooltipString','Left click on the target colors/gray levels. Stop with a right click',...
'Tag','toggle_clickDefine');

uicontrol('Parent',h1, 'Position',[100 127 3 85], 'Style','frame', 'Tag','frame1');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'Position',[111 162 51 15],...
'String','N classes',...
'Style','text');

uicontrol('Parent',h1,...
'Callback',{@classificationfig_uicallback,h1,'push_compute_Callback'},...
'FontName','Helvetica',...
'FontSize',10,...
'Position',[113 126 80 21],...
'String','Compute',...
'Tag','push_compute');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@classificationfig_uicallback,h1,'edit_nNeighbors_Callback'},...
'Position',[62 128 31 21],...
'String','3',...
'Style','edit',...
'Tag','edit_nNeighbors');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'Position',[9 132 52 15],...
'String','Neigbours',...
'Style','text');

uicontrol('Parent',h1,'Position',[102 154 109 3],'Style','frame','Tag','frame2');

uicontrol('Parent',h1,...
'Callback',{@classificationfig_uicallback,h1,'radio_supervised_Callback'},...
'FontName','Helvetica',...
'FontSize',9,...
'Position',[10 186 85 16],...
'String','Supervised',...
'Style','radiobutton',...
'TooltipString','Do supervised image classification',...
'Value',1,...
'Tag','radio_supervised');

uicontrol('Parent',h1,...
'Callback',{@classificationfig_uicallback,h1,'radio_unsupervised_Callback'},...
'FontName','Helvetica',...
'FontSize',9,...
'Position',[109 186 98 16],...
'String','Unsupervised',...
'Style','radiobutton',...
'TooltipString','Do unsupervised image classification',...
'Tag','radio_unsupervised');

uicontrol('Parent',h1,'Position',[0 107 211 3],'Style','frame','Tag','frame_isolate');

uicontrol('Parent',h1,...
'Callback',{@classificationfig_uicallback,h1,'push_getClass_Callback'},...
'FontName','Helvetica',...
'FontSize',9,...
'Position',[66 2 140 21],...
'String','Isolate selected class',...
'TooltipString','Do it.',...
'Enable','off',...
'Tag','push_getClass');

uicontrol('Parent',h1,...
'Callback',{@classificationfig_uicallback,h1,'push_bgColor_Callback'},...
'FontName','Helvetica',...
'FontSize',8,...
'Position',[156 69 50 21],...
'String','BG color',...
'Tag','push_bgColor');

uicontrol('Parent',h1,...
'Callback',{@classificationfig_uicallback,h1,'radio_asColor_Callback'},...
'Position',[80 47 110 15],...
'String','Isolate as color',...
'Style','radiobutton',...
'TooltipString','Extract the original colors corresponding to the selected class(es)',...
'Value',1,...
'Tag','radio_asColor');

uicontrol('Parent',h1,...
'Callback',{@classificationfig_uicallback,h1,'radio_asMask_Callback'},...
'Position',[80 28 110 15],...
'String','Isolate as mask',...
'Style','radiobutton',...
'TooltipString','Extract the selected class(es) as binary image(s)',...
'Tag','radio_asMask');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'HorizontalAlignment','left',...
'Position',[55 74 51 15],...
'String','Classes',...
'Style','text');

uicontrol('Parent',h1,'FontName','Helvetica',...
'FontSize',9,...
'Position',[19 100 140 16],...
'String','Isolate classes (if wished)',...
'Style','text',...
'Tag','text_isolate');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@classificationfig_uicallback,h1,'listbox_classes_Callback'},...
'Position',[4 3 50 85],...
'Style','listbox',...
'TooltipString','Select one or more classes to islolate.',...
'Value',1,...
'Tag','listbox_classes');

function classificationfig_uicallback(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
feval(callback_name,hObject,[],guidata(h1));

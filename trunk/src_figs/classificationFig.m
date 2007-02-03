function varargout = classificationFig(varargin)
%   Perform image classification, supervised an unsupervised, using kmeans
%
% M-File changed by desGUIDE 
% hObject    handle to figure
% varargin   command line arguments to classificationFig (see VARARGIN)
 
hObject = figure('Tag','figure1','Visible','off');
handles = guihandles(hObject);
guidata(hObject, handles);
classificationFig_LayoutFcn(hObject,handles);
handles = guihandles(hObject);

if (~isempty(varargin))
    handles.hCallingFig = varargin{1};
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
handles.nNeighbors = 3;

% Add this figure handle to the carraças list
plugedWin = getappdata(handles.hCallingFig,'dependentFigs');
plugedWin = [plugedWin hObject];
setappdata(handles.hCallingFig,'dependentFigs',plugedWin);

str = sprintf(['With the supervised classification\n'...
    'seed point color is computed by\n'...
    'averaging pixels inside a square\n'...
    'window with this number of points side']);
set(handles.edit_nNeighbors,'TooltipString',str)
    
% Choose default command line output for classificationFig_export
handles.output = hObject;
guidata(hObject, handles);
varargout{1} = handles.output;

set(hObject,'Visible','on');

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
        if (isempty(h)),    return;     end
        xp = get(h,'XData');
        yp = get(h,'YData');
        nClusters = length(xp);
        [c,r] = getpixcoords(handles,xp,yp);
        r = round(r);    c = round(c);
    end

    img = double(get(handles.hImg,'CData')) / 255;      % Get the image
    nx = size(img,2);    ny = size(img,1);
    W = fix(handles.nNeighbors / 2);                    % If W = 3 it means a 3x3 window centered on point i
    
    if (ndims(img) == 2)
        segcolors = img(:);
        to_resize = size(img);
        colors = zeros(nClusters,1);
        if (handles.supervised)
            for (i=1:nClusters)
                colors(i) = mean(mean( img( max(r(i)-W,1):min(r(i)+W,ny), max(c(i)-W,1):min(c(i)+W,nx) ) ));
                %colors(i,:) = img(r(i),c(i));
            end
        end
    else
        %Get individual RGB colors included in segmentation
        rc = img(:,:,1);
        gc = img(:,:,2);
        bc = img(:,:,3);
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
        Idx = dcKMeans(segcolors,nClusters,colors,50);
    else
        %[Idx,colors] = kmeans_j(segcolors,handles.nClasses,'emptyaction','singleton','display','win');
        [Idx,colors] = dcKMeans(segcolors,handles.nClasses,[],50);
    end
    
    if (isempty(Idx)),  return;     end     % Aborted
    
    Idx = uint8(Idx-1);
    Idx = reshape(Idx,to_resize);
    
    if (handles.image_type == 2 || handles.image_type == 20)
        h = mirone(Idx);
        set(h,'ColorMap',colors,'Name','Image classification')
    else
        tmp.X = handles.head(1:2);  tmp.Y = handles.head(3:4);  tmp.head = handles.head;
        tmp.cmap = colors;
        tmp.name = 'Image classification';
        mirone(Idx,tmp);
    end
    
%     figure; imshow(Idx,[]);
%     colormap(colors);
%     tmp = colorbar;
%     set(tmp,'ytick',1:nClusters);

% -------------------------------------------------------------------------------------
function edit_nNeighbors_Callback(hObject, eventdata, handles)
    xx = str2double(get(hObject,'String'));
    if (isnan(xx)),     set(hObject,'String',3);    return;     end
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
   if (Change < 1e-8)	%Have we converged yet?
      break
   end
   OldCentres=Centres;
   waitbar(i/MaxIters)
   if (~ishandle(hWait))        % The Cancel button was hit and the waitbar destroyed
       boing = 1;  break
   end
end

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
function classificationFig_LayoutFcn(h1,handles);

set(h1,...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'MenuBar','none',...
'Name','kmeans classification',...
'NumberTitle','off',...
'PaperPosition',get(0,'defaultfigurePaperPosition'),...
'Position',[520 710 210 90],...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@classificationFig_uicallback,h1,'edit_nClasses_Callback'},...
'Enable','off',...
'Position',[163 39 31 21],...
'String','3',...
'Style','edit',...
'TooltipString','Number of classes to identify',...
'Tag','edit_nClasses');

uicontrol('Parent',h1,...
'Callback',{@classificationFig_uicallback,h1,'toggle_clickDefine_Callback'},...
'Position',[10 40 76 23],...
'String','click-define',...
'Style','togglebutton',...
'TooltipString','Left click on the target colors/gray levels. Stop with a right click',...
'Tag','toggle_clickDefine');

uicontrol('Parent',h1,'Position',[100 35 2 55],'Style','frame');

uicontrol('Parent',h1,'Position',[110 43 51 15],...
'String','N classes','Style','text');

uicontrol('Parent',h1,...
'Callback',{@classificationFig_uicallback,h1,'push_compute_Callback'},...
'FontName','Helvetica',...
'FontSize',10,...
'Position',[112 7 80 23],...
'String','Compute',...
'Tag','push_compute');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@classificationFig_uicallback,h1,'edit_nNeighbors_Callback'},...
'Position',[61 9 31 21],...
'String','3',...
'Style','edit',...
'Tag','edit_nNeighbors');

uicontrol('Parent',h1,'Position',[8 13 52 15],...
'String','Neigbours','Style','text');

uicontrol('Parent',h1,'Position',[0 34 211 2],'Style','frame');

uicontrol('Parent',h1,...
'Callback',{@classificationFig_uicallback,h1,'radio_supervised_Callback'},...
'FontName','Helvetica',...
'FontSize',9,...
'Position',[9 67 85 16],...
'String','Supervised',...
'Style','radiobutton',...
'TooltipString','Do supervised image classification',...
'Value',1,...
'Tag','radio_supervised');

uicontrol('Parent',h1,...
'Callback',{@classificationFig_uicallback,h1,'radio_unsupervised_Callback'},...
'FontName','Helvetica',...
'FontSize',9,...
'Position',[108 67 98 16],...
'String','Unsupervised',...
'Style','radiobutton',...
'TooltipString','Do unsupervised image classification',...
'Tag','radio_unsupervised');

function classificationFig_uicallback(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
feval(callback_name,hObject,[],guidata(h1));

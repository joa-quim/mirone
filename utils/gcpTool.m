function handles = gcpTool(handles,axis_t,X,Y,I)
% Transform a Mirone window into a Ground Control Points selection and
% Image registration tool. If the master image has coordinates this will work
% as Image-to-Map rectification, otherwise it works in the Image-to-Image mode.

%	Copyright (c) 2004-2006 by J. Luis
%
%	This program is free software; you can redistribute it and/or modify
%	it under the terms of the GNU General Public License as published by
%	the Free Software Foundation; version 2 of the License.
%
%	This program is distributed in the hope that it will be useful,
%	but WITHOUT ANY WARRANTY; without even the implied warranty of
%	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%	GNU General Public License for more details.
%
%	Contact info: w3.ualg.pt/~jluis/mirone
% --------------------------------------------------------------------


delete(findobj(handles.figure1,'Tag','NewFigure'));     delete(findobj(handles.figure1,'Tag','ImportGMTgrid'))
delete(findobj(handles.figure1,'Tag','SaveGMTgrid'));   delete(findobj(handles.figure1,'Tag','Preferences'))
delete(findobj(handles.figure1,'Tag','Print'));         delete(findobj(handles.figure1,'Tag','DrawText'))
delete(findobj(handles.figure1,'Tag','DrawGeogCirc'));  delete(findobj(handles.figure1,'Tag','DrawLine'))
delete(findobj(handles.figure1,'Tag','DrawRect'));      delete(findobj(handles.figure1,'Tag','DrawPolyg'))
delete(findobj(handles.figure1,'Tag','DrawArrow'));     delete(findobj(handles.figure1,'Tag','Tesoura'))
delete(findobj(handles.figure1,'Tag','ColorPal'));      delete(findobj(handles.figure1,'Tag','MeasureAzim'))
delete(findobj(handles.figure1,'Tag','Shading'));       delete(findobj(handles.figure1,'Tag','Anaglyph'))
delete(findobj(handles.figure1,'Tag','TerrainMod'));    delete(findobj(handles.figure1,'Tag','MBplaning'))
delete(findobj(handles.figure1,'Tag','FlederPlanar'));  delete(findobj(handles.figure1,'Tag','ImageInfo'))
delete(findobj(handles.figure1,'Tag','Refresh'));

delete(findobj(handles.figure1,'Tag','Image'));         delete(findobj(handles.figure1,'Tag','Tools'))
delete(findobj(handles.figure1,'Tag','Draw'));          delete(findobj(handles.figure1,'Tag','Geophysics'))
delete(findobj(handles.figure1,'Tag','Help'));          delete(findobj(handles.figure1,'Tag','GridTools'))
delete(findobj(handles.figure1,'Tag','Atlas'));         delete(findobj(handles.figure1,'Tag','MagBar'));
delete(findobj(handles.figure1,'Tag','ExtDB'));

% ------------- Cleverer deletion of unwanted uicontrols
h_File = findobj(handles.figure1,'Tag','File');
h1 = findobj(handles.figure1,'Tag','OpenGI');
h2 = get(h_File,'Children');
h3 = setxor(h2,h1);
delete(h3)

h_DataSets = findobj(handles.figure1,'Tag','Datasets');
h2 = get(h_DataSets,'Children');
h21 = findobj(handles.figure1,'Tag','VoidDatasetsCoastLine');
h22 = findobj(handles.figure1,'Tag','VoidDatasetsPB');
h23 = findobj(handles.figure1,'Tag','VoidDatasetsRivers');
h3 = setxor(h2,[h21; h22; h23]);
delete(h3)
% -------------

% Import icons
load ([pwd filesep 'data' filesep 'mirone_icons.mat']);
h_toolbar = findobj(handles.figure1,'Tag','FigureToolBar');
uitoggletool('parent',h_toolbar,'Click',{@InsertPoint_Callback,handles}, ...
   'Tag','singlePoint','cdata',point_ico,'TooltipString','Insert point','Separator','on');
uitoggletool('parent',h_toolbar,'Click',{@InsertPointAndPred_Callback,handles}, ...
   'Tag','PointPred','cdata',single(pointPred_ico),'TooltipString','Insert point and prediction');
uipushtool('parent',h_toolbar,'Click',{@registerImage_Callback,handles},...
    'cdata',R_ico,'TooltipString','Register','Separator','on');

handles.origFig = [];       % We don't need the image copy anymore

posStatusBar = get(findall(handles.figure1,'Type','axes','Tag','sbAxes'),'Pos');
if (strcmp(get(handles.axes1,'Visible'),'on'))
    margAnotVer = 30;      % Inventado
    margAnotHor = 15;      % Inventado
else
    margAnotVer = 0;
    margAnotHor = 0;
end
sldT = 13;          % Slider thickness
marg = 5;           % Horizontal margin between the axes(+slider) and the figure
topMargin = 30;     % To accomudate the uis on top of the images
windowsBar = 60;    % In fact I don't know
LeastImageWidth = 360;  % Working value for a 800 pixels width screen
screen = get(0, 'ScreenSize');
MaxImageWidth = screen(3)/2 - (2*sldT+3*marg + 2*margAnotVer);  % Maximum image width supported in this screen

dims1 = size(get(handles.grd_img,'CData'));
img1H = dims1(1);       img1W = dims1(2);
dims2 = size(I);
img2H = dims2(1);       img2W = dims2(2);

% ------------- Don't make images too small
if (img1W < LeastImageWidth)
    img1H = img1H * LeastImageWidth / img1W;    % Rescale height of img1
    img1W = LeastImageWidth;
elseif (img1W > MaxImageWidth)
    img1H = img1H * MaxImageWidth / img1W;      % Rescale height of img1    
    img1W = MaxImageWidth;
end

if (img2W < LeastImageWidth)
    img2H = img2H * LeastImageWidth / img2W;    % Rescale height of img2
    img2W = LeastImageWidth;
elseif (img2W > MaxImageWidth)
    img2H = img2H * MaxImageWidth / img2W;      % Rescale height of img2    
    img2W = MaxImageWidth;
end
% ---------------

% --------------- Compute new figure size
if (img1W + img2W <= screen(3)+2*sldT+3*marg)
    figH = min(screen(4)-60,min(img1H,img2H)+posStatusBar(end)+margAnotHor) + topMargin + sldT;
    posFig = [0 0 (img1W + img2W + 2*sldT+3*marg + margAnotVer) figH];
else
    posFig = [0 0 screen(3) screen(4)-60];
end

axesW = min(img1W,img2W);
axesH = min(min(img1H,img2H),screen(4)-windowsBar);
aspectImg = axesH / axesW;

set(handles.figure1,'Visible','off')
set(handles.figure1,'Pos',posFig,'Name','GCP Tool')     % Resize the figure
movegui(handles.figure1,'north')                        % And reposition it

% --------------- Reposition the existing axes
set(handles.axes1,'Units','pixels','Tag','axes1')
%set(handles.axes1,'Pos',[(posFig(3)-axesW-sldT-marg) (posStatusBar(4)+margAnotHor+sldT+5) axesW axesH])
set(handles.axes1,'Pos',[(posFig(3)-img1W-sldT-marg) (posStatusBar(4)+margAnotHor+sldT+5) img1W axesH])
xlim = get(handles.axes1,'xlim');        ylim = get(handles.axes1,'ylim');
setappdata(handles.axes1,'ThisImageLims',[xlim ylim])
aspectPixeis = img1H / img1W;
aspectData = dims1(1) / dims1(2);
if (abs(aspectPixeis - aspectImg) > 1e-3)      % This axes was distorted
    if (aspectImg < 1)
        aspectThis = axesH / img1H;
        set(handles.axes1,'ylim',ylim(1)+[0 (ylim(2)-ylim(1))*aspectThis])
    else
        aspectThis = axesW / img1W;
        set(handles.axes1,'xlim',xlim(1)+[0 (xlim(2)-xlim(1))*aspectThis])
    end
end
handles.hMasterImage = handles.grd_img;
% ---------------

% ------------- Create the second axes
fonteName = get(handles.axes1,'FontName');      fonteSize = get(handles.axes1,'FontSize');
handles.axes2 = axes('Units','pixels','FontName',fonteName,'FontSize',fonteSize, ...
    'Pos',[(posFig(1)+marg) (posStatusBar(end)+margAnotHor+sldT+5) img2W axesH]);
handles.hSlaveImage = image(X,Y,I,'Parent',handles.axes2);    set(handles.axes2,'Tag','axes2')
if (strcmp(axis_t,'xy')),           set(handles.axes2,'XDir','normal','YDir','normal')
elseif (strcmp(axis_t,'off')),      set(handles.axes2,'Visible','off')
end
xlim = get(handles.axes2,'xlim');        ylim = get(handles.axes2,'ylim');
setappdata(handles.axes2,'ThisImageLims',[xlim ylim])

aspectPixeis = img2H / img2W;
aspectData = dims2(1) / dims2(2);
if (abs(aspectPixeis - aspectImg) > 1e-3)      % This axes was distorted
    if (aspectImg < 1)
        aspectThis = axesH / img2H;
        set(handles.axes2,'ylim',ylim(1)+[0 (ylim(2)-ylim(1))*aspectThis])
    else
        aspectThis = axesW / img2W;
        set(handles.axes2,'xlim',xlim(1)+[0 (xlim(2)-xlim(1))*aspectThis])
    end
end
% ----------------

% ---------------- Create the sliders
pos1 = get(handles.axes1,'Pos');
handles.slider1_ver = uicontrol('Units','pixels','Style','slider','Pos',[pos1(1)+pos1(3)+1 pos1(2) sldT pos1(4)+1]);
handles.slider1_hor = uicontrol('Units','pixels','Style','slider','Pos',[pos1(1) pos1(2)-sldT-1-margAnotHor pos1(3)+1 sldT]);
set(handles.slider1_hor,'Min',0,'Max',1,'Value',0,'Tag','HOR','Callback',{@slider_Cb,handles.axes1,'SetSliderHor'})
set(handles.slider1_ver,'Min',0,'Max',1,'Value',0,'Tag','VER','Callback',{@slider_Cb,handles.axes1,'SetSliderVer'})
% Register the sliders in the axe's appdata
setappdata(handles.axes1,'SliderAxes',[handles.slider1_hor handles.slider1_ver])
imscroll_j(handles.axes1,'ZoomSetSliders')              % ...
%
pos2 = get(handles.axes2,'Pos');
handles.slider2_ver = uicontrol('Units','pixels','Style','slider','Pos',[pos2(1)+pos2(3)+1 pos2(2) sldT pos2(4)+1]);
handles.slider2_hor = uicontrol('Units','pixels','Style','slider','Pos',[pos2(1) pos2(2)-sldT-1-margAnotHor pos2(3)+1 sldT]);
set(handles.slider2_hor,'Min',0,'Max',1,'Value',0,'Tag','HOR','Callback',{@slider_Cb,handles.axes2,'SetSliderHor'})
set(handles.slider2_ver,'Min',0,'Max',1,'Value',0,'Tag','VER','Callback',{@slider_Cb,handles.axes2,'SetSliderVer'})
% Register the sliders in the axe's appdata
setappdata(handles.axes2,'SliderAxes',[handles.slider2_hor handles.slider2_ver])
imscroll_j(handles.axes2,'ZoomSetSliders')              % ...

% ---------------- Create the popups
handles.popup_transf = uicontrol('Units','Pixels','TooltipString','Type of transformation',...
    'String',{'Affine'; 'Linear conformal'; 'Projective'; 'Polynomial (6 pts)';...
        'Polynomial (10 pts)'; 'Polynomial (16 pts)'; 'Piecewise linear'; 'Loc weighted mean'}, ...
    'Pos',[pos2(1)+pos2(3)-130 pos2(2)+pos2(4)+2 120 21],'Style','popupmenu', 'BackgroundColor','w');

handles.popup_interpMethod = uicontrol('Units','Pixels','TooltipString','Specifies the form of interpolation to use.',...
    'String',{'bilinear'; 'bicubic'; 'nearest';}, 'BackgroundColor','w', ...
    'Pos',[pos1(1)+10 pos1(2)+pos1(4)+2 90 21],'Style','popupmenu');
% ----------------

% ---------------- Create the text image names
uicontrol('Units','Pixels','Style','text','String','Slave Image','FontName','Helvetica','FontSize',10, ...
    'Pos',[pos2(1) pos2(2)+pos2(4)+2 75 16],'HorizontalAlignment','left')
uicontrol('Units','Pixels','Style','text','String','Master Image','FontName','Helvetica','FontSize',10, ...
    'Pos',[pos1(1)+pos1(3)-90 pos1(2)+pos1(4)+2 90 16],'HorizontalAlignment','right')
% ----------------

handles.slavePoints = [];
handles.masterPoints = [];
handles.count = 0;
handles.isCoupled = [];

hand_prev = guidata(handles.figure1);   % The geog type of the Master image is still in appdata
handles.geog = hand_prev.geog;          % but the handles.geog contains the geog state of Slave img
set(handles.figure1,'Visible','on','pointer','arrow')
guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------------------
function slider_Cb(obj,evt,ax,opt)
imscroll_j(ax,opt)

% -----------------------------------------------------------------------------------------
function InsertPoint_Callback(hObject,event,handles)
handles = guidata(handles.figure1);
but = 1;    count = handles.count;
while (but == 1)
    lastClickedAx = getappdata(handles.figure1,'clickedAx');
    [x,y,but] = ginput_pointer(1,'crosshair');
    if (but ~= 1)                                       % Stop insertion. Do eventual cleaning
        wichAxes = get(get(handles.figure1,'CurrentAxes'),'Tag');
        if (~handles.isCoupled)                % Unmatched point point
            delete(handles.hLastPt)
        end
        if (count > 0 && size(handles.slavePoints,1) ~= size(handles.masterPoints,1))
            if (size(handles.slavePoints,1) > size(handles.masterPoints,1))
                handles.slavePoints(end,:) = [];
            else
                handles.masterPoints(end,:) = [];
            end
            count = count - 1;
        end
        setappdata(handles.figure1,'clickedAx',[]);     % Reset the point insertion machine
        %handles.isCoupled = 0;
        set(findobj(handles.figure1,'Tag','singlePoint'),'State','off')
        break
    end
    
    wichAxes = get(get(handles.figure1,'CurrentAxes'),'Tag');
    if (isempty(lastClickedAx))                         % First GCP on this run
        handles.hLastPt = line(x,y,'Marker','o','MarkerFaceColor','y',...
            'MarkerEdgeColor','k','MarkerSize',7,'Tag','GCPSymbol','UserData',count+1);
        handles.isCoupled = 0;
        if (wichAxes(end) == '1')
            handles.masterPoints = [handles.masterPoints; x y];
        else
            handles.slavePoints = [handles.slavePoints; x y];
        end
        count = count + 1;
        handles.count = count;
        set_gcp_uicontext(handles,handles.hLastPt)
    else                                                % Next GCPs
        if (lastClickedAx(end) ~= wichAxes(end))        % Conjugated point, OR ...
            handles.hLastPt = line(x,y,'Marker','o','MarkerFaceColor','y',...
                'MarkerEdgeColor','k','MarkerSize',7,'Tag','GCPSymbol');
            if (wichAxes(end) == '1')
                handles.masterPoints = [handles.masterPoints; x y];
            else
                handles.slavePoints = [handles.slavePoints; x y];
            end
            if (handles.isCoupled == 1)        % Start new pair on the other axes
                handles.isCoupled = 0;
                count = count + 1;
            else                               % Finish pair OR non left-click on the other axes
                handles.isCoupled = 1;
            end
            set(handles.hLastPt,'UserData',count)
            handles.count = count;
            set_gcp_uicontext(handles,handles.hLastPt)
        else                                    % Click on the same axes. Reset point
            set(handles.hLastPt,'XData',x,'YData',y)
            if (wichAxes(end) == '1'),  handles.masterPoints(count,:) = [x y];
            else                        handles.slavePoints(count,:) = [x y];
            end
        end
    end
    setappdata(handles.figure1,'clickedAx',wichAxes)
end
guidata(handles.figure1,handles)

% ---------------------------------------------------------------------------
function InsertPointAndPred_Callback(hObject,event,handles)
handles = guidata(handles.figure1);
tipo = checkTransform(handles);
if (isempty(tipo)),     return;    end      % Error message already issued
but = 1;    count = handles.count;
while (but == 1)
    %lastClickedAx = getappdata(handles.figure1,'clickedAx');
    [x,y,but] = ginput_pointer(1,'crosshair');
    if (but ~= 1),  break;  end
    wichAxes = get(get(handles.figure1,'CurrentAxes'),'Tag');
    if (wichAxes(end) == '1')       % Master -> Slave prediction
        trf = transform_fun('cp2tform',handles.masterPoints, handles.masterPoints, tipo);
        [x_pred,y_pred] = transform_fun('tformfwd',trf,[handles.masterPoints(:,1); x],...
            [handles.masterPoints(:,2); y]);
        handles.masterPoints = [handles.masterPoints; x y];
        handles.slavePoints = [handles.slavePoints; x_pred(end) y_pred(end)];
        h_pt = line(x,y,'Marker','o','MarkerFaceColor','y',...
            'MarkerEdgeColor','k','MarkerSize',7,'Tag','GCPSymbol','UserData',count+1);
        handles.hLastPt = line(x_pred(end),y_pred(end),'Parent',handles.axes2,'Marker','o','MarkerFaceColor','y',...
            'MarkerEdgeColor','k','MarkerSize',7,'Tag','GCPSymbol','UserData',count+1);
    else                            % Slave -> Master prediction
        trf = transform_fun('cp2tform',handles.slavePoints, handles.masterPoints, tipo);
        [x_pred,y_pred] = transform_fun('tformfwd',trf,[handles.slavePoints(:,1); x],...
            [handles.slavePoints(:,2); y]);
        handles.slavePoints = [handles.slavePoints; x y];
        handles.masterPoints = [handles.masterPoints; x_pred(end) y_pred(end)];
        h_pt = line(x,y,'Marker','o','MarkerFaceColor','y',...
            'MarkerEdgeColor','k','MarkerSize',7,'Tag','GCPSymbol','UserData',count+1);
        handles.hLastPt = line(x_pred(end),y_pred(end),'Parent',handles.axes1,'Marker','o','MarkerFaceColor','y',...
            'MarkerEdgeColor','k','MarkerSize',7,'Tag','GCPSymbol','UserData',count+1);
    end
    count = count + 1;
    handles.count = count;
    set_gcp_uicontext(handles,handles.hLastPt)
    set_gcp_uicontext(handles,h_pt)
end
set(findobj(handles.figure1,'Tag','PointPred'),'State','off')
guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------------------
function set_gcp_uicontext(handles,h)
% Set uicontexts for the symbols. h is the handle to the marker (line in fact) object

cmenuHand = uicontextmenu;      set(h, 'UIContextMenu', cmenuHand);
ui_edit_polygon(h)    % Set edition functions
uimenu(cmenuHand, 'Label', 'Remove GCP pair', 'Callback', {@remove_pair,handles,h});
uimenu(cmenuHand, 'Label', 'Remove all GCPs', 'Callback', {@removeALLgcps,handles,h});
uimenu(cmenuHand, 'Label', 'View GCPs table', 'Callback', {@show_gcp,handles,h});
uimenu(cmenuHand, 'Label', 'Show GCP numbers', 'Callback', {@showGCPnumbers,handles},'Tag','GCPlab');
% uimenu(cmenuHand, 'Label', 'Try fine tune GCPs', 'Callback', {@doCPcorr,handles},'Separator','on');
uimenu(cmenuHand, 'Label', 'Save GCPs', 'Callback', {@doWriteGCPs,handles},'Separator','on');

% -----------------------------------------------------------------------------------------
function remove_pair(event,obj,handles,h)
% Remove a pair of conjugated GCPs
handles = guidata(handles.figure1);
id = get(h,'UserData');
hh = findobj(handles.figure1,'Type','line','UserData',id);
delete(hh)
np = size(handles.masterPoints,1);
for (i = id:np)
    hh = findobj(handles.figure1,'Type','line','UserData',i);
    set(hh,'UserData',i-1)
end
handles.masterPoints(id,:) = [];
handles.slavePoints(id,:) = [];
handles.count = handles.count - 1;
guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------------------
function removeALLgcps(event,obj,handles,h)
    delete(findobj(handles.figure1,'Tag','GCPSymbol'))
    handles = guidata(handles.figure1);
    handles.masterPoints = [];
    handles.slavePoints = [];
    handles.count = 0;
    handles.isCoupled = [];
    guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------------------
function handles = getUpdatedCPs(handles)
% We need to update the handles because points may have been mouse edited

hSlaves = findobj(handles.axes2,'Type','line','Tag','GCPSymbol');
hMasters = findobj(handles.axes1,'Type','line','Tag','GCPSymbol');
ordem = get(hMasters,'Userdata');
xSlaves = get(hSlaves,'Xdata');     ySlaves = get(hSlaves,'Ydata');
xMasters = get(hMasters,'Xdata');   yMasters = get(hMasters,'Ydata');

% I think this must be set everywhere when IMG1 has other coordinates than pixeis
% x_lim = get(handles.axes1,'XLim');
% [img1H,img1W,k] = size(get(handles.hSlaveImage,'CData'));
% if ( ~isequal(x_lim(2), img1W) )        % Slave points must be in pixel units
%     xSlaves = axes2pix(img1W,get(handles.hSlaveImage,'XData'),xSlaves);
%     ySlaves = axes2pix(img1H,get(handles.hSlaveImage,'YData'),ySlaves);
% end

if (iscell(ordem))          % Almost allways. When number of GCPs > 1
    ordem = cell2mat(ordem);
    xSlaves = cell2mat(xSlaves);        ySlaves = cell2mat(ySlaves);
    xSlaves = xSlaves(ordem);           ySlaves = ySlaves(ordem);
    xMasters = cell2mat(xMasters);      yMasters = cell2mat(yMasters);
    xMasters = xMasters(ordem);         yMasters = yMasters(ordem);
end

handles.slavePoints = [xSlaves ySlaves];
handles.masterPoints = [xMasters yMasters];
guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------------------
function show_gcp(event,obj,handles,h)
% Display the GCPs in tableGUI

handles = guidata(handles.figure1);
handles = getUpdatedCPs(handles);

[tipo,to_order] = checkTransform(handles);
if (isempty(tipo)),     return;    end      % Error message already issued

if (strncmp(tipo,'poly',4))
    tipo = {'polynomial' to_order-2};   % to_order-2 is a trick to get the right order from popup list order
else
    tipo = {tipo};
end

try
    trf = transform_fun('cp2tform',handles.slavePoints,handles.masterPoints,tipo{:});
    %trf = cp2tform(handles.slavePoints,handles.masterPoints,tipo{:});
catch
    errordlg(lasterr,'Error');    trf = [];     return    
end

[x,y] = transform_fun('tformfwd',trf,handles.slavePoints(:,1),handles.slavePoints(:,2));
if (handles.geog)
    residue = vdist_vectorized(handles.masterPoints(:,2),handles.masterPoints(:,1), y, x);
    str_res = 'Residue (m)';
else
    residue = [handles.masterPoints(:,1) handles.masterPoints(:,2)] - [x y];
    residue = sqrt(residue(:,1).^2 + residue(:,2).^2);
    str_res = 'Residue (?)';
end
gcp = [handles.masterPoints handles.slavePoints residue];
FigName = ['GCP Table - ' tipo{:} ' - Total Res = ' num2str(sum(residue))];
out = tableGUI('array',gcp,'RowNumbers','y','ColNames',{'Master Points - X','Master Points - Y',...
        'Slave Points - X','Slave Points - Y',str_res},'ColWidth',100,'FigName',FigName,'modal','');

% ---------------------------------------------------------------------------
function trf = registerImage_Callback(hObject,event,handles)
handles = guidata(handles.figure1);
[tipo,to_order] = checkTransform(handles);
if (isempty(tipo)),     return;    end      % Error message already issued
handles = getUpdatedCPs(handles);           % Get update GCPs positions

interpola = get(handles.popup_interpMethod,'String');
interpola = interpola{get(handles.popup_interpMethod,'Value')};

if (strncmp(tipo,'poly',4))
    tipo = {'polynomial' to_order-2};   % to_order-2 is a trick to get the right order from popup list order
else
    tipo = {tipo};
end

try
    trf = transform_fun('cp2tform',handles.slavePoints,handles.masterPoints,tipo{:});
catch
    errordlg(lasterr,'Error');    trf = [];     return    
end

if (nargout == 0)   % Register the image
    lims = getappdata(handles.axes2,'ThisImageLims');
	[reg,X,Y] = transform_fun('imtransform',get(handles.hSlaveImage,'CData'),trf,...
                             interpola,'UData',lims(1:2), 'VData',lims(3:4),'size', ...
                             size(get(handles.hSlaveImage,'CData')));
% 	[reg,X,Y] = transform_fun('imtransform',get(handles.hSlaveImage,'CData'),trf,...
%         interpola,'size',size(get(handles.hSlaveImage,'CData')));

    if (handles.image_type == 2)        % Image witout coords
        if (ndims(reg) == 2)
            setappdata(0,'CropedColormap',get(handles.figure1,'Colormap'))
        end
        mirone(reg);
    else
        tmp.geog = handles.geog;
        tmp.X = X;        tmp.Y = Y;
        tmp.head(1:4) = [X Y];        tmp.head(5:7) = [0 255 1];
        tmp.head(8) = diff(X) / (size(reg,2) - 1);
        tmp.head(9) = diff(Y) / (size(reg,1) - 1);
        tmp.name = 'Registered Image';
        if (ndims(reg) == 2)
            tmp.cmap = get(handles.figure1,'Colormap');
        end
        mirone(reg,tmp);
    end
else                % Just return the tform structure
    %
end

% ---------------------------------------------------------------------------
function [tipo,transf] = checkTransform(handles)
% Check that the number of points is enough for the selected transform

transf = get(handles.popup_transf,'Value');
switch transf
    case 1,     tipo = 'affine';
    case 2,     tipo = 'linear conformal';
    case 3,     tipo = 'projective';
    case 4,     tipo = 'polynomial (6 pts)';
    case 5,     tipo = 'polynomial (10 pts)';
    case 6,     tipo = 'polynomial (16 pts)';
    case 7,     tipo = 'piecewise linear';
    case 8,     tipo = 'lwm';
end
n_cps = size(handles.slavePoints,1);
msg = '';
if (transf == 1 && n_cps < 3)
    msg = 'Minimum Control points for affine transform is 3.';
elseif (transf == 2 && n_cps < 2)
    msg = 'Minimum Control points for Linear conformal transform is 2.';
elseif (transf == 3 && n_cps < 4)
    msg = 'Minimum Control points for projective transform is 4.';
elseif (transf == 4 && n_cps < 6)
    msg = 'Minimum Control points for polynomial order 2 transform is 6.';
elseif (transf == 5 && n_cps < 6)
    msg = 'Minimum Control points for polynomial order 3 transform is 10.';
elseif (transf == 6 && n_cps < 6)
    msg = 'Minimum Control points for polynomial order 2 transform is 16.';
elseif (transf == 7 && n_cps < 4)
    msg = 'Minimum Control points for piecewise linear transform is 4.';
elseif (transf == 8 && n_cps < 12)
    msg = 'Minimum Control points for Locolal weightd mean transform is 12.';
end
if (~isempty(msg))
    errordlg([msg ' Either select more CPs or choose another trasform type'],'Error');    tipo = [];
end

%-----------------------------------------------------------------------------------
function doWriteGCPs(event,obj,handles)
% This function is ...
    handles = guidata(handles.figure1);
    str1 = {'*.dat;*.DAT', 'Control points file (*.dat,*.DAT)';'*.*', 'All Files (*.*)'};
    [FileName,PathName] = uiputfile(str1,'GCP file name');
    if isequal(FileName,0);     return;     end
	%Open and write to ASCII file
	if ispc;        fid = fopen([PathName FileName],'wt');
	elseif isunix;  fid = fopen([PathName FileName],'w');
	else    errordlg('Unknown platform.','Error');  return;
	end
	fprintf(fid,'# X(slave)\tY(slave)\tX(master)\tY(master)\n');
	fprintf(fid,'%f\t%f\t%f\t%f\n', [handles.slavePoints handles.masterPoints]');
    fclose(fid);

% %-----------------------------------------------------------------------------------
% function doCPcorr(event,obj,handles)
% 
% hSlaves = findobj(handles.axes2,'Type','line','Tag','GCPSymbol');
% hMasters = findobj(handles.axes1,'Type','line','Tag','GCPSymbol');
% ordem = get(hMasters,'Userdata');
% xSlaves = get(hSlaves,'Xdata');     ySlaves = get(hSlaves,'Ydata');
% xMasters = get(hMasters,'Xdata');   yMasters = get(hMasters,'Ydata');
% 
% if (iscell(ordem))          % Almost allways. When number of GCPs > 1
%     ordem = cell2mat(ordem);
%     xSlaves = cell2mat(xSlaves);        ySlaves = cell2mat(ySlaves);
%     xMasters = cell2mat(xMasters);      yMasters = cell2mat(yMasters);
% end
% 
% slaves_corr = cpcorr_j([xSlaves ySlaves],[xMasters yMasters], ...
%     get(handles.hSlaveImage,'CData'),get(handles.hMasterImage,'CData'));
% 
% diff_pos = sum([xSlaves ySlaves] - slaves_corr, 2);
% I = (diff_pos == 0);
% if (~any(I))
%     warndlg('Could not refine the GCPs positions.','Warning')
%     return
% end
% 
% % Update to the adjusted positions
% for (i = 1:numel(xSlaves))
%     if (I(i)),      continue;   end     % This point has not change
%     set(hSlaves(ordem(i)),'XData',slaves_corr(i,1),'YData',slaves_corr(i,2))
% end

%-----------------------------------------------------------------------------------
function showGCPnumbers(obj,event,handles)
% Plot/desplot the GCPs numbers
handles = guidata(handles.figure1);

if (strcmp(get(obj,'Label'),'Show GCP numbers'))
	% We need to update the handles here because points may have been mouse edited
	hSlaves = findobj(handles.axes2,'Type','line','Tag','GCPSymbol');
	hMasters = findobj(handles.axes1,'Type','line','Tag','GCPSymbol');
	ordem = get(hMasters,'Userdata');
	xSlaves = get(hSlaves,'Xdata');     ySlaves = get(hSlaves,'Ydata');
	xMasters = get(hMasters,'Xdata');   yMasters = get(hMasters,'Ydata');
	
	if (iscell(ordem))          % Almost allways. When number of GCPs > 1
        ordem = cell2mat(ordem);
        xSlaves = cell2mat(xSlaves);        ySlaves = cell2mat(ySlaves);
        xSlaves = xSlaves(ordem);           ySlaves = ySlaves(ordem);
        xMasters = cell2mat(xMasters);      yMasters = cell2mat(yMasters);
        xMasters = xMasters(ordem);         yMasters = yMasters(ordem);
	end
    ordem = ordem(ordem);

	% Estimate the text position shift in order that it doesn't fall over the symbols
    dpis = get(0,'ScreenPixelsPerInch') ;   % screen DPI
    symb_size = 7 / 72 * 2.54;              % Symbol size in cm (circles size is 7 points)
	n_texts = numel(xSlaves);
    
	pos = get(handles.axes1,'Position');    ylim = get(handles.axes1,'Ylim');
    escala = diff(ylim)/(pos(4)*2.54/dpis); % Image units / cm
    dy = symb_size * escala;
    
	handles.hTextMaster = zeros(1,n_texts);
	for i = 1:n_texts
        handles.hTextMaster(i) = text(xMasters(i),yMasters(i)+dy,0,num2str(ordem(i)),'Fontsize',8,'Parent',handles.axes1);
	end
    
	pos = get(handles.axes2,'Position');    ylim = get(handles.axes2,'Ylim');
    escala = diff(ylim)/(pos(4)*2.54/dpis); % Image units / cm
    dy = symb_size * escala;
    
	handles.hTextSlave = zeros(1,n_texts);
	for i = 1:n_texts
        handles.hTextSlave(i) = text(xSlaves(i),ySlaves(i)+dy,0,num2str(ordem(i)),'Fontsize',8,'Parent',handles.axes2);
	end
    
    % Change the uimenus labels to "Hide"
    hM = get(hMasters,'uicontextmenu');
    hS = get(hSlaves,'uicontextmenu');
    if (n_texts > 1)
        for i = 1:n_texts
            set(findobj(hM{i},'Tag','GCPlab'),'Label','Hide GCP numbers')
            set(findobj(hS{i},'Tag','GCPlab'),'Label','Hide GCP numbers')
        end
    else
        set(findobj(hM,'Tag','GCPlab'),'Label','Hide GCP numbers')
        set(findobj(hS,'Tag','GCPlab'),'Label','Hide GCP numbers')
    end
else
    delete(handles.hTextMaster);        handles.hTextMaster = [];
    delete(handles.hTextSlave);         handles.hTextSlave = [];
	hSlaves = findobj(handles.axes2,'Type','line','Tag','GCPSymbol');
	hMasters = findobj(handles.axes1,'Type','line','Tag','GCPSymbol');
    
    % Change the uimenus labels to "Hide"
    hM = get(hMasters,'uicontextmenu');
    hS = get(hSlaves,'uicontextmenu');
    if (numel(hMasters) > 1)
        for i = 1:numel(hMasters)
            set(findobj(hM{i},'Tag','GCPlab'),'Label','Show GCP numbers')
            set(findobj(hS{i},'Tag','GCPlab'),'Label','Show GCP numbers')
        end
    else
        set(findobj(hM,'Tag','GCPlab'),'Label','Show GCP numbers')
        set(findobj(hS,'Tag','GCPlab'),'Label','Show GCP numbers')
    end
end
guidata(handles.figure1,handles)

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

delete(handles.NewFigure);      delete(handles.ImportKnownTypes)
delete(handles.SaveGMTgrid);    delete(handles.Preferences)
delete(handles.Print);          delete(handles.DrawText)
delete(handles.DrawGeogCirc);   delete(handles.DrawLine)
delete(handles.DrawRect);       delete(handles.DrawPolyg)
delete(handles.DrawArrow);      set(handles.Tesoura,'Enable','off') % cannot kill coze test in PanZoom
delete(handles.ColorPal);       delete(handles.Shading);
delete(handles.Anaglyph);       delete(handles.toGE)
delete(handles.MBplaning);      delete(handles.FlederPlanar);
delete(handles.ImageInfo);      delete(handles.Refresh);
delete(handles.Image);          delete(handles.Tools)
delete(handles.Draw);           delete(handles.Geophysics)
delete(handles.Help);           delete(handles.GridTools)
%delete(handles.TerrainMod);
if (ishandle(handles.Projections)),     delete(handles.Projections);   end

% ------------- Cleverer deletion of unwanted uicontrols
h1 = get(handles.File,'Children');
h2 = setxor(h1,handles.OpenGI);
delete(h2)

h2 = get(handles.Datasets,'Children');
h21 = handles.VoidDatasetsCoastLine;
h22 = handles.VoidDatasetsPB;
h23 = handles.VoidDatasetsRivers;
h3 = setxor(h2,[h21; h22; h23]);
delete(h3)
% -------------

% Import icons
load ([pwd filesep 'data' filesep 'mirone_icons.mat']);
h_toolbar = handles.FigureToolBar;
uitoggletool('parent',h_toolbar,'Click',{@InsertPoint_Callback,handles}, ...
   'Tag','singlePoint','cdata',point_ico,'TooltipString','Insert point','Separator','on');
uitoggletool('parent',h_toolbar,'Click',{@InsertPointAndPred_Callback,handles}, ...
   'Tag','PointPred','cdata',pointPred_ico,'TooltipString','Insert point and prediction');
uipushtool('parent',h_toolbar,'Click',{@registerImage_Callback,handles},...
    'cdata',R_ico,'TooltipString','Register','Separator','on');
uipushtool('parent',h_toolbar,'Click',{@registerSIFT_Callback,handles},...
    'cdata',cerejas,'TooltipString','Automatic Registration','Separator','on');

handles.origFig = [];       % We don't need the image copy anymore

posStatusBar = get(findall(handles.figure1,'Type','axes','Tag','sbAxes'),'Pos');
if (strcmp(get(handles.axes1,'Visible'),'on'))
    margAnotVer = 30;      % Inventado
    margAnotHor = 15;      % Inventado
else
    margAnotVer = 0;
    margAnotHor = 0;
end
sldT = 13;              % Slider thickness
marg = 5;               % Horizontal margin between the axes(+slider) and the figure
topMargin = 30;         % To accomudate the uis on top of the images
windowsBar = 60;        % In fact I don't know
LeastImageWidth = 360;  % Working value for a 800 pixels width screen
screen = get(0, 'ScreenSize');
MaxImageWidth = screen(3)/2 - (2*sldT+3*marg + 2*margAnotVer);  % Maximum image width supported in this screen

dims1 = size(get(handles.hImg,'CData'));
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
tmp = getappdata(handles.axes1,'ThisImageLims');
xlim = tmp(1:2);    ylim = tmp(3:4);
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
handles.hMasterImage = handles.hImg;
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

% ------------- Here we must check if both images are indexed and with different cmaps
if (ndims(I) == 2 && ndims(get(handles.hImg,'cdata')) == 2)
    pal = get(handles.figure1,'ColorMap');          % Current colormap (was set according to SlaveImage)
    d_pal = 1;
    if (numel(pal) == numel(handles.origCmap))
        d_pal = pal - handles.origCmap;
    end
    if (any(d_pal(:)))                               % Color maps difer. The only solution it to convert the
        Imaster = get(handles.hMasterImage,'CData'); % Master image to RGB, so that it won't be changed to the
        Imaster = ind2rgb8(Imaster,handles.origCmap);% Slave's colormap
        set(handles.hMasterImage,'CData',Imaster);
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
handles.Mhead = hand_prev.head;         %                   Same with the header
handles.Shead  = handles.head;          %   Just to maintain name consistency
handles.Mimage_type = hand_prev.image_type;
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
uimenu(cmenuHand, 'Label', 'View GCPs table', 'Callback', {@show_gcp,handles},'Sep','on');
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
function resMod = show_gcp(event,obj,handles, opt)
% Display the GCPs in tableGUI, OR return the residues norm (no table show then)
% OPT is used to NOT call getUpdatedCPs, because that function fishes the data
% from figure. Instead, when OPT ~= [], we trust on data stored in handles

if (nargin == 3),   opt = [];   end
handles = guidata(handles.figure1);
if (isempty(opt))
    handles = getUpdatedCPs(handles);       % Points may have been edited
end

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
    resMod = vdist_vectorized(handles.masterPoints(:,2),handles.masterPoints(:,1), y, x);
    str_res = 'Residue (m)';
else
    residue = [handles.masterPoints(:,1) handles.masterPoints(:,2)] - [x y];
    resMod = sqrt(residue(:,1).^2 + residue(:,2).^2);
    str_res = 'Residue (?)';
end

if (nargout == 1)       % Just return the norm of the residues
    return
end

gcp = [handles.masterPoints handles.slavePoints resMod];
FigName = ['GCP Table - ' tipo{:} ' - RMS Res = ' num2str( sqrt(sum(resMod.^2))/sqrt(length(resMod)) )];
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

    if (handles.Mimage_type == 2)        % Master Image without coords
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
function [ties,nPts] = registerSIFTautopano_Callback(img_l,img_r,handles)

    set(handles.figure1,'pointer','watch')
    imwrite(img_l,'lixoleft.jpg','Quality',100);
    imwrite(img_r,'lixoright.jpg','Quality',100);

    tmp{1} = '@echo off';
    tmp{2} = 'generatekeys lixoleft.jpg keyfile_left.xml.gz 1500';
    tmp{3} = 'generatekeys lixoright.jpg keyfile_right.xml.gz 1500';
    tmp{4} = 'autopano --maxmatches 50 --disable-areafilter outputa.pto keyfile_right.xml.gz keyfile_left.xml.gz';
    fname = 'corre_autopano.bat';
    fid = fopen(fname,'wt');
    for (i = 1:4),   fprintf(fid,'%s\n',tmp{i});   end
    fclose(fid);
    %dos([fname ' &']);
    dos(fname);
    %system(fname);
    delete('keyfile_left.xml.gz');      delete('keyfile_right.xml.gz');
    delete('lixoright.jpg');            delete('lixoleft.jpg');
    fid=fopen('outputa.pto','r');
    todos = fread(fid,'*char');         fclose(fid);
    txt=strread(todos,'%s');    
    txt(1:61) = [];
    txt(end-2:end) = [];
    nPts = numel(txt) / 8;
    idx = repmat([true(4,1); false(4,1)],nPts,1);
    txt(idx) = [];
    txt = char(txt);
    txt = txt(:,2:end);                 % Uf1, we finaly have numbers only, though they are still strings
    ties = str2num(txt);
    ties = reshape(ties',4,nPts)';      % Uf2, ties is now a nx4 array
    set(handles.figure1,'pointer','arrow')
        
% ---------------------------------------------------------------------------
function trf = registerSIFT_Callback(hObject,event,handles)

    set(handles.figure1,'pointer','watch')
	handles = guidata(handles.figure1);
	imgM = get(handles.hMasterImage,'CData');
	imgS = get(handles.hSlaveImage,'CData');
        
	if (ndims(imgM) == 3),      imgM = cvlib_mex('color',imgM,'rgb2gray');    end
	if (ndims(imgS) == 3),      imgS = cvlib_mex('color',imgS,'rgb2gray');    end
    [mM,nM] = size(imgM);       [mS,nS] = size(imgS);

    [xy_match,nPts] = registerSIFTautopano_Callback(imgS,imgM,handles);

    if (nPts == 0)
        warndlg('Sorry, I was not able to find any matching points between the two images. Quiting','Warning')
        set(handles.figure1,'pointer','arrow')
        return
    end
    
%     % Test if images have coordinates. If they do convert pixel coords to image coords
%     % Fisrst the Master image
    if (  handles.Mimage_type ~= 2 && handles.Mimage_type ~= 20)
        % Convert pixel to x,y coordinates
        x_inc = handles.Mhead(8);    y_inc = handles.Mhead(9);
        x_min = handles.Mhead(1);    y_min = handles.Mhead(3);
        if (handles.Mhead(7))            % Work in grid registration
            x_min = x_min + x_inc/2;
            y_min = y_min + y_inc/2;
        end
        xy_match(:,1) = (xy_match(:,1)-1) * x_inc + x_min;
        xy_match(:,2) = (xy_match(:,2)-1) * y_inc + y_min;
    end
    % And now the Slave image
    if ( handles.image_type ~= 2 )
        % Convert pixel to x,y coordinates
        x_inc = handles.Shead(8);    y_inc = handles.Shead(9);
        x_min = handles.Shead(1);    y_min = handles.Shead(3);
        if (handles.Shead(7))            % Work in grid registration
            x_min = x_min + x_inc/2;
            y_min = y_min + y_inc/2;
        end
        xy_match(:,3) = (xy_match(:,3)-1) * x_inc + x_min;
        xy_match(:,4) = (xy_match(:,4)-1) * y_inc + y_min;
    end
    
    % Compute the residues of the transformation and sort the keypoints
    % based on the them. This way we will only retain starting by the best points.
    handles.masterPoints = [xy_match(:,1) xy_match(:,2)];
    handles.slavePoints  = [xy_match(:,3) xy_match(:,4)];
    handles.count = size(xy_match,1);
    handles.isCoupled = 1;
    guidata(handles.figure1,handles)            % Need to save because it will be retrieved in show_gcp
    if (nPts >= 3)
        resMod = show_gcp([],[],handles,'No');
        disp(['RMS = ' num2str(norm(resMod) / sqrt(handles.count))])
        [resMod,id] = sort(resMod);    
        xy_match = xy_match(id,:);
    end
    
    if (nPts > 50),   xy_match(51:end,:) = [];    end
    nPts   = size(xy_match,1);
    
    x1M = xy_match(:,1);        y1M = xy_match(:,2);
    x1S = xy_match(:,3);        y1S = xy_match(:,4);
    
    handles.masterPoints = [x1M y1M];
    handles.slavePoints  = [x1S y1S];
    handles.isCoupled = 1;
    handles.count = nPts;
        
    delete(findobj(handles.figure1,'Tag','GCPSymbol'))
    
    for (i = 1:nPts)
        h = line(x1M(i),y1M(i),'Parent',handles.axes1,'Marker','o','MarkerFaceColor','y',...
            'MarkerEdgeColor','k','MarkerSize',6,'Tag','GCPSymbol','UserData',i);
        set_gcp_uicontext(handles,h)
        
        h = line(x1S(i),y1S(i),'Parent',handles.axes2,'Marker','o','MarkerFaceColor','y',...
            'MarkerEdgeColor','k','MarkerSize',6,'Tag','GCPSymbol','UserData',i);
        set_gcp_uicontext(handles,h)
    end
    handles.hLastPt = h;
    set(handles.figure1,'pointer','arrow')
    guidata(handles.figure1,handles)
    
% ---------------------------------------------------------------------------
function [x,y] = axes2axes(handles,x,y,opt)
    axMpos = get(handles.axes1,'Pos');
    axSpos = get(handles.axes2,'Pos');
	dimM = [size(get(handles.hMasterImage,'CData'),2) size(get(handles.hMasterImage,'CData'),1)];
	dimS = [size(get(handles.hSlaveImage,'CData'),2) size(get(handles.hSlaveImage,'CData'),1)];
    xx = x;yy=y;
    if (opt == 1)       % x,y are in Master axes
        x = axMpos(1) + x / dimM(1) * axMpos(3);            % Pixels in the figure reference
        y = axMpos(2) + (dimM(2)-y) / dimM(2) * axMpos(4);  % Pixels in the figure reference
        % Now compute pixels relatively to Slave image
        
        rSx  = dimS(1) / axSpos(3);
        Sllf = axSpos(1);
        x    = rSx * (x - Sllf);
        rSy  = dimS(2) / axSpos(4);
        Surf = axSpos(4);
        y    = rSy * (Surf  - y);
    else                % x,y are in Slave axes
        x = axSpos(1) + x / dimS(1) * axSpos(3);    % Pixels relatively to figure
        y = axSpos(2) + (dimS(2)-y) / dimS(2) * axSpos(4);    % Pixels relatively to figure
        % Now compute pixels relatively to Master image
        
        rMx  = dimM(1) / axMpos(3);
        Mllf = axMpos(1);
        x    = rMx * (x - Mllf);
        rMy  = dimM(2) / axMpos(4);
        Murf = axMpos(4);
        y    = rMy * (Murf  - y);
        
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


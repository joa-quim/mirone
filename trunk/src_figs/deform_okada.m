function varargout = deform_okada(varargin)
% M-File changed by desGUIDE 

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
 
hObject = figure('Tag','figure1','Visible','off');
handles = guihandles(hObject);
guidata(hObject, handles);
deform_okada_LayoutFcn(hObject,handles);
handles = guihandles(hObject);
 
movegui(hObject,'east');

handles.h_calling_fig = [];     % Handles to the calling figure
handles.input_locations = [];   % May contain ground points positions
D2R = pi / 180;

if ~isempty(varargin)
    if (length(varargin) >= 8)
        handles = set_all_faults(handles,varargin{:});
        fault_in = 1;
    else
        handles.h_calling_fig = varargin{1};
        handles.h_fault = varargin{2};
        handles.FaultStrike = varargin{3};
        handles.geog = varargin{4};
        fault_in = 0;
    end
else
    % Tenho de prever este caso
end

handles.n_faults = length(handles.h_fault);
if (handles.n_faults > 1)
    S = [];
    s_format = ['%.' num2str(fix(log10(handles.n_faults))+1) 'd'];
    for (i=1:handles.n_faults)      S = [S; ['Fault ' num2str(i,s_format)]];   end
    set(handles.popup_fault,'String',cellstr(S))
    %set(handles.popup_fault,'String',cellstr(num2str(1:handles.n_faults,'%.2d')'))
    set(handles.h_fault(1),'LineStyle','--');   % set the top fault one with a dashed line type
    refresh(handles.h_calling_fig);             % otherwise, ML BUG
else
    set(handles.popup_fault,'Visible','off')
    delete(findobj(hObject,'Style','text','Tag','fault_number'))    % Otherwise it would reborn in Pro look
end

fault_x = get(handles.h_fault,'XData');     fault_y = get(handles.h_fault,'YData');
if (handles.n_faults > 1)
    for (k=1:handles.n_faults)  nvert(k) = size(fault_x{k},2) - 1;  end
else
    nvert = size(fault_x,2) - 1;
end

if (any(nvert > 1))
    set(handles.popup_segment,'Visible','on');    S = [];
    % Even if we have more than one fault, the segments popup will start with only the first fault's segments
    s_format = ['%.' num2str(fix(log10(nvert(1)))+1) 'd'];
    for (i=1:nvert(1))     S = [S; ['Segment ' num2str(i,s_format)]];   end
    set(handles.popup_segment,'String',cellstr(S))
else
    set(handles.popup_segment,'Visible','off')
    delete(findobj(hObject,'Style','text','Tag','fault_segment'))    % Otherwise it would reborn in Pro look
end

handles.fault_x = fault_x;
handles.fault_y = fault_y;
handles.nvert = nvert;
handles.hide_planes(1:handles.n_faults) = 0;
handles.dms_xinc = 0;           handles.dms_yinc = 0;
handles.one_or_zero = 1;        % For Grid Registration grids, which are the most common cases

handles.FaultLength = LineLength(handles.h_fault,handles.geog);
if (~fault_in)
	% Make them all cell arrays to simplify logic
	if (~iscell(handles.FaultLength))   handles.FaultLength = {handles.FaultLength};   end
	if (~iscell(handles.FaultStrike))   handles.FaultStrike = {handles.FaultStrike};   end
	if (~iscell(handles.fault_x))       handles.fault_x = {handles.fault_x};    handles.fault_y = {handles.fault_y};   end
	handles.DislocStrike = handles.FaultStrike;
	
	for (k=1:handles.n_faults)
		handles.FaultDip{k}(1:nvert(k)) = 45;       handles.FaultWidth{k}(1:nvert(k)) = NaN;
		handles.FaultDepth{k}(1:nvert(k)) = NaN;	handles.FaultTopDepth{k}(1:nvert(k)) = 0;
		handles.DislocSlip{k}(1:nvert(k)) = NaN;	handles.DislocRake{k}(1:nvert(k)) = NaN;
        handles.ux{k}(1:nvert(k)) = 1;              handles.uy{k}(1:nvert(k)) = 0;
        handles.uz{k}(1:nvert(k)) = 0;
	end
	
	z1 = num2str(handles.FaultLength{1}(1));    z2 = num2str(handles.FaultStrike{1}(1),'%.1f');
	z3 = num2str(handles.FaultDip{1}(1),'%.1f');
	set(handles.edit_FaultLength,'String',z1,'Enable','off')
	set(handles.edit_FaultStrike,'String',z2,'Enable','off')
	set(handles.edit_FaultDip,'String',z3)
	set(handles.edit_DislocStrike,'String',z2)
	set(handles.edit_DislocSlip,'String','')
	set(handles.edit_DislocRake,'String','')

	% Set a default unit dislocation as a thrust motion
	set(handles.edit_ux,'String',0)
	set(handles.edit_uy,'String',1)
	set(handles.edit_uz,'String',0)
	
	% Default the top depth fault to zero
	set(handles.edit_FaultTopDepth,'String','0')
else
    set(handles.edit_FaultLength,'String',num2str(handles.FaultLength{1}(1)));
end

% Set a default view vector (for interferograms s must be different)
sx = cos((90-handles.FaultStrike{1})*D2R);   handles.sx = sx;
sy = sin((90-handles.FaultStrike{1})*D2R);   handles.sy = sy;       handles.sz(1:nvert) = 0;
set(handles.edit_sx,'String',num2str(sx(1)))
set(handles.edit_sy,'String',num2str(sy(1)))
set(handles.edit_sz,'String','0')

%-----------
% Fill in the grid limits boxes (in case user wants to compute a grid)
% But also try to guess if we are dealing with other (m or km) than geogs
head = getappdata(handles.h_calling_fig,'GMThead');
handles.is_meters = 0;     handles.is_km = 0;
if (~isempty(head))
    if (~handles.geog)      % Try to guess if user units are km or meters
        dx = head(2) - head(1);   dy = head(4) - head(3);
        len = sqrt(dx.*dx + dy.*dy);         % Distance in user unites
        if (len > 1e5)      % If grid's diagonal > 1e5 consider we have meters
            handles.is_meters = 1;     handles.is_km = 0;
            set(handles.popup_GridCoords,'Value',2)
        else
            handles.is_meters = 0;     handles.is_km = 1;
            set(handles.popup_GridCoords,'Value',3)
        end
    end
else
    delete(hObject);    return
end
x1 = num2str(head(1),'%.10f');      x2 = num2str(head(2),'%.10f');
y1 = num2str(head(3),'%.10f');      y2 = num2str(head(4),'%.10f');
% But remove any possible trailing zeros
x1 = ddewhite(x1,'0');              x2 = ddewhite(x2,'0');
y1 = ddewhite(y1,'0');              y2 = ddewhite(y2,'0');
set(handles.edit_Xmin,'String',x1); set(handles.edit_Xmax,'String',x2)
set(handles.edit_Ymin,'String',y1); set(handles.edit_Ymax,'String',y2)
handles.x_min = head(1);            handles.x_max = head(2);
handles.y_min = head(3);            handles.y_max = head(4);

[m,n] = size(getappdata(handles.h_calling_fig,'dem_z'));

% Fill in the x,y_inc and nrow,ncol boxes
set(handles.edit_Nrows,'String',num2str(m))
set(handles.edit_Ncols,'String',num2str(n))
% Compute default xinc, yinc based on map limits
yinc = (head(4) - head(3)) / (m-1);   xinc = (head(2) - head(1)) / (n-1);
set(handles.edit_Yinc,'String',num2str(yinc,10))
set(handles.edit_Xinc,'String',num2str(xinc,10))
%-----------

handles.nrows = m;      handles.ncols = n;
handles.x_inc = xinc;   handles.y_inc = yinc;

%------------ Give a Pro look (3D) to the frame boxes  -------------------------------
bgcolor = get(0,'DefaultUicontrolBackgroundColor');
framecolor = max(min(0.65*bgcolor,[1 1 1]),[0 0 0]);
set(0,'Units','pixels');    set(hObject,'Units','pixels')    % Pixels are easier to reason with
h_f = findobj(hObject,'Style','Frame');
for i=1:length(h_f)
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
% Recopy the text fields on top of previously created frames (uistack is to slow)
h_t = findobj(hObject,'Style','Text');
for i=1:length(h_t)
    usr_d = get(h_t(i),'UserData');
    t_size = get(h_t(i),'Position');   t_str = get(h_t(i),'String');    fw = get(h_t(i),'FontWeight');
    bgc = get (h_t(i),'BackgroundColor');   fgc = get (h_t(i),'ForegroundColor');
    t_just = get(h_t(i),'HorizontalAlignment');     t_tag = get (h_t(i),'Tag');
    uicontrol('Parent',hObject, 'Style','text', 'Position',t_size,'String',t_str,'Tag',t_tag,...
        'BackgroundColor',bgc,'ForegroundColor',fgc,'FontWeight',fw,...
        'UserData',usr_d,'HorizontalAlignment',t_just);
end
delete(h_t)
%------------- END Pro look (3D) -------------------------------------------------------

% Choose default command line output for deform_okada_export
handles.output = hObject;
guidata(hObject, handles);

% UIWAIT makes deform_okada_export wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% ------------------------------------------------------------------------------------
set(hObject,'Visible','on');
% NOTE: If you make uiwait active you have also to uncomment the next three lines
% handles = guidata(hObject);
% out = deform_okada_OutputFcn(hObject, [], handles);
% varargout{1} = out;

% --- Outputs from this function are returned to the command line.
function varargout = deform_okada_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% Get default command line output from handles structure
varargout{1} = handles.output;

% ------------------------------------------------------------------------------------
function edit_FaultLength_Callback(hObject, eventdata, handles)
% Cannot be changed

% ------------------------------------------------------------------------------------
function edit_FaultWidth_Callback(hObject, eventdata, handles)
% Actualize the "FaultWidth" field
xx = str2double(get(hObject,'String'));
if (handles.n_faults > 1)   fault = get(handles.popup_fault,'Value');
else                        fault = 1;      end
if (xx < 0)         % If user tried to give a negative width
    xx = -xx;
    set(hObject,'String',num2str(xx))
end
dip = str2double(get(handles.edit_FaultDip,'String'));
top_d = str2double(get(handles.edit_FaultTopDepth,'String'));
depth = top_d + xx * cos((90-dip)*pi/180);
set(handles.edit_FaultDepth,'String',num2str(depth));
seg = get(handles.popup_segment,'Value');
handles.FaultWidth{fault}(seg) = xx;
handles.FaultDepth{fault}(seg) = depth;

% Update the patch that represents the surface projection of the fault plane
xx = [handles.fault_x{fault}(seg); handles.fault_x{fault}(seg+1)];
yy = [handles.fault_y{fault}(seg); handles.fault_y{fault}(seg+1)];

D2R = pi / 180;
off = handles.FaultWidth{fault}(seg) * cos(handles.FaultDip{fault}(seg)*D2R);
strk = handles.FaultStrike{fault}(seg);

if (handles.geog)
    rng = off / 6371 / D2R;
    [lat1,lon1] = circ_geo(yy(1),xx(1),rng,strk+90,1);
    [lat2,lon2] = circ_geo(yy(2),xx(2),rng,strk+90,1);
else
    if (handles.is_meters)  off = off * 1e3;    end
    lon1 = xx(1) + off * cos(strk*D2R);     lon2 = xx(2) + off * cos(strk*D2R);
    lat1 = yy(1) - off * sin(strk*D2R);     lat2 = yy(2) - off * sin(strk*D2R);
end
x = [xx(1) xx(2) lon2 lon1];    y = [yy(1) yy(2) lat2 lat1];
hp = getappdata(handles.h_fault(fault),'PatchHand');
try,    set(hp(seg),'XData',x,'YData',y,'FaceColor',[.8 .8 .8],'EdgeColor','k','LineWidth',1);  end
guidata(hObject, handles);

% ------------------------------------------------------------------------------------
function edit_FaultStrike_Callback(hObject, eventdata, handles)
% Cannot be changed

% ------------------------------------------------------------------------------------
function edit_FaultDip_Callback(hObject, eventdata, handles)
% Actualize the "FaultDip" field
xx = str2double(get(hObject,'String'));
if (handles.n_faults > 1)   fault = get(handles.popup_fault,'Value');
else                        fault = 1;      end
top_d = str2double(get(handles.edit_FaultTopDepth,'String'));
W = str2double(get(handles.edit_FaultWidth,'String'));
depth = top_d + W * cos((90-xx)*pi/180);
set(handles.edit_FaultDepth,'String',num2str(depth));
seg = get(handles.popup_segment,'Value');
handles.FaultDip{fault}(seg) = xx;
handles.FaultDepth{fault}(seg) = depth;

% Update the patch that represents the surface projection of the fault plane
xx = [handles.fault_x{fault}(seg); handles.fault_x{fault}(seg+1)];
yy = [handles.fault_y{fault}(seg); handles.fault_y{fault}(seg+1)];

D2R = pi / 180;
off = handles.FaultWidth{fault}(seg) * cos(handles.FaultDip{fault}(seg)*D2R);
strk = handles.FaultStrike{fault}(seg);

if (handles.geog)
    rng = off / 6371 / D2R;
    [lat1,lon1] = circ_geo(yy(1),xx(1),rng,strk+90,1);
    [lat2,lon2] = circ_geo(yy(2),xx(2),rng,strk+90,1);
else
    if (handles.is_meters)  off = off * 1e3;    end
    lon1 = xx(1) + off * cos(strk*D2R);     lon2 = xx(2) + off * cos(strk*D2R);
    lat1 = yy(1) - off * sin(strk*D2R);     lat2 = yy(2) - off * sin(strk*D2R);
end
x = [xx(1) xx(2) lon2 lon1];    y = [yy(1) yy(2) lat2 lat1];
hp = getappdata(handles.h_fault(fault),'PatchHand');
try,    set(hp(seg),'XData',x,'YData',y,'FaceColor',[.8 .8 .8],'EdgeColor','k','LineWidth',1);  end
guidata(hObject, handles);

% ------------------------------------------------------------------------------------
function edit_FaultDepth_Callback(hObject, eventdata, handles)
% Actualize the "FaultTopDepth" field
xx = str2double(get(hObject,'String'));
if (handles.n_faults > 1)   fault = get(handles.popup_fault,'Value');
else                        fault = 1;      end
if (xx < 0)         % If user tried to give a negative depth
    xx = -xx;
    set(hObject,'String',num2str(xx))
end
W = str2double(get(handles.edit_FaultWidth,'String'));
dip = str2double(get(handles.edit_FaultDip,'String'));
top_d = xx - W * cos((90-dip)*pi/180);
set(handles.edit_FaultTopDepth,'String',num2str(top_d));
seg = get(handles.popup_segment,'Value');
handles.FaultDepth{fault}(seg) = xx;
handles.FaultTopDepth{fault}(seg) = top_d;
guidata(hObject, handles);

% ------------------------------------------------------------------------------------
function edit_FaultTopDepth_Callback(hObject, eventdata, handles)
% Actualize the "FaultDepth" field
xx = str2double(get(hObject,'String'));
if (handles.n_faults > 1)   fault = get(handles.popup_fault,'Value');
else                        fault = 1;      end
if (xx < 0)         % If user tried to give a negative depth
    xx = -xx;
    set(hObject,'String',num2str(xx))
end
W = str2double(get(handles.edit_FaultWidth,'String'));
dip = str2double(get(handles.edit_FaultDip,'String'));
depth = xx + W * cos((90-dip)*pi/180);
set(handles.edit_FaultDepth,'String',num2str(depth));
seg = get(handles.popup_segment,'Value');
handles.FaultTopDepth{fault}(seg) = xx;
handles.FaultDepth{fault}(seg) = depth;
guidata(hObject, handles);

% ------------------------------------------------------------------------------------
function popup_segment_Callback(hObject, eventdata, handles)
seg = get(hObject,'Value');
if (handles.n_faults > 1)   fault = get(handles.popup_fault,'Value');
else                        fault = 1;      end

% Fault parameters
set(handles.edit_FaultLength,'String',num2str(handles.FaultLength{fault}(seg)))
set(handles.edit_FaultStrike,'String',num2str(handles.FaultStrike{fault}(seg),'%.1f'))

if (isnan(handles.FaultWidth{fault}(seg)))     str = '';
else    str = num2str(handles.FaultWidth{fault}(seg));     end
set(handles.edit_FaultWidth,'String',str)

set(handles.edit_FaultDip,'String',num2str(handles.FaultDip{fault}(seg),'%.1f'))
set(handles.edit_FaultTopDepth,'String',num2str(handles.FaultTopDepth{fault}(seg)))

if (isnan(handles.FaultDepth{fault}(seg)))     str = '';
else    str = num2str(handles.FaultDepth{fault}(seg));     end
set(handles.edit_FaultDepth,'String',str)

% Dislocation parameters
set(handles.edit_DislocStrike,'String',num2str(handles.DislocStrike{fault}(seg),'%.1f'))
if (isnan(handles.DislocSlip{fault}(seg)))     str = '';
else    str = num2str(handles.DislocSlip{fault}(seg));     end
set(handles.edit_DislocSlip,'String',str)
if (isnan(handles.DislocRake{fault}(seg)))     str = '';
else    str = num2str(handles.DislocRake{fault}(seg),'%.1f');     end
set(handles.edit_DislocRake,'String',str)

set(handles.edit_ux,'String',num2str(handles.ux{fault}(seg)))
set(handles.edit_uy,'String',num2str(handles.uy{fault}(seg)))
set(handles.edit_uz,'String',num2str(handles.uz{fault}(seg)))

% -----------------------------------------------------------------------------------------
function popup_fault_Callback(hObject, eventdata, handles)
fault = get(hObject,'Value');
S = [];
s_format = ['%.' num2str(fix(log10(handles.nvert(fault)))+1) 'd'];
for (i=1:handles.nvert(fault))     S = [S; ['Segment ' num2str(i,s_format)]];   end
set(handles.popup_segment,'String',cellstr(S),'Value',1)    
seg = 1;    % Make current the first segment

% Identify the currently active fault by setting its linestyle to dash
set(handles.h_fault,'LineStyle','-')
set(handles.h_fault(fault),'LineStyle','--')

% Set the hide planes checkbox with the correct value for this fault
if (handles.hide_planes(fault))
    set(handles.checkbox_hideFaultPlanes,'Value',1)
else
    set(handles.checkbox_hideFaultPlanes,'Value',0)
end

% Fault parameters
set(handles.edit_FaultLength,'String',num2str(handles.FaultLength{fault}(seg)))
set(handles.edit_FaultStrike,'String',num2str(handles.FaultStrike{fault}(seg),'%.1f'))

if (isnan(handles.FaultWidth{fault}(seg)))     str = '';
else    str = num2str(handles.FaultWidth{fault}(seg));     end
set(handles.edit_FaultWidth,'String',str)

set(handles.edit_FaultDip,'String',num2str(handles.FaultDip{fault}(seg),'%.1f'))
set(handles.edit_FaultTopDepth,'String',num2str(handles.FaultTopDepth{fault}(seg)))

if (isnan(handles.FaultDepth{fault}(seg)))     str = '';
else    str = num2str(handles.FaultDepth{fault}(seg));     end
set(handles.edit_FaultDepth,'String',str)

% Dislocation parameters
set(handles.edit_DislocStrike,'String',num2str(handles.DislocStrike{fault}(seg),'%.1f'))
if (isnan(handles.DislocSlip{fault}(seg)))     str = '';
else    str = num2str(handles.DislocSlip{fault}(seg));     end
set(handles.edit_DislocSlip,'String',str)
if (isnan(handles.DislocRake{fault}(seg)))     str = '';
else    str = num2str(handles.DislocRake{fault}(seg),'%.1f');     end
set(handles.edit_DislocRake,'String',str)
if (isnan(handles.ux{fault}(seg)))     str = '';
else    str = num2str(handles.ux{fault}(seg),'%.1f');     end
set(handles.edit_ux,'String',str)
if (isnan(handles.uy{fault}(seg)))     str = '';
else    str = num2str(handles.uy{fault}(seg),'%.1f');     end
set(handles.edit_uy,'String',str)
if (isnan(handles.uz{fault}(seg)))     str = '0';
else    str = num2str(handles.uz{fault}(seg),'%.1f');     end
set(handles.edit_uz,'String',str)
refresh(handles.h_calling_fig);         % otherwise, ML BUG

% ---------------------------------------------------------------
function popup_GridCoords_Callback(hObject, eventdata, handles)
xx = get(hObject,'Value');
if (xx == 1)        handles.geog = 1;       handles.is_meters = 0;  handles.is_km = 0;
elseif (xx == 2)    handles.is_meters = 1;  handles.is_geog = 0;    handles.is_km = 0;
elseif (xx == 3)    handles.is_km = 1;      handles.is_geog = 0;    handles.is_meters = 0;
end
guidata(hObject,handles)

% ------------------------------------------------------------------------------------
function edit_DislocStrike_Callback(hObject, eventdata, handles)
D2R = pi / 180;
if (handles.n_faults > 1)   fault = get(handles.popup_fault,'Value');
else                        fault = 1;      end
xx = str2double(get(hObject,'String'));
f_strike = str2double(get(handles.edit_FaultStrike,'String'));      % Fault strike
slip = str2double(get(handles.edit_DislocSlip,'String'));           % Dislocation slip
ux = slip * cos((f_strike-xx)*D2R);
uz = slip * sin((f_strike-xx)*D2R);
set(handles.edit_ux,'String',num2str(ux))
set(handles.edit_uz,'String',num2str(uz))
seg = get(handles.popup_segment,'Value');
handles.DislocStrike{fault}(seg) = xx;
handles.ux{fault}(seg) = ux;   handles.uz{fault}(seg) = uz;
guidata(hObject, handles);

% ------------------------------------------------------------------------------------
function edit_DislocRake_Callback(hObject, eventdata, handles)
D2R = pi / 180;
xx = str2double(get(hObject,'String'));
if (handles.n_faults > 1)   fault = get(handles.popup_fault,'Value');
else                        fault = 1;      end
f_strike = str2double(get(handles.edit_FaultStrike,'String'));      % Fault strike
d_strike = str2double(get(handles.edit_DislocStrike,'String'));     % Dislocation strike
slip = str2double(get(handles.edit_DislocSlip,'String'));           % Dislocation slip
ux = slip * cos(xx*D2R) * cos((f_strike - d_strike)*D2R);
uy = slip * sin(xx*D2R) * cos((f_strike - d_strike)*D2R);
set(handles.edit_ux,'String',num2str(ux))
set(handles.edit_uy,'String',num2str(uy))
seg = get(handles.popup_segment,'Value');
handles.DislocRake{fault}(seg) = xx;
handles.ux{fault}(seg) = ux;   handles.uy{fault}(seg) = uy;
guidata(hObject, handles);

% ------------------------------------------------------------------------------------
function edit_DislocSlip_Callback(hObject, eventdata, handles)
xx = str2double(get(hObject,'String'));
if (handles.n_faults > 1)   fault = get(handles.popup_fault,'Value');
else                        fault = 1;      end
seg = get(handles.popup_segment,'Value');
if (isnan(xx))
    set(hObject,'String','')
    handles.DislocSlip{fault}(seg) = NaN;
else
    handles.DislocSlip{fault}(seg) = xx;
end
guidata(hObject, handles);

% ------------------------------------------------------------------------------------
function edit_ux_Callback(hObject, eventdata, handles)
if ( isempty(get(hObject,'String')) ),  set(hObject,'String','0');   return;     end
if (handles.n_faults > 1)   fault = get(handles.popup_fault,'Value');
else                        fault = 1;      end
ux = str2double(get(hObject,'String'));
uy = str2double(get(handles.edit_uy,'String'));
uz = str2double(get(handles.edit_uz,'String'));
slip = sqrt(ux^2 + uy^2 + uz^2);
set(handles.edit_DislocSlip,'String',num2str(slip))
seg = get(handles.popup_segment,'Value');
handles.ux{fault}(seg) = ux;
handles.DislocSlip{fault}(seg) = slip;
guidata(hObject, handles);

% ------------------------------------------------------------------------------------
function edit_uy_Callback(hObject, eventdata, handles)
if ( isempty(get(hObject,'String')) ),  set(hObject,'String','0');   return;     end
if (handles.n_faults > 1)   fault = get(handles.popup_fault,'Value');
else                        fault = 1;      end
uy = str2double(get(hObject,'String'));
ux = str2double(get(handles.edit_ux,'String'));
uz = str2double(get(handles.edit_uz,'String'));
slip = sqrt(ux^2 + uy^2 + uz^2);
set(handles.edit_DislocSlip,'String',num2str(slip))
seg = get(handles.popup_segment,'Value');
handles.uy{fault}(seg) = uy;
handles.DislocSlip{fault}(seg) = slip;
guidata(hObject, handles);

% ------------------------------------------------------------------------------------
function edit_uz_Callback(hObject, eventdata, handles)
if ( isempty(get(hObject,'String')) ),  set(hObject,'String','0');   return;     end
if (handles.n_faults > 1)   fault = get(handles.popup_fault,'Value');
else                        fault = 1;      end
uz = str2double(get(hObject,'String'));
ux = str2double(get(handles.edit_ux,'String'));
uy = str2double(get(handles.edit_uy,'String'));
slip = sqrt(ux^2 + uy^2 + uz^2);
set(handles.edit_DislocSlip,'String',num2str(slip))
seg = get(handles.popup_segment,'Value');
handles.uz{fault}(seg) = uz;
handles.DislocSlip{fault}(seg) = slip;
guidata(hObject, handles);

% ------------------------------------------------------------------------------------
function edit_Xmin_Callback(hObject, eventdata, handles)
xx = get(hObject,'String');     val = test_dms(xx);
if ~isempty(val)            % when dd:mm or dd:mm:ss was given
    x_min = 0;
    if str2double(val{1}) > 0
        for i = 1:length(val)   x_min = x_min + str2double(val{i}) / (60^(i-1));    end
    else
        for i = 1:length(val)   x_min = x_min - abs(str2double(val{i})) / (60^(i-1));   end
    end
    handles.x_min = x_min;
    if ~isempty(handles.x_max) & x_min >= handles.x_max
        errordlg('West Longitude >= East Longitude ','Error in Longitude limits')
        set(hObject,'String','');   guidata(hObject, handles);  return
    end
    nc = get(handles.edit_Ncols,'String');
    if ~isempty(handles.x_max) & ~isempty(nc)       % x_max and ncols boxes are filled
        % Compute Ncols, but first must recompute x_inc
        x_inc = ivan_the_terrible((handles.x_max - x_min),round(abs(str2double(nc))),1);
        xx = floor((handles.x_max - str2double(xx)) / (str2double(get(handles.edit_Xinc,'String')))+0.5) + handles.one_or_zero;
        set(handles.edit_Xinc,'String',num2str(x_inc,10))
    elseif ~isempty(handles.x_max)      % x_max box is filled but ncol is not, so put to the default (100)
        x_inc = ivan_the_terrible((handles.x_max - x_min),100,1);
        set(handles.edit_Xinc,'String',num2str(x_inc,10))
        set(handles.edit_Ncols,'String','100')
    end
else                % box is empty, so clear also x_inc and ncols
    set(handles.edit_Xinc,'String','');     set(handles.edit_Ncols,'String','');
    set(hObject,'String','');
end
guidata(hObject, handles);

% ------------------------------------------------------------------------------------
function edit_Xmax_Callback(hObject, eventdata, handles)
xx = get(hObject,'String');     val = test_dms(xx);
if ~isempty(val)
    x_max = 0;
    if str2double(val{1}) > 0
        for i = 1:length(val)   x_max = x_max + str2double(val{i}) / (60^(i-1));    end
    else
        for i = 1:length(val)   x_max = x_max - abs(str2double(val{i})) / (60^(i-1));   end
    end
    handles.x_max = x_max;
    if ~isempty(handles.x_min) & x_max <= handles.x_min 
        errordlg('East Longitude <= West Longitude','Error in Longitude limits')
        set(hObject,'String','');   guidata(hObject, handles);  return
    end
    nc = get(handles.edit_Ncols,'String');
    if ~isempty(handles.x_min) & ~isempty(nc)       % x_max and ncols boxes are filled
        % Compute Ncols, but first must recompute x_inc
        x_inc = ivan_the_terrible((x_max - handles.x_min),round(abs(str2double(nc))),1);
        xx = floor((handles.x_min - str2double(xx)) / (str2double(get(handles.edit_Xinc,'String')))+0.5) + handles.one_or_zero;
        set(handles.edit_Xinc,'String',num2str(x_inc,10))
    elseif ~isempty(handles.x_min)      % x_min box is filled but ncol is not, so put to the default (100)
        x_inc = ivan_the_terrible((x_max - handles.x_min),100,1);
        set(handles.edit_Xinc,'String',num2str(x_inc,10))
        set(handles.edit_Ncols,'String','100')
    end
else                % box is empty, so clear also x_inc and ncols
    set(handles.edit_Xinc,'String','');     set(handles.edit_Ncols,'String','');
    set(hObject,'String','');
end
guidata(hObject, handles);

% ------------------------------------------------------------------------------------
function edit_Xinc_Callback(hObject, eventdata, handles)
dms = 0;
xx = get(hObject,'String');     val = test_dms(xx);
if isempty(val),    return;     end
% If it survived then ...
if length(val) > 1    dms = 1;      end         % inc given in dd:mm or dd:mm:ss format
x_inc = 0;
for i = 1:length(val)   x_inc = x_inc + str2double(val{i}) / (60^(i-1));    end
if ~isempty(handles.x_min) & ~isempty(handles.x_max)
    % Make whatever x_inc given compatible with GMT_grd_RI_verify
    x_inc = ivan_the_terrible((handles.x_max - handles.x_min), x_inc,2);
    if ~dms         % case of decimal unities
        set(hObject,'String',num2str(x_inc,8))
        ncol = floor((handles.x_max - handles.x_min) / x_inc + 0.5) + handles.one_or_zero;
    else            % inc was in dd:mm or dd:mm:ss format
        ncol = floor((handles.x_max - handles.x_min) / x_inc + 0.5) + handles.one_or_zero;
        ddmm = dec2deg(x_inc);
        set(hObject,'String',ddmm)
    end
    set(handles.edit_Ncols,'String',num2str(ncol))
end
handles.dms_xinc = dms;     handles.x_inc = str2double(xx);
guidata(hObject, handles);

% ------------------------------------------------------------------------------------
function edit_Ncols_Callback(hObject, eventdata, handles)
xx = get(hObject,'String');
if isempty(xx)          % Idiot user attempt. Reset ncols.
    set(hObject,'String',handles.ncols);    return;
end
if ~isempty(get(handles.edit_Xmin,'String')) & ~isempty(get(handles.edit_Xmax,'String')) & ...
        ~isempty(get(handles.edit_Xinc,'String')) & ~isempty(xx)
    x_inc = ivan_the_terrible((handles.x_max - handles.x_min),round(abs(str2double(xx))),1);
    if handles.dms_xinc        % x_inc was given in dd:mm:ss format
        ddmm = dec2deg(x_inc);
        set(handles.edit_Xinc,'String',ddmm)
    else                    % x_inc was given in decimal format
        set(handles.edit_Xinc,'String',num2str(x_inc,10));
    end
    handles.ncols = str2double(xx);
end
guidata(hObject, handles);

% ------------------------------------------------------------------------------------
function edit_Ymin_Callback(hObject, eventdata, handles)
% Read value either in decimal or in the dd:mm or dd_mm:ss formats and do some tests
xx = get(hObject,'String');     val = test_dms(xx);
if ~isempty(val)
    y_min = 0;
    if str2double(val{1}) > 0
        for i = 1:length(val)   y_min = y_min + str2double(val{i}) / (60^(i-1));    end
    else
        for i = 1:length(val)   y_min = y_min - abs(str2double(val{i})) / (60^(i-1));   end
    end
    handles.y_min = y_min;
    if ~isempty(handles.y_max) & y_min >= handles.y_max
        errordlg('South Latitude >= North Latitude','Error in Latitude limits')
        set(hObject,'String','');   guidata(hObject, handles);  return
    end
    nr = get(handles.edit_Nrows,'String');
    if ~isempty(handles.y_max) & ~isempty(nr)       % y_max and nrows boxes are filled
        % Compute Nrowss, but first must recompute y_inc
        y_inc = ivan_the_terrible((handles.y_max - y_min),round(abs(str2double(nr))),1);
        xx = floor((handles.y_max - str2double(xx)) / (str2double(get(handles.edit_Yinc,'String')))+0.5) + handles.one_or_zero;
        set(handles.edit_Yinc,'String',num2str(y_inc,10))
    elseif ~isempty(handles.y_max)      % y_max box is filled but nrows is not, so put to the default (100)
        y_inc = ivan_the_terrible((handles.y_max - y_min),100,1);
        set(handles.edit_Yinc,'String',num2str(y_inc,10))
        set(handles.edit_Nrows,'String','100')
    end
else                % box is empty, so clear also y_inc and nrows
    set(handles.edit_Yinc,'String','');     set(handles.edit_Nrows,'String','');
    set(hObject,'String','');
end
guidata(hObject, handles);

% ------------------------------------------------------------------------------------
function edit_Ymax_Callback(hObject, eventdata, handles)
% Read value either in decimal or in the dd:mm or dd_mm:ss formats and do some tests
xx = get(hObject,'String');     val = test_dms(xx);
if ~isempty(val)
    y_max = 0;
    if str2double(val{1}) > 0
        for i = 1:length(val)   y_max = y_max + str2double(val{i}) / (60^(i-1));    end
    else
        for i = 1:length(val)   y_max = y_max - abs(str2double(val{i})) / (60^(i-1));   end
    end
    handles.y_max = y_max;
    if ~isempty(handles.y_min) & y_max <= handles.y_min 
        errordlg('North Latitude <= South Latitude','Error in Latitude limits')
        set(hObject,'String','');   guidata(hObject, handles);  return
    end
    nr = get(handles.edit_Nrows,'String');
    if ~isempty(handles.y_min) & ~isempty(nr)       % y_min and nrows boxes are filled
        % Compute Nrows, but first must recompute y_inc
        y_inc = ivan_the_terrible((y_max - handles.y_min),round(abs(str2double(nr))),1);
        xx = floor((handles.y_min - str2double(xx)) / (str2double(get(handles.edit_Yinc,'String')))+0.5) + handles.one_or_zero;
        set(handles.edit_Yinc,'String',num2str(y_inc,10))
    elseif ~isempty(handles.y_min)      % y_min box is filled but nrows is not, so put to the default (100)
        y_inc = ivan_the_terrible((y_max - handles.y_min),100,1);
        set(handles.edit_Yinc,'String',num2str(y_inc,10))
        set(handles.edit_Nrows,'String','100')
    end
else                % This box is empty, so clear also y_inc and nrows
    set(handles.edit_Yinc,'String','');     set(handles.edit_Nrows,'String','');
    set(hObject,'String','');
end
guidata(hObject, handles);

% ------------------------------------------------------------------------------------
function edit_Yinc_Callback(hObject, eventdata, handles)
dms = 0;
xx = get(hObject,'String');     val = test_dms(xx);
if isempty(val)
    set(hObject, 'String', '');    return
end
% If it survived then ...
if length(val) > 1    dms = 1;      end         % inc given in dd:mm or dd:mm:ss format
y_inc = 0;
for i = 1:length(val)   y_inc = y_inc + str2double(val{i}) / (60^(i-1));    end
if ~isempty(handles.y_min) & ~isempty(handles.y_max)
    % Make whatever y_inc given compatible with GMT_grd_RI_verify
    y_inc = ivan_the_terrible((handles.y_max - handles.y_min), y_inc,2);
    if ~dms         % case of decimal unities
        set(hObject,'String',num2str(y_inc,10))
        nrow = floor((handles.y_max - handles.y_min) / y_inc + 0.5) + handles.one_or_zero;
    else            % inc was in dd:mm or dd:mm:ss format
        nrow = floor((handles.y_max - handles.y_min) / y_inc + 0.5) + handles.one_or_zero;
        ddmm = dec2deg(y_inc);
        set(hObject,'String',ddmm)
    end
    set(handles.edit_Nrows,'String',num2str(nrow))
end
handles.dms_yinc = dms;     handles.y_inc = str2double(xx);
guidata(hObject, handles);

% ------------------------------------------------------------------------------------
function edit_Nrows_Callback(hObject, eventdata, handles)
xx = get(hObject,'String');
if isempty(xx)          % Idiot user attempt. Reset nrows.
    set(hObject,'String',handles.nrows);    return;
end
if ~isempty(get(handles.edit_Ymin,'String')) & ~isempty(get(handles.edit_Ymax,'String')) & ...
        ~isempty(get(handles.edit_Yinc,'String'))
    y_inc = ivan_the_terrible((handles.y_max - handles.y_min),round(abs(str2double(xx))),1);
    if handles.dms_yinc        % y_inc was given in dd:mm:ss format
        ddmm = dec2deg(y_inc);
        set(handles.edit_Yinc,'String',ddmm)
    else                    % y_inc was given in decimal format
        set(handles.edit_Yinc,'String',num2str(y_inc,10));
    end
    handles.nrows = str2double(xx);
end
guidata(hObject, handles);

% ------------------------------------------------------------------------------------
function checkbox_ToggleXY_Callback(hObject, eventdata, handles)

% ------------------------------------------------------------------------------------
function pushbutton_Help_R_Callback(hObject, eventdata, handles)
message = {'That''s prety obvious to guess what this option does. You select an area,'
    'the grid spacing or the number of rows/columns and the deformation will'
    'be computed at all nodes of that grid.'};
helpdlg(message,'Help on deformation grid');

% ------------------------------------------------------------------------------------
function checkbox_Option_H_Callback(hObject, eventdata, handles)

% ------------------------------------------------------------------------------------
function edit_nHeaders_Callback(hObject, eventdata, handles)

% ------------------------------------------------------------------------------------
function pushbutton_Help_H_Callback(hObject, eventdata, handles)
message = {'If you have a file with x,y positions, then the deformation will be computed at those postions'}
helpdlg(message,'Little Help');

% ------------------------------------------------------------------------------------
function edit_InputFile_Callback(hObject, eventdata, handles)
fname = get(hObject,'String');
if isempty(fname),  handles.input_locations = [];   return;   end
hFig = gcf;
[bin,n_column,multi_seg,n_headers] = guess_file(fname);
% If msgbox exist we have to move it from behind the main window. So get it's handle
hMsgFig = gcf;
if (hFig ~= hMsgFig)        uistack(hMsgFig,'top');   end   % If msgbox exists, bring it forward
% If error in reading file
if isempty(bin) & isempty(n_column) & isempty(multi_seg) & isempty(n_headers)
    errordlg(['Error reading file ' fname],'Error');    return
end
if multi_seg ~= 0   % multisegments are not spported
    errordlg('Multisegment files are yet not supported.','Error');   return
end
if (bin == 0)   % ASCII
    if n_column < 2
        errordlg('File error. Your file doesn''t have at least 2 columns','Error'); return
    end
    handles.input_locations = read_xy(fname,n_column,n_headers);
    if (hFig ~= hMsgFig);       figure(hFig);   end     % gain access to the drawing figure
    [nr,nc] = size(handles.input_locations);
    if (nr == 0)
        errordlg('Your file is empty.','Chico Clever');   return
    end
    if (n_headers > 0)      % We have headers in file (ai!, ai!)
        set(handles.checkbox_Option_H,'Value',1)
        set(handles.edit_nHeaders,'String',num2str(n_headers))
    end
else        % BINARY
    errordlg('Sorry, reading binary files is not yet programed','Error');   return
end
guidata(hObject,handles)

% ------------------------------------------------------------------------------------
function pushbutton_InputFile_Callback(hObject, eventdata, handles)
if (~isempty(handles.h_calling_fig))                    % If we know the handle to the calling fig
    cfig_handles = guidata(handles.h_calling_fig);      % get handles of the calling fig
    last_dir = cfig_handles.last_dir;
    home = cfig_handles.home_dir;
else
    last_dir = [];
end

if (~isempty(last_dir)),    cd(last_dir);   end
[FileName,PathName] = uigetfile({'*.dat;*.DAT;*.xy', 'Maregraph location (*.dat,*.DAT,*.xy)';'*.*', 'All Files (*.*)'},'Select Maregraphs position');
if (~isempty(last_dir)),    cd(home);   end
if isequal(FileName,0);     return;     end
fname = [PathName FileName];
hFig = gcf;
[bin,n_column,multi_seg,n_headers] = guess_file(fname);
% If msgbox exist we have to move it from behind the main window. So get it's handle
hMsgFig = gcf;
if (hFig ~= hMsgFig)        uistack(hMsgFig,'top');   end   % If msgbox exists, bring it forward
% If error in reading file
if isempty(bin) & isempty(n_column) & isempty(multi_seg) & isempty(n_headers)
    errordlg(['Error reading file ' fname],'Error');    return
end
if multi_seg ~= 0   % multisegments are not spported
    errordlg('Multisegment files are yet not supported.','Error');   return
end
if (bin == 0)   % ASCII
    if (n_column < 2)
        errordlg('File error. Your file doesn''t have at least 2 columns','Error'); return
    end
    handles.input_locations = read_xy(fname,n_column,n_headers);
    if (hFig ~= hMsgFig);       figure(hFig);   end     % gain access to the drawing figure
    [nr,nc] = size(handles.input_locations);
    if (nr == 0)
        errordlg('Your file is empty.','Chico Clever');   return
    end
    if (n_headers > 0)      % We have headers in file (ai!, ai!)
        set(handles.checkbox_Option_H,'Value',1)
        set(handles.edit_nHeaders,'String',num2str(n_headers))
    end
else        % BINARY
    errordlg('Sorry, reading binary files is not yet programed','Error');   return
end
set(handles.edit_InputFile,'String',fname)
guidata(hObject,handles)

% ------------------------------------------------------------------------------------
function edit_sx_Callback(hObject, eventdata, handles)
if ( isempty(get(hObject,'String')) ),  set(hObject,'String','0');   return;     end

% ------------------------------------------------------------------------------------
function edit_sy_Callback(hObject, eventdata, handles)
if ( isempty(get(hObject,'String')) ),  set(hObject,'String','0');   return;     end

% ------------------------------------------------------------------------------------
function edit_sz_Callback(hObject, eventdata, handles)
if ( isempty(get(hObject,'String')) ),  set(hObject,'String','0');   return;     end

% ------------------------------------------------------------------------------------
function pushbutton_compute_Callback(hObject, eventdata, handles)
% If cartesian coordinates, they must be in meters
if (any(isnan(cat(2,handles.FaultWidth{:}))))
    errordlg('One or more segments where not set with the fault''s Width','Error');    return
end
if (any(isnan(cat(2,handles.FaultDepth{:}))))
    errordlg('One or more segments where not set with the fault''s Depth','Error');    return
end
if (any(isnan(cat(2,handles.DislocSlip{:}))))
    errordlg('One or more segments where not set with the movement''s slip','Error');    return
end

if( all(cat(2,handles.ux{:}) == 0) & all(cat(2,handles.uy{:}) == 0) & all(cat(2,handles.uz{:}) == 0) )
    errordlg('No movement along faults(s), nothing to compute there.','Error')
    return
end

sx = str2double(get(handles.edit_sx,'String'));
sy = str2double(get(handles.edit_sy,'String'));
sz = str2double(get(handles.edit_sz,'String'));
s = [sx sy sz];

if( sx == 0 & sy == 0 & sz == 0 )
    errordlg('Looking vector is looking nowhere.','Error')
    return
elseif (sqrt(sx^2 + sy^2 + sz^2) < 0.99)
    warndlg('Norm of looking vector is less than 1. It means you are loosing deformation.','Warning')
end

% Get grid params
xmin = str2double(get(handles.edit_Xmin,'String'));     xmax = str2double(get(handles.edit_Xmax,'String'));
ymin = str2double(get(handles.edit_Ymin,'String'));     ymax = str2double(get(handles.edit_Ymax,'String'));
nrow = str2double(get(handles.edit_Nrows,'String'));    ncol = str2double(get(handles.edit_Ncols,'String'));

x = handles.fault_x;    y = handles.fault_y;
if (~iscell(x))         x = {x};    y = {y};    end
fig_xlim = [xmin xmax];   fig_ylim = [ymin ymax];
to_km = 1;      % The conversion from m->km will be done inside range_change

if (handles.geog)    
    opt_M = '-M';
else
    if (handles.is_meters)      to_km = 1000;   end
    for (i=1:handles.n_faults)
        x{i} = x{i} / to_km;   y{i} = y{i} / to_km;
        handles.FaultLength{i} = handles.FaultLength{i} / to_km;
    end
    opt_M = '';
end

for (i=1:handles.n_faults)
    % I have to do fish the patch coords because range_change seams to use not
    % the fault trace coords but the coordinates of the fault at its depth 
    hp = getappdata(handles.h_fault(i),'PatchHand');
    xp = get(hp,'XData');    yp = get(hp,'YData');
    if (iscell(xp))
        x{i} = [];    y{i} = [];
        for (k=1:length(xp))
            x{i} = [x{i}; xp{k}(4)/to_km];
            y{i} = [y{i}; yp{k}(4)/to_km];
        end
    else
        x{i} = [xp(4)]/to_km;   y{i} = [yp(4)]/to_km;
    end
end

if (isempty(handles.input_locations))   % If ground positions were not given, compute a grid
    E = linspace(fig_xlim(1),fig_xlim(2),ncol)/to_km;
    N = linspace(fig_ylim(1),fig_ylim(2),nrow)/to_km;
    N = N(:);               % From the rngchn example, y coords are in a column vector

    % Compute deformation
    %U = range_change(x,y,strike,depth,dip,ux,uy,uz,L,W,E,N,s);
    if (handles.n_faults > 1)           % We have multiple faults
        U = zeros(nrow,ncol);
        h = waitbar(0,'Computing deformation');
        for (i=1:handles.n_faults)
            waitbar(i/handles.n_faults)
            U0 = range_change(x{i}(:),y{i}(:),handles.FaultStrike{i}(:),handles.FaultDepth{i}(:),handles.FaultDip{i}(:),...
                handles.ux{i}(:),handles.uy{i}(:),handles.uz{i}(:),handles.FaultLength{i}(:),handles.FaultWidth{i}(:),...
                E,N,s,opt_M);
            U = U0 + U;
        end
        close(h);    clear U0;
    else                                % We have only one fault
        U = range_change(x{1}(:),y{1}(:),handles.FaultStrike{1}(:),handles.FaultDepth{1}(:),handles.FaultDip{1}(:),...
            handles.ux{1}(:),handles.uy{1}(:),handles.uz{1}(:),handles.FaultLength{1}(:),handles.FaultWidth{1}(:),...
            E,N,s,opt_M);
    end

    z_max = max(U(:));     z_min = min(U(:));
    dx = str2double(get(handles.edit_Xinc,'String'));
    dy = str2double(get(handles.edit_Yinc,'String'));

    if (handles.geog)
        tmp = [xmin xmax ymin ymax z_min z_max 0];
        head.head = [tmp dx dy];
        head.X = linspace(xmin,xmax,ncol);
        tmp = linspace(ymin,ymax,nrow);
        head.Y = tmp';
    else
        E = E * to_km;   N = N * to_km;     % Convert to grid coords
        head.head = [E(1) E(end) N(1) N(end) z_min z_max 0 E(2)-E(1) N(2)-N(1)];
        head.X = E;     head.Y = N';
    end

    if get(handles.radiobutton_deformation,'Value')
        new_window = mirone(U,head,'Deformation',handles.h_calling_fig);
    else
        new_window = mirone(U,head,'Interfero',28.4);
    end
else        % Ground positions were given
    E = handles.input_locations(:,1)/to_km;
    N = handles.input_locations(:,2)/to_km;
    if (handles.n_faults > 1)           % We have multiple faults
        U0 = [];
        for (i=1:handles.n_faults)
            if (get(handles.checkbox_ToggleXY,'Value'))
                U0 = range_change(x{i}(:),y(:),handles.FaultStrike(:),handles.FaultDepth(:),handles.FaultDip(:),...
                    handles.ux(:),handles.uy(:),handles.uz(:),handles.FaultLength(:),handles.FaultWidth(:),...
                    N,E,s,opt_M);
            else
                U0 = range_change(x(:),y(:),handles.FaultStrike(:),handles.FaultDepth(:),handles.FaultDip(:),...
                    handles.ux(:),handles.uy(:),handles.uz(:),handles.FaultLength(:),handles.FaultWidth(:),...
                    E,N,s,opt_M);
            end
            U = U0 + U;
        end
        clear U0;
    else                                % We have only one fault
        if (get(handles.checkbox_ToggleXY,'Value'))
            U = range_change(x{1}(:),y{1}(:),handles.FaultStrike{1}(:),handles.FaultDepth{1}(:),handles.FaultDip{1}(:),...
                handles.ux(:),handles.uy{1}(:),handles.uz{1}(:),handles.FaultLength{1}(:),handles.FaultWidth{1}(:),...
                N,E,s,opt_M);
        else
            U = range_change(x{1}(:),y{1}(:),handles.FaultStrike{1}(:),handles.FaultDepth{1}(:),handles.FaultDip{1}(:),...
                handles.ux(:),handles.uy{1}(:),handles.uz{1}(:),handles.FaultLength{1}(:),handles.FaultWidth{1}(:),...
                E,N,s,opt_M);
        end
    end
    [FileName,PathName] = uiputfile({'*.dat;*.DAT', 'Data files (*.dat,*.DAT)';'*.*', 'All Files (*.*)'},'Select output deformation file');
    if isequal(FileName,0);     return;      end     % User gave up
    double2ascii([PathName FileName],[U(:,1) U(:,2) U(:,3)],'%f\t%f\t%f');     %NAO SEI SE E ASSIM. FALTA TESTAR
end

% ------------------------------------------------------------------------------------
function [lonlim,zone] = utmorigin(lon)
% Returns the UTM longitude limits and the Zone. Note that there is no way of telling
% if it is a North or South zone, but that should be easy by knowing the latitude.
lons = [-180:6:180]';	
ind = find(lons <= lon);  lonsidx = ind(max(ind));
		
if (lonsidx < 1 | lonsidx > 61)     lonsidx = [];
elseif (lonsidx == 61)              lonsidx = 60;    end

zone = [num2str(lonsidx)];

if (length(zone) == 1)
    lonsidx = str2num(zone);
elseif (length(zone) == 2)
    num = str2num(zone);
    if isempty(num)
        lonsidx = str2num(zone(1));
    else
        lonsidx = num;
    end
end
latlim = [];  lonlim = [];
lonlims = [(-180:6:174)' (-174:6:180)'];
lonlim = lonlims(lonsidx,:);

% ------------------------------------------------------------------------------------
function radiobutton_deformation_Callback(hObject, eventdata, handles)
if get(hObject,'Value')
    set(handles.radiobutton_interfero,'Value',0)
    set(handles.edit_sx,'String',num2str(handles.sx))
    set(handles.edit_sy,'String',num2str(handles.sy))
    set(handles.edit_sz,'String',num2str(handles.sz))
else
    set(hObject,'Value',1)
end

% ------------------------------------------------------------------------------------
function radiobutton_interfero_Callback(hObject, eventdata, handles)
% poe por defeito o s do ERS
if get(hObject,'Value')
    set(handles.radiobutton_deformation,'Value',0)
    set(handles.edit_sx,'String','0.333')
    set(handles.edit_sy,'String','-0.07')
    set(handles.edit_sz,'String','0.94')
else
    set(hObject,'Value',1)
end

% ------------------------------------------------------------------------------------
function pushbutton_cancel_Callback(hObject, eventdata, handles)
delete(handles.figure1)

% -----------------------------------------------------------------------------------------
function len = LineLength(h,geog)
x = get(h,'XData');     y = get(h,'YData');
len = [];
if (~iscell(x))
	if (geog)
        D2R = pi/180;    earth_rad = 6371;
        x = x * D2R;    y = y * D2R;
        lat_i = y(1:length(y)-1);   lat_f = y(2:length(y));     clear y;
        lon_i = x(1:length(x)-1);   lon_f = x(2:length(x));     clear x;
        tmp = sin(lat_i).*sin(lat_f) + cos(lat_i).*cos(lat_f).*cos(lon_f-lon_i);
        clear lat_i lat_f lon_i lon_f;
        len = [len; acos(tmp) * earth_rad];         % Distance in km
	else
        dx = diff(x);   dy = diff(y);
        len = [len; sqrt(dx.*dx + dy.*dy)];         % Distance in user unites
	end
else
	if (geog)
        D2R = pi/180;    earth_rad = 6371;
        for (k=1:length(x))
            xx = x{k} * D2R;    yy = y{k} * D2R;
            lat_i = yy(1:length(yy)-1);   lat_f = yy(2:length(yy));
            lon_i = xx(1:length(xx)-1);   lon_f = xx(2:length(xx));
            tmp = sin(lat_i).*sin(lat_f) + cos(lat_i).*cos(lat_f).*cos(lon_f-lon_i);
            len{k} = acos(tmp) * earth_rad;         % Distance in km
        end
	else
        for (k=1:length(x))
            xx = x{k};      yy = y{k};
            dx = diff(xx);  dy = diff(yy);
            len{k} = sqrt(dx.*dx + dy.*dy);         % Distance in user unites
        end
	end
end

% -----------------------------------------------------------------------------------------
function xy = read_xy(file,n_col,n_head)
% build the format string to read the data n_columns
xy = [];    format = [];    fid = fopen(file,'r');
for (i=1:n_col),    format = [format '%f '];    end
% Jump header lines
for (i = 1:n_head),    tline = fgetl(fid);  end

todos = fread(fid,'*char');
xy = sscanf(todos,format,[n_col inf])';    % After hours strugling agains this FILHO DA PUTA, I may have found
fclose(fid);

% -----------------------------------------------------------------------------------------
function checkbox_hideFaultPlanes_Callback(hObject, eventdata, handles)
if (handles.n_faults > 1)   fault = get(handles.popup_fault,'Value');
else                        fault = 1;      end
hp = getappdata(handles.h_fault(fault),'PatchHand');
if (get(hObject,'Value'))
    try,    set(hp,'Visible','off');    end
    handles.hide_planes(fault) = 1;
else
    try,    set(hp,'Visible','on');     end
    handles.hide_planes(fault) = 0;
end
guidata(hObject,handles)

% -----------------------------------------------------------------------------------------
function handles = set_all_faults(handles,varargin);
handles.h_calling_fig = varargin{1};
handles.geog = varargin{2};
handles.h_fault = varargin{3};
handles.FaultStrike = varargin{4};
handles.DislocStrike = varargin{4};
handles.FaultTopDepth = varargin{5};
FaultSlip = varargin{6};
handles.FaultDip  = varargin{7};
handles.DislocRake = varargin{8};

n_fault = length(handles.h_fault);      n_seg = length(handles.FaultStrike{1});

for (k=1:n_fault)
    handles.FaultWidth{k} = repmat(12,n_seg,1);
    %handles.DislocSlip{k}(1:n_seg) = FaultSlip{k}(1:n_seg)/100;
    handles.DislocSlip{k}(1:n_seg) = FaultSlip{k}(1:n_seg);
    handles.uz{k}(1:n_seg) = 0;
end

for (k=1:n_fault)   % Keep it scalar to avoid dimensional multiplication problems
    for (i=1:n_seg)
        handles.FaultDepth{k}(i) = handles.FaultTopDepth{k}(i) + handles.FaultWidth{k}(i) * ...
            sin(handles.FaultDip{k}(i) * pi/180);
        handles.ux{k}(i) = handles.DislocSlip{k}(i) * cos(handles.DislocRake{k}(i) * pi/180);
        handles.uy{k}(i) = handles.DislocSlip{k}(i) * sin(handles.DislocRake{k}(i) * pi/180);
    end
end

set(handles.edit_FaultWidth,'String','8');      % FAULT WIDTH - SHOULD NOT BE FIX
set(handles.edit_FaultStrike,'String',num2str(handles.FaultStrike{1}(1)));
set(handles.edit_FaultDip,'String',num2str(handles.FaultDip{1}(1)));
set(handles.edit_FaultTopDepth,'String',num2str(handles.FaultTopDepth{1}(1)));
set(handles.edit_FaultDepth,'String',num2str(handles.FaultDepth{1}(1)));
set(handles.edit_DislocSlip,'String',num2str(handles.DislocSlip{1}(1)));
set(handles.edit_DislocStrike,'String',num2str(handles.FaultStrike{1}(1)));
set(handles.edit_DislocRake,'String',num2str(handles.DislocRake{1}(1)));
set(handles.edit_ux,'String',num2str(handles.ux{1}(1)));
set(handles.edit_uy,'String',num2str(handles.uy{1}(1)));
set(handles.edit_uz,'String','0');


% --- Creates and returns a handle to the GUI figure. 
function deform_okada_LayoutFcn(h1,handles);

set(h1,'Units','pixels',...
'PaperUnits',get(0,'defaultfigurePaperUnits'),...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'MenuBar','none',...
'Name','Okada deformation',...
'NumberTitle','off',...
'Position',[520 415 536 385],...
'Renderer',get(0,'defaultfigureRenderer'),...
'RendererMode','manual',...
'Resize','off',...
'Tag','figure1',...
'HandleVisibility','callback',...
'UserData',[]);

h2 = uicontrol('Parent',h1,'Position',[320 224 211 131],...
'String',{''},'Style','frame','Tag','frame4');

h3 = uicontrol('Parent',h1,'Position',[10 224 181 131],...
'String',{''},'Style','frame','Tag','frame3');

h4 = uicontrol('Parent',h1,'Enable','inactive','Position',[11 15 350 67],...
'String',{''},'Style','frame','Tag','frame2');

h5 = uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@deform_okada_uicallback,h1,'edit_FaultLength_Callback'},...
'Position',[20 311 71 21],...
'Style','edit',...
'TooltipString','Fault length (km)',...
'Tag','edit_FaultLength');

h6 = uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@deform_okada_uicallback,h1,'edit_FaultWidth_Callback'},...
'Position',[110 311 71 21],...
'Style','edit',...
'TooltipString','Fault width (km)',...
'Tag','edit_FaultWidth');

h7 = uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@deform_okada_uicallback,h1,'edit_FaultStrike_Callback'},...
'Position',[20 271 71 21],...
'Style','edit',...
'TooltipString','Fault strike (degrees)',...
'Tag','edit_FaultStrike');

h8 = uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@deform_okada_uicallback,h1,'edit_FaultDip_Callback'},...
'Position',[110 271 71 21],...
'Style','edit',...
'TooltipString','Fault dip (degrees)',...
'Tag','edit_FaultDip');

h9 = uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@deform_okada_uicallback,h1,'edit_FaultDepth_Callback'},...
'Position',[20 232 71 21],...
'Style','edit',...
'TooltipString','Depth of the base of fault''s plane',...
'Tag','edit_FaultDepth');

h10 = uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@deform_okada_uicallback,h1,'edit_FaultTopDepth_Callback'},...
'Position',[110 231 71 21],...
'Style','edit',...
'TooltipString','Alternatively, give depth to the fault''s top ',...
'Tag','edit_FaultTopDepth');

h11 = uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@deform_okada_uicallback,h1,'edit_DislocStrike_Callback'},...
'Position',[330 311 51 21],...
'Style','edit',...
'Tag','edit_DislocStrike');

h12 = uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@deform_okada_uicallback,h1,'edit_DislocRake_Callback'},...
'Position',[400 311 51 21],...
'Style','edit',...
'TooltipString','Displacement angle clock-wise from horizontal',...
'Tag','edit_DislocRake');

h13 = uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@deform_okada_uicallback,h1,'edit_DislocSlip_Callback'},...
'Position',[470 311 51 21],...
'Style','edit',...
'TooltipString','Total displacement',...
'Tag','edit_DislocSlip');

h14 = uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@deform_okada_uicallback,h1,'edit_ux_Callback'},...
'Position',[330 271 51 21],...
'Style','edit',...
'TooltipString','Left-lateral displacement along the fault plane (along strike)',...
'Tag','edit_ux');

h15 = uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@deform_okada_uicallback,h1,'edit_uy_Callback'},...
'Position',[400 271 51 21],...
'Style','edit',...
'TooltipString','Displacement up-dip the fault plane (across strike)',...
'Tag','edit_uy');

h16 = uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@deform_okada_uicallback,h1,'edit_uz_Callback'},...
'Position',[470 271 51 21],...
'Style','edit',...
'TooltipString','fault tensile slip',...
'Tag','edit_uz');

h17 = uicontrol('Parent',h1,'Enable','inactive','Position',[10 109 350 93],...
'String',{''},'Style','frame','Tag','frame1');

h18 = uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@deform_okada_uicallback,h1,'edit_sx_Callback'},...
'Position',[330 231 51 21],...
'Style','edit',...
'TooltipString','Component of unit vector along North coords',...
'Tag','edit_sx');

h19 = uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@deform_okada_uicallback,h1,'edit_sy_Callback'},...
'Position',[400 231 51 21],...
'Style','edit',...
'TooltipString','Component of unit vector along East coords',...
'Tag','edit_sy');

h20 = uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@deform_okada_uicallback,h1,'edit_sz_Callback'},...
'Position',[470 231 51 21],...
'Style','edit',...
'TooltipString','Component of unit vector along Vertical coords',...
'Tag','edit_sz');

h21 = uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@deform_okada_uicallback,h1,'edit_Xmin_Callback'},...
'HorizontalAlignment','left',...
'Position',[76 162 71 21],...
'Style','edit',...
'TooltipString','X min value',...
'Tag','edit_Xmin');

h22 = uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@deform_okada_uicallback,h1,'edit_Xmax_Callback'},...
'HorizontalAlignment','left',...
'Position',[152 162 71 21],...
'Style','edit',...
'TooltipString','X max value',...
'Tag','edit_Xmax');

h23 = uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@deform_okada_uicallback,h1,'edit_Xinc_Callback'},...
'HorizontalAlignment','left',...
'Position',[228 162 71 21],...
'Style','edit',...
'TooltipString','DX grid spacing',...
'Tag','edit_Xinc');

h24 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@deform_okada_uicallback,h1,'edit_Ncols_Callback'},...
'Position',[304 162 45 21],...
'Style','edit',...
'TooltipString','Number of columns in the grid',...
'Tag','edit_Ncols');

h25 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@deform_okada_uicallback,h1,'edit_Ymin_Callback'},...
'HorizontalAlignment','left',...
'Position',[76 136 71 21],...
'Style','edit',...
'TooltipString','Y min value',...
'Tag','edit_Ymin');

h26 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@deform_okada_uicallback,h1,'edit_Ymax_Callback'},...
'HorizontalAlignment','left',...
'Position',[152 136 71 21],...
'Style','edit',...
'TooltipString','Y max value',...
'Tag','edit_Ymax');

h27 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@deform_okada_uicallback,h1,'edit_Yinc_Callback'},...
'HorizontalAlignment','left',...
'Position',[228 136 71 21],...
'Style','edit',...
'TooltipString','DY grid spacing',...
'Tag','edit_Yinc');

h28 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@deform_okada_uicallback,h1,'edit_Nrows_Callback'},...
'Position',[304 136 45 21],...
'Style','edit',...
'TooltipString','Number of columns in the grid',...
'Tag','edit_Nrows');

h29 = uicontrol('Parent',h1,...
'BackgroundColor',[0.831372559070587 0.815686285495758 0.7843137383461],...
'Callback',{@deform_okada_uicallback,h1,'pushbutton_Help_R_Callback'},...
'FontWeight','bold',...
'ForegroundColor',[0 0 1],...
'Position',[289 114 61 18],...
'String','?',...
'Tag','pushbutton_Help_R');

h30 = uicontrol('Parent',h1,...
'Callback',{@deform_okada_uicallback,h1,'checkbox_Option_H_Callback'},...
'CData',[],...
'Position',[22 54 65 15],...
'String','Headers?',...
'Style','checkbox',...
'TooltipString','Are there any header lines in the input file?',...
'Tag','checkbox_Option_H',...
'UserData',[]);

h31 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@deform_okada_uicallback,h1,'edit_nHeaders_Callback'},...
'HorizontalAlignment','left',...
'Position',[171 50 31 20],...
'String','1',...
'Style','edit',...
'TooltipString','How many?',...
'Tag','edit_nHeaders');

h32 = uicontrol('Parent',h1,...
'Callback',{@deform_okada_uicallback,h1,'checkbox_ToggleXY_Callback'},...
'Position',[221 53 75 19],...
'String','Toggle x,y',...
'Style','checkbox',...
'TooltipString','Toggle x and y columns',...
'Tag','checkbox_ToggleXY');

h33 = uicontrol('Parent',h1,'Enable','inactive','Position',[18 167 55 15],...
'String','X Direction','Style','text','Tag','text1');

h34 = uicontrol('Parent',h1,'Enable','inactive','Position',[17 141 55 15],...
'String','Y Direction','Style','text','Tag','text2');

h35 = uicontrol('Parent',h1,'Enable','inactive','Position',[169 184 41 13],...
'String','Max','Style','text','Tag','text3');

h36 = uicontrol('Parent',h1,'Enable','inactive','Position',[91 185 41 13],...
'String','Min','Style','text','Tag','text4');

h37 = uicontrol('Parent',h1,'Enable','inactive','Position',[246 185 41 13],...
'String','Spacing','Style','text','Tag','text5');

h38 = uicontrol('Parent',h1,'Enable','inactive','Position',[302 185 51 13],...
'String','# of lines','Style','text','Tag','text6','UserData',[]);

h39 = uicontrol('Parent',h1,'Enable','inactive','Position',[30 195 121 15],...
'String','Griding Line Geometry','Style','text','Tag','text7');

h40 = uicontrol('Parent',h1,...
'BackgroundColor',[0.831372559070587 0.815686285495758 0.7843137383461],...
'Callback',{@deform_okada_uicallback,h1,'pushbutton_Help_H_Callback'},...
'FontWeight','bold',...
'ForegroundColor',[0 0 1],...
'Position',[331 50 22 22],...
'String','?',...
'Tag','pushbutton_Help_H');

h41 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@deform_okada_uicallback,h1,'edit_InputFile_Callback'},...
'HorizontalAlignment','left',...
'Position',[21 22 310 22],...
'Style','edit',...
'TooltipString','File name with x, y positions where to compute deformation',...
'Tag','edit_InputFile');

h42 = uicontrol('Parent',h1,...
'Callback',{@deform_okada_uicallback,h1,'pushbutton_InputFile_Callback'},...
'Position',[331 21 23 23],...
'Tag','pushbutton_InputFile');

h43 = uicontrol('Parent',h1,'Enable','inactive','Position',[32 74 145 15],...
'String','Input Ground Positions File','Style','text','Tag','text8');

h44 = uicontrol('Parent',h1,'HorizontalAlignment','left','Position',[102 53 67 15],...
'String','Nº of headers','Style','text','TooltipString','How many?','Tag','text9');

h45 = uicontrol('Parent',h1,'Enable','inactive','Position',[36 333 41 13],...
'String','Length','Style','text','Tag','text10');

h46 = uicontrol('Parent',h1,'Enable','inactive','Position',[125 334 41 13],...
'String','Width','Style','text','Tag','text11');

h47 = uicontrol('Parent',h1,'Enable','inactive','Position',[34 293 41 13],...
'String','Strike','Style','text','Tag','text12');

h48 = uicontrol('Parent',h1,'Enable','inactive','Position',[124 293 41 13],...
'String','Dip','Style','text','Tag','text13');

h49 = uicontrol('Parent',h1,'Enable','inactive','Position',[108 252 75 16],...
'String','Depth to Top','Style','text',...
'TooltipString','Depth to the top of the fault (>= 0)',...
'Tag','text14');

h50 = uicontrol('Parent',h1,'Enable','inactive','Position',[335 333 41 13],...
'String','Strike','Style','text','Tag','text15');

h51 = uicontrol('Parent',h1,'Enable','inactive','Position',[404 333 41 13],...
'String','Rake','Style','text','Tag','text16');

h52 = uicontrol('Parent',h1,'Enable','inactive','Position',[474 333 41 13],...
'String','Slip','Style','text','Tag','text17');

h53 = uicontrol('Parent',h1,...
'Callback',{@deform_okada_uicallback,h1,'radiobutton_deformation_Callback'},...
'Position',[420 115 79 15],...
'String','Deformation',...
'Style','radiobutton',...
'Value',1,...
'Tag','radiobutton_deformation');

h54 = uicontrol('Parent',h1,...
'Callback',{@deform_okada_uicallback,h1,'radiobutton_interfero_Callback'},...
'Position',[420 89 83 15],...
'String','Interferogram',...
'Style','radiobutton',...
'Tag','radiobutton_interfero');

h55 = uicontrol('Parent',h1,...
'Callback',{@deform_okada_uicallback,h1,'pushbutton_compute_Callback'},...
'FontWeight','bold',...
'Position',[420 53 71 23],...
'String','Compute',...
'Tag','pushbutton_compute');

h56 = uicontrol('Parent',h1,'Enable','inactive','Position',[34 253 41 16],...
'String','Depth','Style','text',...
'TooltipString','Depth to the top of the fault (>= 0)','Tag','text18');

h57 = uicontrol('Parent',h1,'Enable','inactive','Position',[346 293 21 13],...
'String','u1','Style','text','Tag','text19');

h58 = uicontrol('Parent',h1,'Enable','inactive','Position',[415 293 21 13],...
'String','u2','Style','text','Tag','text20');

h59 = uicontrol('Parent',h1,'Enable','inactive','Position',[485 293 21 13],...
'String','u3','Style','text','Tag','text21');

h60 = uicontrol('Parent',h1,'Enable','inactive','Position',[345 253 21 13],...
'String','Sn','Style','text','Tag','text22');

h61 = uicontrol('Parent',h1,'Enable','inactive','Position',[414 253 21 13],...
'String','Se','Style','text','Tag','text23');

h62 = uicontrol('Parent',h1,'Enable','inactive','Position',[484 253 21 13],...
'String','Sz','Style','text','Tag','text24');

h63 = uicontrol('Parent',h1,...
'Callback',{@deform_okada_uicallback,h1,'pushbutton_cancel_Callback'},...
'FontWeight','bold',...
'Position',[420 22 71 23],...
'String','Cancel',...
'Tag','pushbutton_cancel');

h64 = uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@deform_okada_uicallback,h1,'popup_segment_Callback'},...
'Position',[210 319 91 22],...
'Style','popupmenu',...
'TooltipString','Set parameters with respect to this segment',...
'Value',1,'Tag','popup_segment');

h65 = uicontrol('Parent',h1,'Enable','inactive','Position',[225 342 57 15],...
'String','Segment','Style','text','Tag','fault_segment');

h66 = uicontrol('Parent',h1,'Enable','inactive','Position',[53 348 85 15],...
'String','Fault Geometry','Style','text','Tag','text26');

h67 = uicontrol('Parent',h1,'Enable','inactive','Position',[373 348 111 15],...
'String','Dislocation Geometry','Style','text','Tag','text27');

h68 = uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@deform_okada_uicallback,h1,'popup_fault_Callback'},...
'Position',[210 273 91 22],...
'Style','popupmenu',...
'TooltipString','Toggle between faults',...
'Value',1,'Tag','popup_fault');

h69 = uicontrol('Parent',h1,'Enable','inactive','Position',[236 295 32 15],...
'String','Faults','Style','text','Tag','fault_number');

h70 = uicontrol('Parent',h1,...
'Callback',{@deform_okada_uicallback,h1,'checkbox_hideFaultPlanes_Callback'},...
'Position',[406 168 104 17],...
'String','Hide fault plane',...
'Style','checkbox','Tag','checkbox_hideFaultPlanes');

uicontrol('Parent',h1,'Position',[225 248 55 15],'ForegroundColor',[1 0 0],...
'String','CONFIRM','Style','text','Tag','text22');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@deform_okada_uicallback,h1,'popup_GridCoords_Callback'},...
'String', {'Geogs' 'Meters' 'Kilometers'},...
'Position',[210 225 91 22],'Style','popupmenu',...
'TooltipString','GRID COORDINATES: IT IS YOUR RESPONSABILITY THAT THIS IS CORRECT',...
'Value',1,'Tag','popup_GridCoords');

function deform_okada_uicallback(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
feval(callback_name,hObject,[],guidata(h1));

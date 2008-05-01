function varargout = compute_euler(varargin)
% M-File changed by desGUIDE 
% varargin   command line arguments to compute_euler (see VARARGIN)

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
compute_euler_LayoutFcn(hObject,handles);
handles = guihandles(hObject);
movegui(hObject,'west')
 
% Initialize those
handles.h_calling_fig = [];
handles.LonRange = 30;
handles.LatRange = 30;
handles.AngRange = 4;
handles.Nintervals = 20;
handles.isoca1 = [];
handles.isoca2 = [];
handles.pLon_ini = [];
handles.pLat_ini = [];
handles.pAng_ini = [];
handles.do_graphic = 0;
handles.path_continent = [pwd filesep 'continents' filesep];
set(handles.slider_wait,'Max',handles.Nintervals^2)

if (length(varargin) == 1)
    handles.h_calling_fig = varargin{1};        % This the Mirone's fig handle
else
    errordlg('COMPUTE EULER: wrong number of arguments.','Error')
    delete(hObject);    return
end

handMir = guidata(handles.h_calling_fig);
if (handMir.no_file)
    errordlg('You didn''t even load a file. What are you expecting then?','ERROR')
    delete(hObject);    return
end
if (~handMir.geog)
    errordlg('This operation is currently possible only for geographic type data','ERROR')
    delete(hObject);    return
end

%--------------- Give a Pro look (3D) to the frame boxes -------------------------
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

% Recopy the text fields on top of previously created frames (uistack is to damn slow)
h_t = findobj(hObject,'Style','Text');
for i=1:length(h_t)
    usr_d = get(h_t(i),'UserData');
    t_size = get(h_t(i),'Position');   t_str = get(h_t(i),'String');    fw = get(h_t(i),'FontWeight');
    bgc = get (h_t(i),'BackgroundColor');   fgc = get (h_t(i),'ForegroundColor');
    t_just = get(h_t(i),'HorizontalAlignment');     t_tag = get (h_t(i),'Tag');
    uicontrol('Parent',hObject, 'Style','text', 'Position',t_size,'String',t_str,'Tag',t_tag, ...
        'BackgroundColor',bgc,'ForegroundColor',fgc,'FontWeight',fw,...
        'UserData',usr_d,'HorizontalAlignment',t_just);
end
delete(h_t)
%------------------- END Pro look (3D) ----------------------------------------------------------

% Choose default command line output for compute_euler_export
handles.output = hObject;
guidata(hObject, handles);

set(hObject,'Visible','on');
% NOTE: If you make uiwait active you have also to uncomment the next three lines
% handles = guidata(hObject);
% out = compute_euler_OutputFcn(hObject, [], handles);
% varargout{1} = out;

% --- Outputs from this function are returned to the command line.
function varargout = compute_euler_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% handles    structure with handles and user data (see GUIDATA)
% Get default command line output from handles structure
varargout{1} = handles.output;

% -------------------------------------------------------------------------------------
function edit_first_file_Callback(hObject, eventdata, handles)
fname = get(hObject,'String');
if isempty(fname)    return;    end
% Let the pushbutton_first_file_Callback do all the work
pushbutton_first_file_Callback(hObject,[],guidata(gcbo),fname)

% -------------------------------------------------------------------------------------
function pushbutton_first_file_Callback(hObject, eventdata, handles,opt)
if (nargin == 4)    fname = opt;
else                opt = [];
end

if (~isempty(handles.h_calling_fig))                    % If we know the handle to the calling fig
    cfig_handles = guidata(handles.h_calling_fig);      % get handles of the calling fig
    last_dir = cfig_handles.last_dir;
    home = cfig_handles.home_dir;
else
    last_dir = [];
end

if (isempty(opt))    % Otherwise we already know fname from the 4th input argument
    if (~isempty(last_dir)),    cd(last_dir);   end
    [FileName,PathName] = uigetfile({'*.dat;*.DAT', 'Mag file (*.dat,*.DAT)';'*.*', 'All Files (*.*)'},'Select file');
    pause(0.01);
    if (~isempty(last_dir)),    cd(home);   end
    if isequal(FileName,0);     return;     end
    fname = [PathName FileName];
end
handles.isoca1 = le_fiche(fname);      % A maluca
handles.isoca1(:,2) = geog2auth(handles.isoca1(:,2));   % Convert to authalic lats
guidata(hObject, handles);

% -------------------------------------------------------------------------------------
function edit_second_file_Callback(hObject, eventdata, handles)
fname = get(hObject,'String');
if isempty(fname)    return;    end
% Let the pushbutton_first_file_Callback do all the work
pushbutton_second_file_Callback(hObject,[],guidata(gcbo),fname)

% -------------------------------------------------------------------------------------
function pushbutton_second_file_Callback(hObject, eventdata, handles, opt)
if (nargin == 4)    fname = opt;
else                opt = [];
end

if (~isempty(handles.h_calling_fig))                    % If we know the handle to the calling fig
    cfig_handles = guidata(handles.h_calling_fig);      % get handles of the calling fig
    last_dir = cfig_handles.last_dir;
    home = cfig_handles.home_dir;
else
    last_dir = [];
end

if (isempty(opt))    % Otherwise we already know fname from the 4th input argument
    if (~isempty(last_dir)),    cd(last_dir);   end
    [FileName,PathName] = uigetfile({'*.dat;*.DAT', 'Mag file (*.dat,*.DAT)';'*.*', 'All Files (*.*)'},'Select file');
    pause(0.01);
    if (~isempty(last_dir)),    cd(home);   end
    if isequal(FileName,0);     return;     end
    fname = [PathName FileName];
end
handles.isoca2 = le_fiche(fname);      % The fixed line
handles.isoca2(:,2) = geog2auth(handles.isoca2(:,2));   % Convert to authalic lats
guidata(hObject, handles);

% -------------------------------------------------------------------------------------
function edit_pLon_ini_Callback(hObject, eventdata, handles)
    handles.pLon_ini = str2double(get(hObject,'String'));
    guidata(hObject, handles);

% -------------------------------------------------------------------------------------
function edit_pLat_ini_Callback(hObject, eventdata, handles)
    handles.pLat_ini = str2double(get(hObject,'String'));
    guidata(hObject, handles);

% -------------------------------------------------------------------------------------
function edit_pAng_ini_Callback(hObject, eventdata, handles)
    handles.pAng_ini = str2double(get(hObject,'String'));
    guidata(hObject, handles);

% -------------------------------------------------------------------------------------
function pushbutton_polesList_Callback(hObject, eventdata, handles)
fid = fopen([handles.path_continent 'lista_polos.dat'],'rt');
c = fread(fid,'*char').';
fclose(fid);
s = strread(c,'%s','delimiter','\n');

[s,v] = choosebox('Name','One Euler list',...
                    'PromptString','List of poles:',...
                    'SelectString','Selected poles:',...
                    'ListSize',[380 300],'ListString',s);

if (v == 1)         % Finite pole
    handles.pLon_ini = s(1);
    handles.pLat_ini = s(2);
    handles.pAng_ini = s(3);
    set(handles.edit_pLon_ini, 'String', num2str(s(1)))
    set(handles.edit_pLat_ini, 'String', num2str(s(2)))
    set(handles.edit_pAng_ini, 'String', num2str(s(3)))
    guidata(hObject,handles)
elseif (v == 2)     % Stage poles
    %set(handles.edit_polesFile,'String',s)
end

% -------------------------------------------------------------------------------------
function edit_LonRange_Callback(hObject, eventdata, handles)
handles.LonRange = str2double(get(hObject,'String'));
guidata(hObject, handles);

% -------------------------------------------------------------------------------------
function edit_Nintervals_Callback(hObject, eventdata, handles)
handles.Nintervals = str2double(get(hObject,'String'));
set(handles.slider_wait,'Max',handles.Nintervals^2)
guidata(hObject, handles);

% -------------------------------------------------------------------------------------
function edit_LatRange_Callback(hObject, eventdata, handles)
handles.LatRange = str2double(get(hObject,'String'));
guidata(hObject, handles);

% -------------------------------------------------------------------------------------
function edit_AngRange_Callback(hObject, eventdata, handles)
handles.AngRange = str2double(get(hObject,'String'));
guidata(hObject, handles);

% -------------------------------------------------------------------------------------
function togglebutton_pickLines_Callback(hObject, eventdata, handles)
if (get(hObject,'Value'))
    % Test if we have potential target lines and their type
    h_mir_lines = findobj(handles.h_calling_fig,'Type','line');     % Fish all objects of type line in Mirone figure
    if (isempty(h_mir_lines))                                       % We don't have any lines
        str = ['If you hited this button on purpose, than you deserve the following insult.',...
                'You #!|"*!%!?~^)--$&.',... 
                'THERE ARE NO LINES IN THAT FIGURE.'];
        errordlg(str,'Chico Clever');   set(hObject,'Value',0);     return;
    end
    if (length(h_mir_lines) == 1)                                    % We don't have at least two lines
        str = ['If you hited this button on purpose, than you deserve the following insult.',...
                'You -$&#!*!%!?~^)-.|"/',... 
                'THERE IS ONLY ONE LINE IN THAT FIGURE.'];
        errordlg(str,'Chico Clever');   set(hObject,'Value',0);     return;
    end
    
    % The above test is not enough. For exemple, coastlines are not eligible neither,
    % but is very cumbersome to test all the possibilities of pure non-eligible lines.
    set(handles.h_calling_fig,'pointer','crosshair')
    h_line = get_polygon(handles.h_calling_fig);          % Get first line handle
    if (~isempty(h_line))
        x = get(h_line,'XData');    y = get(h_line,'YData');
        %y = geog2auth(y);                       % Convert to authalic latitudes
        handles.isoca1 = [x(:) y(:)];
        set(handles.edit_first_file,'String','Got left line','FontAngle','italic')
    else
        handles.isoca1 = [];
        set(handles.edit_first_file,'String','')
    end
    h_line = get_polygon(handles.h_calling_fig);          % Get second line handle
    if (~isempty(h_line))
        x = get(h_line,'XData');    y = get(h_line,'YData');
        %y = geog2auth(y);                       % Convert to authalic latitudes
        handles.isoca2 = [x(:) y(:)];
        set(handles.edit_second_file,'String','Got right line','FontAngle','italic')
    else
        handles.isoca2 = [];
        set(handles.edit_second_file,'String','')
    end
    set(handles.h_calling_fig,'pointer','arrow')
    if (isempty(handles.isoca1) | isempty(handles.isoca2))
        set(hObject,'Value',0)
        handles.do_graphic = 0;
    else
        handles.do_graphic = 1;
    end
    set(hObject,'Value',0)
    figure(handles.figure1)         % Bring this figure to front again
else        % What should I do?
    %handles.do_graphic = 0;
end
guidata(hObject, handles);

% -------------------------------------------------------------------------------
function slider_wait_Callback(hObject, eventdata, handles)
% Nothing to do here. Moreover, the slider 'Enable' prop is 'off'

% -------------------------------------------------------------------------------
function pushbutton_stop_Callback(hObject, eventdata, handles)
    set(handles.slider_wait,'Value',0)  % We have to do this first
    set(handles.slider_wait,'Max',1)    % This will signal the fit_pEuler function to stop

% -------------------------------------------------------------------------------------
function pushbutton_cancel_Callback(hObject, eventdata, handles)
    delete(handles.figure1)

% -------------------------------------------------------------------------------------
function pushbutton_compute_Callback(hObject, eventdata, handles)
% OK. See if we have all the information needed to compute the Euler pole

set(handles.edit_BFresidue,'String','');        set(handles.edit_InitialResidue,'String','')
set(handles.edit_pLon_fim,'String','');         set(handles.edit_pLat_fim,'String','')
set(handles.edit_pAng_fim,'String','')

if (isempty(handles.isoca1) | isempty(handles.isoca2))
    errordlg('Compute Euler pole with what? It would help if you provide me TWO lines.','Chico Clever')
    return
end
if (isempty(handles.pLon_ini) | isempty(handles.pLat_ini) | isempty(handles.pAng_ini))
    errordlg(['I need a first guess of the Euler pole you are seeking for.' ...
    'Pay attention to the "Starting Pole Section"'],'Error')
    return
end

calca_pEuler(handles)

% -------------------------------------------------------------------------------
function numeric_data = le_fiche(fname)
[bin,n_column,multi_seg,n_headers] = guess_file(fname);
% If error in reading file
if isempty(bin) && isempty(n_column) && isempty(multi_seg) && isempty(n_headers)
    errordlg(['Error reading file ' fname],'Error');
    numeric_data = [];
    return
end
if (isempty(n_headers)),     n_headers = NaN;    end
if (multi_seg)
    [numeric_data,multi_segs_str,headerlines] = text_read(fname,NaN,n_headers,'>');
else
    [numeric_data,multi_segs_str,headerlines] = text_read(fname,NaN,n_headers);
end

% -------------------------------------------------------------------------------
function calca_pEuler(handles)

D2R = pi / 180;
if (handles.do_graphic)     % Create a empty line handle
    h_line = line('parent',get(handles.h_calling_fig,'CurrentAxes'),'XData',[],'YData',[], ...
        'LineStyle','-.','LineWidth',2,'Tag','Fitted Line','Userdata',1);
    %draw_funs(h_line,'isochron',{'Fitted Line'})
end

% Compute the initial misfit with Bruno merit function
% [rlon,rlat] = rot_euler(handles.isoca1(:,1),handles.isoca1(:,2),handles.pLon_ini,handles.pLat_ini,handles.pAng_ini);
% [teta,lambda] = mincoordspher(handles.pLon_ini*D2R,handles.pLat_ini*D2R,handles.isoca2(:,1)*D2R, ...
%     handles.isoca2(:,2)*D2R,rlon*D2R,rlat*D2R);
% sum2 = sum( log(cosh(1*teta).*cosh(lambda)) );
% [rlon,rlat] = rot_euler(handles.isoca2(:,1),handles.isoca2(:,2),handles.pLon_ini,handles.pLat_ini,-handles.pAng_ini);
% [teta,lambda] = mincoordspher(handles.pLon_ini*D2R,handles.pLat_ini*D2R,handles.isoca1(:,1)*D2R, ...
%     handles.isoca1(:,2)*D2R,rlon*D2R,rlat*D2R);
% sum1 = sum( log(cosh(1*teta).*cosh(lambda)) );
% area0 = ((sum1 + sum2) / 2) / D2R;

% Compute the initial misfit with the polygArea merit function
% [rlon,rlat] = rot_euler(handles.isoca2(:,1),handles.isoca2(:,2),handles.pLon_ini,handles.pLat_ini,-handles.pAng_ini);
% sum1 = polygArea(handles.isoca1(:,1)*D2R, handles.isoca1(:,2)*D2R ,rlon*D2R, rlat*D2R);
% [rlon,rlat] = rot_euler(handles.isoca1(:,1),handles.isoca1(:,2),handles.pLon_ini,handles.pLat_ini,handles.pAng_ini);
% sum2 = polygArea(handles.isoca2(:,1)*D2R, handles.isoca2(:,2)*D2R, rlon*D2R, rlat*D2R);
% area0 = (sum1 + sum2) / 2;

% Compute the initial misfit with Miguel merit function
[rlon,rlat] = rot_euler(handles.isoca1(:,1),handles.isoca1(:,2),handles.pLon_ini,handles.pLat_ini,handles.pAng_ini);
dist1 = distmin(handles.isoca2(:,1)*D2R,handles.isoca2(:,2)*D2R,rlon*D2R,rlat*D2R);
sum1 = sum( dist1 ) / length(dist1);
dist2 = distmin(rlon*D2R,rlat*D2R,handles.isoca2(:,1)*D2R,handles.isoca2(:,2)*D2R);
sum2 = sum( dist2 ) / length(dist2);
area0 = (sum1 + sum2) / 2;

set(handles.edit_InitialResidue,'String',num2str(area0,'%.6f'));    pause(0.01)
set(h_line,'XData',rlon,'YData',rlat)

% Now comes the semi-brute force aproach to compute the pole
n = handles.Nintervals;         n1 = round(n/2);    n2 = n - n1;
dLon = handles.LonRange / 2;
dLat = handles.LatRange / 2;
dAng = handles.AngRange / 2;
p_lon = (handles.pLon_ini + linspace(-dLon,dLon,n)) * D2R;
p_lat = (handles.pLat_ini + linspace(-dLat,dLat,n)) * D2R;
% p_lat = (handles.pLat_ini + linspace(0,-dLat,n1));
% dx = abs(p_lat(2) - p_lat(1));
% p_lat = [p_lat [handles.pLat_ini + dx + linspace(0,dLat,n2)]] * D2R;
p_omeg = (handles.pAng_ini + linspace(-dAng,dAng,n)) * D2R;
%[p_lon,p_lat,p_omega,area_f] = fit_pEuler(handles,p_lon,p_lat,p_omeg,area0*D2R,h_line);
[p_lon,p_lat,p_omega,area_f] = fit_pEuler(handles,p_lon,p_lat,p_omeg,area0,h_line);

figure(handles.h_calling_fig)       % Bring the Mirone figure to front
draw_funs(h_line,'isochron',{'Fitted Line'})

if (area_f >= area0)
    msgbox(['I am quite good, but this time I couldn''t find a better pole than the one ',...
            'you gave me (Ah! probably it was me who computed it in a previous re-incarnation)'],'Report')
end

% -------------------------------------------------------------------------------
function [lon_bf,lat_bf,omega_bf,area_f] = fit_pEuler(handles,p_lon,p_lat,p_omeg,area0,h_line)
% This is the function that does the real work of computing the best Euler pole
% that fits the lines "isoca1" "isoca2"
lon_bf = [];        lat_bf = [];
omega_bf = [];      area_f = [];
%area0 = 1e6;
D2R = pi / 180;
n = length(p_lon);
isoca1 = handles.isoca1 * D2R;      % A maluca
isoca2 = handles.isoca2 * D2R;      % A fixa

%isoca2_lon_ct = isoca2(:,1) .* cos(isoca2(:,2));    % This doesn't change, so compute it only once here
isoca2_lon_ct = isoca2(:,1);

i_m = 1;    j_m = 1;    k_m = 1;    ij = 0;
for (i=1:n)                % Loop over lon
    if (get(handles.slider_wait,'Max') == 1)        % The STOP button was pushed. So, stop
        set(handles.slider_wait,'Max',handles.Nintervals^2)     % Reset it for the next run
        area_f = 0;     return;
    end
    for (j=1:n)            % Loop over lat
        ij = ij + 1;
        if (get(handles.slider_wait,'Max') == 1)    % The STOP button was pushed. So, stop
            set(handles.slider_wait,'Max',handles.Nintervals^2)     % Reset it for the next run
            area_f = 0;     return;
        end
        set(handles.slider_wait,'Value',ij);     pause(0.001);   % Otherwise the slider risks to not be updated
        count_ang = 0;
        for (k=1:n)        % Loop over omega            
%             [rlon,rlat] = rot_euler(isoca2(:,1),isoca2(:,2),p_lon(i),p_lat(j),-p_omeg(k),'radians');
%             [teta,lambda] = mincoordspher(p_lon(i),p_lat(j),isoca1(:,1),isoca1(:,2),rlon,rlat);
%             sum1 = sum( log(cosh(1*teta).*cosh(lambda)) );
%             [rlon,rlat] = rot_euler(isoca1(:,1),isoca1(:,2),p_lon(i),p_lat(j),p_omeg(k),'radians');
%             [teta,lambda] = mincoordspher(p_lon(i),p_lat(j),isoca2(:,1),isoca2(:,2),rlon,rlat);
%             sum2 = sum( log(cosh(1*teta).*cosh(lambda)) );
            
            % Isto e uma tentativa com a funcao custo do Miguel
            [rlon,rlat] = rot_euler(isoca1(:,1),isoca1(:,2),p_lon(i),p_lat(j),p_omeg(k),'radians');
            %dist1 = distmin(isoca2(:,1),isoca2(:,2),rlon,rlat);
            dist1 = distmin(isoca2_lon_ct,isoca2(:,2),rlon,rlat);
            sum1 = sum( dist1 ) / length(dist1);
            %dist2 = distmin(rlon,rlat,isoca2(:,1),isoca2(:,2));
            dist2 = distmin(rlon,rlat,isoca2_lon_ct,isoca2(:,2));
            sum2 = sum( dist2 ) / length(dist2);
            %sum2 = sum1;

%             % Isto e outra tentativa
%             [rlon,rlat] = rot_euler(isoca2(:,1),isoca2(:,2),p_lon(i),p_lat(j),-p_omeg(k),'radians');
%             sum1 = polygArea(isoca1(:,1),isoca1(:,2),rlon,rlat);
%             [rlon,rlat] = rot_euler(isoca1(:,1),isoca1(:,2),p_lon(i),p_lat(j),p_omeg(k),'radians');
%             sum2 = polygArea(isoca2(:,1),isoca2(:,2),rlon,rlat);
            
            area = (sum1 + sum2) / 2;
            if (area < area0)
                area0 = area;
                i_m = i;    j_m = j;    k_m = k;
                count_ang = 0;
                if (handles.do_graphic)
                    set(h_line,'XData',rlon/D2R,'YData',rlat/D2R)
                    pause(0.001)
                end
                set(handles.edit_pLon_fim,'String',num2str(p_lon(i) / D2R,'%.2f'))
                set(handles.edit_pLat_fim,'String',num2str(p_lat(j) / D2R,'%.2f'))
                set(handles.edit_pAng_fim,'String',num2str(p_omeg(k) / D2R,'%.2f'))
                %set(handles.edit_BFresidue,'String',num2str(area / D2R,'%.5f'))
                set(handles.edit_BFresidue,'String',num2str(area,'%.6f'))
            elseif (count_ang > 3)      % We had more than 3 attemps with icreasing both the angle and
                count_ang = 0;          % the error in the fit. So don't wast more time in this direction
                continue
            elseif (area - area0 > 0)           % Fit did not improve
                count_ang = count_ang + 1;      % Count number of sucessive fit non-improving with angle
            end
        end
    end
end

set(handles.slider_wait,'Value',0)      % Reset it for the next run
lon_bf = p_lon(i_m) / D2R;
lat_bf = p_lat(j_m) / D2R;
omega_bf = p_omeg(k_m) / D2R;
%area_f = area0 / D2R;
area_f = area0;

% -------------------------------------------------------------------------------
function area = polygArea(lon,lat,r_lon,r_lat)
% a = polygArea(x,y) returns twice the signed area of the polygons
% Reference: 
% http://geometryalgorithms.com/Archive/algorithm_0101/algorithm_0101.htm
% D2R = pi / 180;
% DEG2KM = 2 * pi * 6371.0087714 / 360.0;         % GRS-80 sphere km/degree
y = [lat; r_lat(end:-1:1)];
x = [lon; r_lon(end:-1:1)];
x = x .* cos(y);
x = x - mean(x);
n = length(x);
i = [2:n 1];
j = [3:n 1 2];
area = abs(sum(x(i) .* (y(j) - y([1:n])))) * 1e6;

% -------------------------------------------------------------------------------
function [teta,lambda] = mincoordspher(p_lon,p_lat,lon,lat,r_lon,r_lat)
% Given that this function is called many times, the angles (the arguments) are
% expected to have already been converted to radians.

D2R = pi / 180;
cosrf = sin(p_lat).*sin(lat) + cos(p_lat).*cos(lat).*cos(lon-p_lon);
sinrf = sqrt(1 - cosrf.^2);
rf = atan2(sinrf,cosrf);
cthf = sin(lat) - cosrf.*sin(p_lat);
sthf = cos(lat).*cos(p_lat).*sin(lon-p_lon);
tetaf = atan2(sthf,cthf);

cosrn = sin(p_lat).*sin(r_lat) + cos(p_lat).*cos(r_lat).*cos(r_lon-p_lon);
sinrn = sqrt(1 - cosrn.^2);
rn = atan2(sinrn,cosrn);
cthn = sin(r_lat) - cosrn.*sin(p_lat);
sthn = cos(r_lat).*cos(p_lat).*sin(r_lon-p_lon);
tetan = atan2(sthn,cthn);

% Now the *f and *n have different sizes
n_size = length(r_lon);
f_size = length(lon);
teta = repmat(0,f_size,1);      lambda = teta;
for (k=1:f_size)
    teta(k) = min( abs(rn - rf(k)) );
    lambda(k) = min( abs(rn .* (tetan - tetaf(k))) );
end
% teta = teta / D2R;
% lambda = lambda / D2R;
teta = teta;
lambda = lambda;


% -------------------------------------------------------------------------------
% subroutine distmin(tme,tmn,xvert,yvert,nvert,dist0)
function dist = distmin__(lon,lat,r_lon,r_lat)
D2R = pi / 180;
r_lon = r_lon .* cos(r_lat);      % VERY strange. This should improve the fit
lon = lon .* cos(lat);            % But, on the contrary, it degrades it ??
epsilon = 1e-3;
%eps2 = 1e-6;
dist = zeros(length(lon),1);
for (k=1:length(lon))
    dist0 = 1e20;
    for (j=1:length(r_lon)-1)
        x1 = r_lon(j);      x2 = r_lon(j+1);
        y1 = r_lat(j);      y2 = r_lat(j+1);
        dd = 1e20;
        while (dd > epsilon)
            d1 = sqrt((lon(k)-x1)^2+(lat(k)-y1)^2);
            d2 = sqrt((lon(k)-x2)^2+(lat(k)-y2)^2);
            dd = sqrt((x1-x2)^2+(y1-y2)^2);
            xm = (x1+x2) * 0.5;
            ym = (y1+y2) * 0.5;
            if (d1 < d2)
                x2 = xm;    y2 = ym;
            else
                x1 = xm;    y1 = ym;
            end
        end
        dist1 = (d1+d2) * 0.5;
        if(dist1 < dist0) dist0 = dist1;  end   % When true, the closest point on the other line was found
    end
    dist(k) = dist0;
end
% for (jj=1:nvert-1)
%     x1=xvert(jj);
%     y1=yvert(jj);
%     x2=xvert(jj+1);
%     y2=yvert(jj+1);
% 11  d1=sqrt((tme-x1).^2+(tmn-y1).^2);
%     d2=sqrt((tme-x2).^2+(tmn-y2).^2);
%     dd=sqrt((x1-x2).^2+(y1-y2).^2);
%     while ()
%    if(dd > eps)
%        xm=(x1+x2)/2;
%        ym=(y1+y2)/2;
%        if (d1 < d2)
%           x2=xm;
%           y2=ym;
%        else
%           x1=xm;
%           y1=ym;
%        end
%        go to 11
%    end
%    dist1=(d1+d2)/2;
%    if(dist0 > dist1) dist0=dist1;  end
% end

% -------------------------------------------------------------------------------
function dist = distmin(lon,lat,r_lon,r_lat)
r_lon = r_lon .* cos(r_lat);      % VERY strange. This should improve the fit
lon = lon .* cos(lat);            % But, on the contrary, it degrades it ??
r_lon = r_lon(:)';      r_lat = r_lat(:)';      % Make sure they are row vectors
eps1 = 1e-3;            eps2 = eps1 * eps1;
n_pt = numel(lon);      n_pt_rot = numel(r_lon);
dist = zeros(n_pt,1);
%Dsts = zeros(n_pt,n_pt_rot);    % (i,j) element will hold distance from vertex i of line (lon,lat)
                                % to vertex j of line (r_lon,r_lat)
outliners = false(n_pt,1);      % To remove vertex where lines do not intersect
for (k=1:n_pt)
    dist0 = 1e20;
    Dsts = sqrt((lon(k)-r_lon).^2+(lat(k)-r_lat).^2);
    [D,ind] = min(Dsts);
    if (ind == 1 || ind == n_pt_rot && n_pt > 4)    % This point is outside the lines intersection. Flag it to die.
        outliners(k) = true;
        continue
    end
    %for (j = ind-1:min(ind+1,n_pt_rot)-1)
	d1 = 0;		d2 = 0;		xm = 0;		ym = 0;
	for (j = ind-1:ind)
        x1 = r_lon(j);      x2 = r_lon(j+1);
        y1 = r_lat(j);      y2 = r_lat(j+1);
        dd = 1e20;
        while (dd > eps2)
            d1 = sqrt((lon(k)-x1)^2+(lat(k)-y1)^2);
            d2 = sqrt((lon(k)-x2)^2+(lat(k)-y2)^2);
            dd = (x1-x2)^2+(y1-y2)^2;
            xm = (x1+x2) * 0.5;
            ym = (y1+y2) * 0.5;
            if (d1 < d2)
                x2 = xm;    y2 = ym;
            else
                x1 = xm;    y1 = ym;
            end
        end
        dist1 = (d1+d2) * 0.5;
        if(dist1 < dist0) dist0 = dist1;  end   % When true, the closest point on the other line was found        
    end
    % 'dist' has the least distance of each vertex of line (lon,lat) to rotated line (r_lon,r_lat)
    dist(k) = dist0;
end

dist(outliners) = [];       % Delete points not belonging to the lines intersection zone
if (isempty(dist)),     dist = 0;   end     % Don't let it go empty

% --------------------------------------------------------------------------------
function lat = geog2auth(lat0)
% Convert geodetic (geographic) latitudes to authalic latitudes. 
% Degrees are assumed in input
lat0 = lat0 * pi / 180;
flat = 1/298.257222101;         % GRS80 flatness
ecc = sqrt(2*flat - flat.^2);   % exccentricity
f1 = ecc^2 /3 + 31*ecc^4 / 180 + 59*ecc^6 / 560;
f2 = 17*ecc^4 / 360 + 61*ecc^6 / 1260;
f3 = 383*ecc^6 / 45360;

lat = lat0 - f1*sin(2*lat0) + f2*sin(4*lat0) - f3*sin(6*lat0);
lat = lat * 180 / pi;  %  Convert back to degrees

% --------------------------------------------------------------------------------
function lat = auth2geog(lat0)
% Convert authalic  latitudes to geodetic (geographic) latitudes. 
% Degrees are assumed in input
lat0 = lat0 * pi / 180;
flat = 1/298.257222101;         % GRS80 flatness
ecc = sqrt(2*flat - flat.^2);   % exccentricity
f1 = ecc^2 / 3 + 31*ecc^4 / 180 + 517*ecc^6 / 5040;
f2 = 23*ecc^4 / 360 + 251*ecc^6 / 3780;
f3 = 761*ecc^6 / 45360;

lat = lat0 + f1*sin(2*lat0) + f2*sin(4*lat0) + f3*sin(6*lat0);
lat = lat * 180 / pi;  %  Convert back to degrees
     
% --- Creates and returns a handle to the GUI figure. 
function compute_euler_LayoutFcn(h1,handles);
set(h1,...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'MenuBar','none',...
'Name','compute_euler',...
'NumberTitle','off',...
'Position',[520 427 522 373],...
'RendererMode','manual',...
'Resize','off',...
'Tag','figure1',...
'UserData',[]);

h2 = uicontrol('Parent',h1,'Position',[10 49 501 110],'String',{''},'Style','frame','Tag','frame3');
h3 = uicontrol('Parent',h1,'Position',[10 182 501 64],'String',{''},'Style','frame','Tag','frame2');
h4 = uicontrol('Parent',h1,'Position',[10 266 501 101],'String',{''},'Style','frame','Tag','frame1');

h5 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@compute_euler_uicallback,h1,'edit_first_file_Callback'},...
'Position',[20 316 211 21],...
'Style','edit','Tag','edit_first_file');

h6 = uicontrol('Parent',h1,...
'Callback',{@compute_euler_uicallback,h1,'pushbutton_first_file_Callback'},...
'FontSize',10,...
'FontWeight','bold',...
'Position',[231 316 21 21],...
'String','...','Tag','pushbutton_first_file');

h7 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@compute_euler_uicallback,h1,'edit_second_file_Callback'},...
'Position',[267 316 211 21],...
'Style','edit','Tag','edit_second_file');

h8 = uicontrol('Parent',h1,...
'Callback',{@compute_euler_uicallback,h1,'pushbutton_second_file_Callback'},...
'FontSize',10,...
'FontWeight','bold',...
'Position',[478 316 21 21],...
'String','...','Tag','pushbutton_second_file');

h9 = uicontrol('Parent',h1,...
'Callback',{@compute_euler_uicallback,h1,'togglebutton_pickLines_Callback'},...
'Position',[186 279 141 23],...
'String','Pick lines from Figure',...
'Style','togglebutton',...
'TooltipString','Allows you to mouse select the two lines from a Mirone figure',...
'Tag','togglebutton_pickLines');

h10 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@compute_euler_uicallback,h1,'edit_pLon_ini_Callback'},...
'Position',[20 201 61 21],...
'Style','edit',...
'TooltipString','Start Euler pole Longitude',...
'Tag','edit_pLon_ini');

h11 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@compute_euler_uicallback,h1,'edit_pLat_ini_Callback'},...
'Position',[110 201 61 21],...
'Style','edit',...
'TooltipString','Start Euler pole Latitude',...
'Tag','edit_pLat_ini');

h12 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@compute_euler_uicallback,h1,'edit_pAng_ini_Callback'},...
'Position',[200 201 61 21],...
'Style','edit',...
'TooltipString','Start Euler pole Angle',...
'Tag','edit_pAng_ini');

h13 = uicontrol('Parent',h1,...
'Callback',{@compute_euler_uicallback,h1,'pushbutton_polesList_Callback'},...
'Position',[289 200 121 23],...
'String','Poles selector',...
'TooltipString','Select a pole from the default list',...
'Tag','pushbutton_polesList');

h14 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@compute_euler_uicallback,h1,'edit_LonRange_Callback'},...
'Position',[104 125 41 21],...
'String','30',...
'Style','edit',...
'TooltipString','The pole will be searched arround it''s starting longitude +/- half this range',...
'Tag','edit_LonRange');

h15 = uicontrol('Parent',h1,...
'HorizontalAlignment','left',...
'Position',[20 128 84 16],...
'String','Longitude Range',...
'Style','text',...
'Tag','text1');

h16 = uicontrol('Parent',h1,...
'Position',[162 118 70 15],...
'String','N of Intervals',...
'Style','text',...
'Tag','text2');

h17 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@compute_euler_uicallback,h1,'edit_LatRange_Callback'},...
'Position',[104 98 41 21],...
'String','30',...
'Style','edit',...
'TooltipString','The pole will be searched arround it''s starting latitude +/- half this range',...
'Tag','edit_LatRange');

h18 = uicontrol('Parent',h1,...
'HorizontalAlignment','left',...
'Position',[22 99 75 18],...
'String','Latitude Range',...
'Style','text',...
'Tag','text3');

h19 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@compute_euler_uicallback,h1,'edit_AngRange_Callback'},...
'Position',[104 70 41 21],...
'String','4',...
'Style','edit',...
'TooltipString','The pole will be searched arround it''s starting angle +/- half this range',...
'Tag','edit_AngRange');

h20 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@compute_euler_uicallback,h1,'edit_Nintervals_Callback'},...
'Position',[171 95 41 21],...
'String','20',...
'Style','edit',...
'TooltipString','The range parameters are divided into this number of intervals steps',...
'Tag','edit_Nintervals');

h21 = uicontrol('Parent',h1,...
'HorizontalAlignment','left',...
'Position',[21 74 72 15],...
'String','Angular Range',...
'Style','text',...
'Tag','text4');

h22 = uicontrol('Parent',h1,...
'FontSize',10,...
'FontWeight','bold',...
'Position',[136 282 41 16],...
'String','OR',...
'Style','text',...
'Tag','text9');

h23 = uicontrol('Parent',h1,...
'FontSize',10,...
'Position',[97 339 81 16],...
'String','First Line',...
'Style','text',...
'Tag','text10');

h24 = uicontrol('Parent',h1,...
'FontSize',10,...
'Position',[348 339 81 16],...
'String','Second Line',...
'Style','text',...
'Tag','text11');

h25 = uicontrol('Parent',h1,...
'FontSize',10,...
'Position',[248 237 154 16],...
'String','Starting Pole Section',...
'Style','text',...
'Tag','text13');

h26 = uicontrol('Parent',h1,...
'FontSize',10,...
'Position',[247 357 113 16],...
'String','Data Section',...
'Style','text',...
'Tag','text14');

h27 = uicontrol('Parent',h1,...
'Position',[24 225 51 15],...
'String','Longitude',...
'Style','text',...
'Tag','text15');

h28 = uicontrol('Parent',h1,...
'Position',[116 225 51 15],...
'String','Latitude',...
'Style','text',...
'Tag','text16');

h29 = uicontrol('Parent',h1,...
'Position',[205 225 51 15],...
'String','Angle',...
'Style','text',...
'Tag','text17');

h30 = uicontrol('Parent',h1,...
'FontSize',10,...
'Position',[248 150 122 16],...
'String','Compute Section',...
'Style','text',...
'Tag','text18');

h31 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Position',[352 125 61 21],...
'Style','edit',...
'TooltipString','Computed Euler pole Longitude',...
'Tag','edit_pLon_fim');

h32 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Position',[352 97 61 21],...
'Style','edit',...
'TooltipString','Computed Euler pole Latitude',...
'Tag','edit_pLat_fim');

h33 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Position',[352 69 61 21],...
'Style','edit',...
'TooltipString','Computed Euler pole Angle',...
'Tag','edit_pAng_fim');

h34 = uicontrol('Parent',h1,...
'HorizontalAlignment','left',...
'Position',[298 128 51 15],...
'String','Longitude',...
'Style','text',...
'Tag','text19');

h35 = uicontrol('Parent',h1,...
'HorizontalAlignment','left',...
'Position',[298 99 47 17],...
'String','Latitude',...
'Style','text',...
'Tag','text20');

h36 = uicontrol('Parent',h1,'HorizontalAlignment','left',...
'Position',[298 73 46 15],'String','Angle','Style','text','Tag','text21');

h37 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Position',[431 69 61 21],...
'Style','edit','TooltipString','Residue of the merit function',...
'Tag','edit_BFresidue');

h38 = uicontrol('Parent',h1,...
'Position',[433 92 57 15],'String','BF Residue','Style','text','Tag','text22');

h39 = uicontrol('Parent',h1,...
'Callback',{@compute_euler_uicallback,h1,'pushbutton_cancel_Callback'},...
'Position',[335 10 76 23],...
'String','Cancel','Tag','pushbutton_cancel');

h40 = uicontrol('Parent',h1,...
'Callback',{@compute_euler_uicallback,h1,'pushbutton_compute_Callback'},...
'Position',[435 10 76 23],...
'String','Compute','Tag','pushbutton_compute');

h41 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Position',[430 113 61 21],...
'Style','edit',...
'TooltipString','Starting residue (starting pole) of the merit function',...
'Tag','edit_InitialResidue');


h42 = uicontrol('Parent',h1,...
'Position',[432 136 57 15],...
'String','St Residue',...
'Style','text',...
'Tag','text23');

h43 = uicontrol('Parent',h1,...
'BackgroundColor',[0.899999976158142 0.899999976158142 0.899999976158142],...
'Callback',{@compute_euler_uicallback,h1,'slider_wait_Callback'},...
'Enable','inactive',...
'Position',[10 14 231 16],...
'String',{  '' },...
'Style','slider',...
'Tag','slider_wait');


h44 = uicontrol('Parent',h1,...
'Callback',{@compute_euler_uicallback,h1,'pushbutton_stop_Callback'},...
'Position',[241 12 37 19],...
'String','STOP','Tag','pushbutton_stop');

function compute_euler_uicallback(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
feval(callback_name,hObject,[],guidata(h1));

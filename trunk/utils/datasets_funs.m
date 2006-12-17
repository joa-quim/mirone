function  datasets_funs(opt,varargin)
% This contains the Mirone's 'Datasets' funtions (except, for the time being, the 'pscoast')

switch opt(1:3)
    case 'Hot'
        DatasetsHotspots(varargin{:})
    case 'Vol'
        DatasetsVolcanoes(varargin{:})
    case 'Tid'
        DatasetsTides(varargin{:})
    case 'Iso'
        DatasetsIsochrons(varargin{:})
    case 'Pla'
        DatasetsPlateBound_PB_All(varargin{:})
    case 'Cit'
        DatasetsCities(varargin{:})
    case 'ODP'
        DatasetsODP_DSDP(varargin{:})
end

% --------------------------------------------------------------------
function DatasetsHotspots(handles)
% Read hotspot.dat which has 4 columns (lon lat name age)
if (aux_funs('msg_dlg',3,handles));     return;      end    % Test geog & no_file
fid = fopen([handles.path_data 'hotspots.dat'],'r');
tline = fgetl(fid);             % Jump the header line
todos = fread(fid,'*char');     fclose(fid);
[hot.x hot.y hot.name hot.age] = strread(todos,'%f %f %s %f');     % Note: hot.name is a cell array of chars
clear todos;

% Get rid of Fogspots that are outside the map limits
[x,y,indx,indy] = aux_funs('in_map_region',handles,hot.x,hot.y,0,[]);
hot.name(indx) = [];   hot.age(indx) = [];
hot.name(indy) = [];   hot.age(indy) = [];
n_hot = length(x);    h_hotspot = zeros(1,n_hot);
for (i = 1:n_hot)
    h_hotspot(i) = line(x(i),y(i),'Marker','p','MarkerFaceColor','r',...
        'MarkerEdgeColor','k','MarkerSize',10,'Tag','hotspot','Userdata',i);
end
draw_funs(h_hotspot,'hotspot',hot)

% --------------------------------------------------------------------
function DatasetsVolcanoes(handles)
% Read volcanoes.dat which has 6 columns (lat lon name ...)
if (aux_funs('msg_dlg',3,handles));     return;      end    % Test geog & no_file
fid = fopen([handles.path_data 'volcanoes.dat'],'r');
todos = fread(fid,'*char');
[volc.y volc.x volc.name region volc.desc volc.dating] = strread(todos,'%f %f %s %s %s %s');
fclose(fid);    clear region todos

% Get rid of Volcanoes that are outside the map limits
[x,y,indx,indy] = aux_funs('in_map_region',handles,volc.x,volc.y,0,[]);
volc.name(indx) = [];       volc.desc(indx) = [];       volc.dating(indx) = [];
volc.name(indy) = [];       volc.desc(indy) = [];       volc.dating(indy) = [];
n_volc = length(x);    h_volc = zeros(1,n_volc);
for (i = 1:n_volc)
    h_volc(i) = line(x(i),y(i),'Marker','^','MarkerFaceColor','y',...
        'MarkerEdgeColor','k','MarkerSize',8,'Tag','volcano','Userdata',i);
end
draw_funs(h_volc,'volcano',volc)

% --------------------------------------------------------------------
function DatasetsTides(handles)
if (aux_funs('msg_dlg',3,handles));     return;      end    % Test geog & no_file
load([handles.path_data 't_xtide.mat']);
% Get rid of Tide stations that are outside the map limits
[x,y] = aux_funs('in_map_region',handles,xharm.longitude,xharm.latitude,0,[]);
h_tides = line(x,y,'Marker','^','MarkerFaceColor','y','MarkerEdgeColor','k','MarkerSize',6,...
    'LineStyle','none','Tag','TideStation');
draw_funs(h_tides,'TideStation',[])

% --------------------------------------------------------------------
function DatasetsIsochrons(handles, opt)
% Read multisegment isochrons.dat which has 3 columns (lat lon id)
if (aux_funs('msg_dlg',3,handles));     return;      end    % Test geog & no_file
if (nargin == 2)            % Read a ascii multi-segment with info file
    str1 = {'*.dat;*.DAT', 'Data files (*.dat,*.DAT)';'*.*', 'All Files (*.*)'};
    
    cd(handles.last_dir)
    [FileName,PathName] = uigetfile(str1,'Select File');
    if (PathName ~= 0),         handles.last_dir = PathName;    end
    pause(0.01);        cd(handles.home_dir);       % allways go home
    if (FileName == 0),     return;     end

    if (~isempty(strfind([PathName FileName],' ')))
        errordlg('If you had RTFM you should know that names (path included) with white spaces are totaly FORBIDEN here.','ERROR')
        return
    end
    guidata(handles.figure1,handles)    
    tag = 'Unnamed';        fname = [PathName FileName];
else
    tag = 'isochron';       fname = [handles.path_data 'isochrons.dat'];
end
xx = get(handles.axes1,'Xlim');           yy = get(handles.axes1,'Ylim');
set(handles.figure1,'pointer','watch')
[bin,n_column,multi_seg,n_headers] = guess_file(fname);
if (n_column == 1 && multi_seg == 0)        % Take it as a file names list
    fid = fopen(fname);
    c = char(fread(fid))';      fclose(fid);
    names = strread(c,'%s','delimiter','\n');   clear c fid;
else
    names = {fname};
end

if (nargin == 1),   ix = 2;     iy = 1;
else                ix = 1;     iy = 2;     end
tol = 0.5;

for (k=1:length(names))
    fname = names{k};
    j = strfind(fname,filesep);
    if (isempty(j)),    fname = [PathName fname];   end
    [numeric_data,multi_segs_str] = text_read(fname,NaN,NaN,'>');
	n_isoc = 0;     n_segments = length(numeric_data);
	h_isoc = ones(n_segments,1)*NaN;   % This is the maximum we can have
	n_clear = false(n_segments,1);
	for i=1:n_segments
        % Get rid of points that are outside the map limits
        [tmpx,tmpy] = aux_funs('in_map_region',handles,numeric_data{i}(:,ix),numeric_data{i}(:,iy),tol,[xx yy]);
        if (~isempty(tmpx))
            n_isoc = n_isoc + 1;
            h_isoc(i) = line(tmpx,tmpy,'Linewidth',handles.DefLineThick,'Color',handles.DefLineColor,'Tag',tag,'Userdata',n_isoc);
            setappdata(h_isoc(i),'LineInfo',multi_segs_str{i})  % To work with the sessions and will likely replace old mechansim
        else
            n_clear(i) = 1;             % Store indexes for clearing vanished segments info
        end
	end
	multi_segs_str(n_clear) = [];       % Clear the unused info
	
	ind = isnan(h_isoc);    h_isoc(ind) = [];      % Clear unused rows in h_isoc (due to over-dimensioning)
	draw_funs(h_isoc,'isochron',multi_segs_str)
end
set(handles.figure1,'pointer','arrow')

% --------------------------------------------------------------------
function DatasetsPlateBound_PB_All(handles)
% Read and plot the of the modified (by me) Peter Bird's Plate Boundaries
if (aux_funs('msg_dlg',3,handles));     return;      end    % Test geog & no_file
set(handles.figure1,'pointer','watch')
load([handles.path_data 'PB_boundaries.mat'])

% ------------------
% Get rid of boundary segments that are outside the map limits
xx = get(handles.axes1,'Xlim');      yy = get(handles.axes1,'Ylim');
tol = 0.5;
% ------------------ OTF class
n = length(OTF);    k = [];
for i = 1:n
    ind = find(OTF(i).x_otf < xx(1)-tol | OTF(i).x_otf > xx(2)+tol);
    OTF(i).x_otf(ind) = [];     OTF(i).y_otf(ind) = [];
    if isempty(OTF(i).x_otf),   k = [k i];  end         % k is a counter to erase out-of-map segments
end
OTF(k) = [];
n = length(OTF);    k = [];
for i = 1:n
    ind = find(OTF(i).y_otf < yy(1)-tol | OTF(i).y_otf > yy(2)+tol);
    OTF(i).x_otf(ind) = [];     OTF(i).y_otf(ind) = [];
    if isempty(OTF(i).x_otf),   k = [k i];  end
end
OTF(k) = [];
% ------------------ OSR class
n = length(OSR);    k = [];
for i = 1:n
    ind = find(OSR(i).x_osr < xx(1)-tol | OSR(i).x_osr > xx(2)+tol);
    OSR(i).x_osr(ind) = [];     OSR(i).y_osr(ind) = [];
    if isempty(OSR(i).x_osr),   k = [k i];  end
end;    OSR(k) = [];
n = length(OSR);    k = [];
for i = 1:n
    ind = find(OSR(i).y_osr < yy(1)-tol | OSR(i).y_osr > yy(2)+tol);
    OSR(i).x_osr(ind) = [];     OSR(i).y_osr(ind) = [];
    if isempty(OSR(i).x_osr),   k = [k i];  end
end
OSR(k) = [];
% ------------------ CRB class
n = length(CRB);    k = [];
for i = 1:n
    ind = find(CRB(i).x_crb < xx(1)-tol | CRB(i).x_crb > xx(2)+tol);
    CRB(i).x_crb(ind) = [];     CRB(i).y_crb(ind) = [];
    if isempty(CRB(i).x_crb),   k = [k i];  end
end
CRB(k) = [];
n = length(CRB);    k = [];
for i = 1:n
    ind = find(CRB(i).y_crb < yy(1)-tol | CRB(i).y_crb > yy(2)+tol);
    CRB(i).x_crb(ind) = [];     CRB(i).y_crb(ind) = [];
    if isempty(CRB(i).x_crb),   k = [k i];  end
end
CRB(k) = [];
% ------------------ CTF class
n = length(CTF);    k = [];
for i = 1:n
    ind = find(CTF(i).x_ctf < xx(1)-tol | CTF(i).x_ctf > xx(2)+tol);
    CTF(i).x_ctf(ind) = [];     CTF(i).y_ctf(ind) = [];
    if isempty(CTF(i).x_ctf),   k = [k i];  end
end
CTF(k) = [];
n = length(CTF);    k = [];
for i = 1:n
    ind = find(CTF(i).y_ctf < yy(1)-tol | CTF(i).y_ctf > yy(2)+tol);
    CTF(i).x_ctf(ind) = [];     CTF(i).y_ctf(ind) = [];
    if isempty(CTF(i).x_ctf),   k = [k i];  end
end
CTF(k) = [];
% ------------------ CCB class
n = length(CCB);    k = [];
for i = 1:n
    ind = find(CCB(i).x_ccb < xx(1)-tol | CCB(i).x_ccb > xx(2)+tol);
    CCB(i).x_ccb(ind) = [];     CCB(i).y_ccb(ind) = [];
    if isempty(CCB(i).x_ccb),   k = [k i];  end
end
CCB(k) = [];
n = length(CCB);    k = [];
for i = 1:n
    ind = find(CCB(i).y_ccb < yy(1)-tol | CCB(i).y_ccb > yy(2)+tol);
    CCB(i).x_ccb(ind) = [];     CCB(i).y_ccb(ind) = [];
    if isempty(CCB(i).x_ccb),   k = [k i];  end
end
CCB(k) = [];
% ------------------ OCB class
n = length(OCB);    k = [];
for i = 1:n
    ind = find(OCB(i).x_ocb < xx(1)-tol | OCB(i).x_ocb > xx(2)+tol);
    OCB(i).x_ocb(ind) = [];     OCB(i).y_ocb(ind) = [];
    if isempty(OCB(i).x_ocb),   k = [k i];  end
end
OCB(k) = [];
n = length(OCB);    k = [];
for i = 1:n
    ind = find(OCB(i).y_ocb < yy(1)-tol | OCB(i).y_ocb > yy(2)+tol);
    OCB(i).x_ocb(ind) = [];     OCB(i).y_ocb(ind) = [];
    if isempty(OCB(i).x_ocb),   k = [k i];  end
end
OCB(k) = [];
% ------------------ SUB class
n = length(SUB);    k = [];
for i = 1:n
    ind = find(SUB(i).x_sub < xx(1)-tol | SUB(i).x_sub > xx(2)+tol);
    SUB(i).x_sub(ind) = [];     SUB(i).y_sub(ind) = [];
    if isempty(SUB(i).x_sub),   k = [k i];  end
end
SUB(k) = [];
n = length(SUB);    k = [];
for i = 1:n
    ind = find(SUB(i).y_sub < yy(1)-tol | SUB(i).y_sub > yy(2)+tol);
    SUB(i).x_sub(ind) = [];     SUB(i).y_sub(ind) = [];
    if isempty(SUB(i).x_sub),   k = [k i];  end
end
SUB(k) = [];

% ------------------ Finally do the ploting ------------------------------------
% Plot the OSR class
n = length(OSR);    h_PB_All_OSR = zeros(n,1);
for i = 1:n
    line(OSR(i).x_osr,OSR(i).y_osr,'Linewidth',3,'Color','k','Tag','PB_All','Userdata',i);
    h_PB_All_OSR(i) = line(OSR(i).x_osr,OSR(i).y_osr,'Linewidth',2,'Color','r','Tag','PB_All','Userdata',i);
end
% Plot the OTF class
n = length(OTF);    h_PB_All_OTF = zeros(n,1);
for i = 1:n
    line(OTF(i).x_otf,OTF(i).y_otf,'Linewidth',3,'Color','k','Tag','PB_All','Userdata',i);
    h_PB_All_OTF(i) = line(OTF(i).x_otf,OTF(i).y_otf,'Linewidth',2,'Color','g','Tag','PB_All','Userdata',i);
end
% Plot the CRB class
n = length(CRB);    h_PB_All_CRB = zeros(n,1);
for i = 1:n
    line(CRB(i).x_crb,CRB(i).y_crb,'Linewidth',3,'Color','k','Tag','PB_All','Userdata',i);
    h_PB_All_CRB(i) = line(CRB(i).x_crb,CRB(i).y_crb,'Linewidth',2,'Color','b','Tag','PB_All','Userdata',i);
end
% Plot the CTF class
n = length(CTF);    h_PB_All_CTF = zeros(n,1);
for i = 1:n
    line(CTF(i).x_ctf,CTF(i).y_ctf,'Linewidth',3,'Color','k','Tag','PB_All','Userdata',i);
    h_PB_All_CTF(i) = line(CTF(i).x_ctf,CTF(i).y_ctf,'Linewidth',2,'Color','y','Tag','PB_All','Userdata',i);
end
% Plot the CCB class
n = length(CCB);    h_PB_All_CCB = zeros(n,1);
for i = 1:n
    line(CCB(i).x_ccb,CCB(i).y_ccb,'Linewidth',3,'Color','k','Tag','PB_All','Userdata',i);
    h_PB_All_CCB(i) = line(CCB(i).x_ccb,CCB(i).y_ccb,'Linewidth',2,'Color','m','Tag','PB_All','Userdata',i);
end
% Plot the OCB class
n = length(OCB);    h_PB_All_OCB = zeros(n,1);
for i = 1:n
    line(OCB(i).x_ocb,OCB(i).y_ocb,'Linewidth',3,'Color','k','Tag','PB_All','Userdata',i);
    h_PB_All_OCB(i) = line(OCB(i).x_ocb,OCB(i).y_ocb,'Linewidth',2,'Color','c','Tag','PB_All','Userdata',i);
end
% Plot the SUB class
n = length(SUB);    h_PB_All_SUB = zeros(n,1);
for i = 1:n
    line(SUB(i).x_sub,SUB(i).y_sub,'Linewidth',3,'Color','k','Tag','PB_All','Userdata',i);
    h_PB_All_SUB(i) = line(SUB(i).x_sub,SUB(i).y_sub,'Linewidth',2,'Color','c','Tag','PB_All','Userdata',i);
end

% Join all line handles into a single variable
h.OSR = h_PB_All_OSR;    h.OTF = h_PB_All_OTF;    h.CRB = h_PB_All_CRB;    h.CTF = h_PB_All_CTF;
h.CCB = h_PB_All_CCB;    h.OCB = h_PB_All_OCB;    h.SUB = h_PB_All_SUB;
% Join all data into a single variable
data.OSR = OSR;    data.OTF = OTF;    data.CRB = CRB;    data.CTF = CTF;
data.CCB = CCB;    data.OCB = OCB;    data.SUB = SUB;
draw_funs(h,'PlateBound_All_PB',data);      set(handles.figure1,'pointer','arrow')

% --------------------------------------------------------------------
function DatasetsCities(handles,opt)
if (aux_funs('msg_dlg',3,handles));     return;      end    % Test geog & no_file
if strcmp(opt,'major')
    fid = fopen([handles.path_data 'wcity_major.dat'],'r');
    tag = 'City_major';
elseif strcmp(opt,'other')
    fid = fopen([handles.path_data 'wcity.dat'],'r');
    tag = 'City_other';
end
todos = fread(fid,'*char');     fclose(fid);
[city.x city.y city.name] = strread(todos,'%f %f %s');     % Note: city.name is a cell array of chars
% Get rid of Cities that are outside the map limits
[x,y,indx,indy] = aux_funs('in_map_region',handles,city.x,city.y,0,[]);
city.name(indx) = [];       city.name(indy) = [];
n_city = length(x);

if (n_city == 0),   return;     end     % No cities inside area. Return.
h_city = line(x,y,'LineStyle','none','Marker','o','MarkerFaceColor','k',...
    'MarkerEdgeColor','w','MarkerSize',6,'Tag',tag);
draw_funs(h_city,'DrawSymbol')                  % Set symbol's uicontextmenu

% Estimate the text position shift in order that it doesn't fall over the city symbol 
pos = get(handles.figure1,'Position');
x_lim = get(handles.axes1,'xlim');
z1 = 7 / pos(3);
dx = z1 * (x_lim(2) - x_lim(1));

city.name = strrep(city.name,'_',' ');          % Replace '_' by ' '
textHand = zeros(1,n_city);
for i = 1:n_city                                % Plot the City names
    textHand(i) = text(x(i)+dx,y(i),0,city.name{i},'Tag',tag);
    draw_funs(textHand(i),'DrawText')           % Set text's uicontextmenu
end

% --------------------------------------------------------------------
function DatasetsODP_DSDP(handles,opt)
if (aux_funs('msg_dlg',3,handles));     return;      end    % Test geog & no_file
set(handles.figure1,'pointer','watch')
fid = fopen([handles.path_data 'DSDP_ODP.dat'],'r');
todos = fread(fid,'*char');
[ODP.x ODP.y zz ODP.leg ODP.site ODP.z ODP.penetration] = strread(todos,'%f %f %s %s %s %s %s');
fclose(fid);    clear todos zz

% Get rid of Sites that are outside the map limits
[ODP.x,ODP.y,indx,indy] = aux_funs('in_map_region',handles,ODP.x,ODP.y,0,[]);

ODP.leg(indx) = [];     ODP.site(indx) = [];    ODP.z(indx) = [];   ODP.penetration(indx) = [];
ODP.leg(indy) = [];     ODP.site(indy) = [];    ODP.z(indy) = [];   ODP.penetration(indy) = [];

% If there no sites left, return
if isempty(ODP.x)
    set(handles.figure1,'pointer','arrow');    msgbox('Warning: There are no sites inside this area.','Warning');    return;
end

% Find where in file is the separation of DSDP from ODP legs
ind = find(str2double(ODP.leg) >= 100);
if ~isempty(ind),   ind = ind(1);   end
if (strcmp(opt,'ODP'))      % If only ODP sites were asked remove DSDP from data structure
    ODP.x(1:ind-1) = [];    ODP.y(1:ind-1) = [];    ODP.z(1:ind-1) = [];
    ODP.leg(1:ind-1) = [];  ODP.site(1:ind-1) = []; ODP.penetration(1:ind-1) = [];
elseif (strcmp(opt,'DSDP'))
    ODP.x(ind:end) = [];    ODP.y(ind:end) = [];    ODP.z(ind:end) = [];
    ODP.leg(ind:end) = [];  ODP.site(ind:end) = []; ODP.penetration(ind:end) = [];
end

n_sites = length(ODP.x);    h_sites = zeros(n_sites,1);
if (strcmp(opt,'DSDP'))
    if (n_sites == 0)           % If there are no sites, give a warning and exit
        set(handles.figure1,'pointer','arrow');        msgbox('Warning: There are no DSDP sites inside this area.','Warning');    return;
    end
    for i = 1:n_sites
        h_sites(i) = line(ODP.x(i),ODP.y(i),'Marker','o','MarkerFaceColor','g',...
            'MarkerEdgeColor','k','MarkerSize',8,'Tag','DSDP','Userdata',i);
    end
    draw_funs(h_sites,'ODP',ODP)
elseif (strcmp(opt,'ODP'))
    if (n_sites == 0)           % If there are no sites, give a warning and exit
        set(handles.figure1,'pointer','arrow');        msgbox('Warning: There are no ODP sites inside this area.','Warning');    return;
    end
    for i = 1:n_sites
        h_sites(i) = line(ODP.x(i),ODP.y(i),'Marker','o','MarkerFaceColor','r',...
            'MarkerEdgeColor','k','MarkerSize',8,'Tag','ODP','Userdata',i);
    end
    draw_funs(h_sites,'ODP',ODP)
else
    h_sites = zeros(length(1:ind-1),1);
    for i = 1:ind-1
        h_sites(i) = line(ODP.x(i),ODP.y(i),'Marker','o','MarkerFaceColor','g',...
            'MarkerEdgeColor','k','MarkerSize',8,'Tag','DSDP','Userdata',i);
    end
    draw_funs(h_sites,'ODP',ODP)
    h_sites = zeros(length(ind:n_sites),1);
    for (i = 1:length(ind:n_sites))
        j = i + ind - 1;
        h_sites(i) = line(ODP.x(j),ODP.y(j),'Marker','o','MarkerFaceColor','r',...
            'MarkerEdgeColor','k','MarkerSize',8,'Tag','ODP','Userdata',j);
    end
    draw_funs(h_sites,'ODP',ODP)
end
set(handles.figure1,'pointer','arrow')

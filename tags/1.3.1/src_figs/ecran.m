function varargout = ecran(varargin)
% M-File changed by desGUIDE 
% varargin   command line arguments to ecran (see VARARGIN)

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

%#function select_cols

hObject = figure('Tag','figure1','Visible','off');
handles = guihandles(hObject);
guidata(hObject, handles);
ecran_LayoutFcn(hObject,handles);
handles = guihandles(hObject);

global x_fname_in x_fname_out nx_col nx_head x_toggle home_dir
if isempty(home_dir),   handles.d_path = [pwd filesep 'data' filesep];
else                    handles.d_path = [home_dir filesep 'data' filesep];   end

% Reposition the window on screen
movegui(hObject,'east');
set(hObject,'RendererMode','auto')

handles.n_plot = 0;         % Counter of the number of lines. Used for line color painting

if isempty(nx_head) ,   handles.n_head = 0;             % When called with no file name
else                    handles.n_head = nx_head;       end

if (~nargin),   varargin(1) = {[]};   end
for i=nargin:9                         % So that varargin{1:9} allways exists.
    varargin{i+1} = [];
end

if  (strcmp(varargin{1},'reuse') && nargin < 3)
    errordlg('Error calling ecran: Minimum arguments are "type",X,Y','Error');      return
end
if ~isempty(varargin{1}) && ~ischar(varargin{1})
    errordlg('Error calling ecran: First arguments must be a string.','Error');     return
end
if (strcmp(varargin{1},'Image') || strcmp(varargin{1},'grdtrack')) && nargin < 5
    errordlg('Error calling ecran: Minimum arguments are "type",X,Y,Z','Error');    return
end

% Choose the default ploting mode
if isempty(varargin{1})          % When the file will be read latter
    set(handles.checkbox_geog,'Visible','off')          % Hide those
    set(handles.popup_selectPlot,'Visible','off');    set(handles.popup_selectSave,'Visible','off');
    set(hObject,'Name','XY view')
elseif strcmp(varargin{1},'Image')
    handles.data(:,1) = varargin{2};    handles.data(:,2) = varargin{3};    handles.data(:,3) = varargin{4};
    set(handles.popup_selectPlot,'String','Distance along profile (data units)');
    set(handles.popup_selectSave,'String',{'Save Profile on disk';'Distance,Z (data units -> ascii)';
        'Distance,Z (data units -> binary)';'X,Y,Z (data units -> ascii)';'X,Y,Z (data units -> binary)';
        'Distance,Z (data units -> mat file)'});
    % Some ToolTips
    set(handles.checkbox_geog,'TooltipString',sprintf(['Check this if your data is in geographical coordinates.\n' ...
        'You will than be able to see and save the profile in km (or m) vs z.']));
    set(handles.popup_selectPlot,'TooltipString',sprintf('Select different ways of seeing the profile'));
    set(handles.popup_selectSave,'TooltipString',sprintf('Choose how to save the profile'));
    xd = diff(handles.data(:,1));              yd = diff(handles.data(:,2));
    tmp = sqrt(xd.*xd + yd.*yd);               rd = zeros(size(handles.data,1),1);
    for i=2:length(handles.data(:,1));         rd(i) = rd(i-1) + tmp(i-1);        end
    handles.dist = rd;    % This one is by default, so save it in case user wants to save it to file
    plot(rd,handles.data(:,3));      axis tight;      zoom_j on;
    set(hObject,'Name',varargin{5})
    %handles.data(:,1) = get(h,'XData')';        handles.data(:,2) = get(h,'YData')';
    
elseif strcmp(varargin{1},'grdtrack')
    set(handles.popup_selectPlot,'String','Distance along profile (data units)');
    set(handles.popup_selectSave,'String',{'Save Profile on disk';'Distance,Z (data units -> ascii)';
        'Distance,Z (data units -> binary)';'X,Y,Z (data units -> ascii)';'X,Y,Z (data units -> binary)';
        'Distance,Z (data units -> mat file)'});
    % Some ToolTips
    set(handles.checkbox_geog,'TooltipString',sprintf(['Check this if your data is in geographical coordinates.\n' ...
        'You will than be able to see and save the profile in km (or m) vs z.']));
    set(handles.popup_selectPlot,'TooltipString',sprintf('Select different ways of seeing the profile'));
    set(handles.popup_selectSave,'TooltipString',sprintf('Choose how to save the profile'));
    data = read_xy(handles,x_fname_in,nx_col);
    if isempty(x_toggle) || ~isempty(findstr(x_toggle,':o')),   handles.data(:,1) = data(:,1); handles.data(:,2) = data(:,2);
    else                                                        handles.data(:,1) = data(:,2); handles.data(:,2) = data(:,1);
    end
    handles.data(:,3) = data(:,nx_col); % nx_col because only last column contains the interpolated z value
    xd = diff(data(:,1));               yd = diff(data(:,2));
    tmp = sqrt(xd.*xd + yd.*yd);        rd(1:length(data(:,1))) = 0;
    for i=2:length(data(:,1));          rd(i) = rd(i-1) + tmp(i-1);        end
    handles.dist = rd;                  % This one is by default, so save it in case user wants to save it to file
    plot(rd,data(:,nx_col));          axis tight;     zoom_j on;
    set(hObject,'Name',varargin{5})
    
elseif strcmp(varargin{1},'reuse')                      % Case of auto-referenced call
    set(handles.checkbox_geog,'Visible','off')          % Hide those
    set(handles.popup_selectPlot,'Visible','off');    set(handles.popup_selectSave,'Visible','off');
    handles.data(:,1) = varargin{2};        handles.data(:,2) = varargin{3};
    if ~isempty(varargin{9}) && strcmp(varargin{9},'semilogy')
        semilogy(handles.data(:,1),handles.data(:,2));          axis tight;     zoom_j on;
    elseif ~isempty(varargin{9}) && strcmp(varargin{9},'semilogx')
        semilogx(handles.data(:,1),handles.data(:,2));          axis tight;     zoom_j on;
    else
        plot(handles.data(:,1),handles.data(:,2));              axis tight;     zoom_j on;
    end
    if ~isempty(varargin{5}),    set(hObject,'Name',varargin{5});    end     % Figure Name
    if ~isempty(varargin{6}),    xlabel(varargin{6});                end     % XLabel
    if ~isempty(varargin{7}),    ylabel(varargin{7});                end     % YLabel
    if ~isempty(varargin{8}),    title(varargin{8});                 end     % Title

elseif strcmp(varargin{1},'xy')                         % Case of call from a generic gmt xy program
    set(handles.checkbox_geog,'Visible','off')          % Hide those
    set(handles.popup_selectPlot,'Visible','off');    set(handles.popup_selectSave,'Visible','off');
    data = read_xy(handles,x_fname_in,nx_col);
    if isempty(x_toggle) || ~isempty(findstr(x_toggle,':o')),    h=plot(data(:,1),data(:,2));
    else                                                        h=plot(data(:,2),data(:,1));    end
    axis tight;     zoom_j on;    set(hObject,'Name',varargin{5})
    handles.data(:,1) = get(h,'XData')';        handles.data(:,2) = get(h,'YData')';
    hold on

    if exist(x_fname_out)
        data_out = read_xy(handles,x_fname_out,2);
        plot(data_out(:,1),data_out(:,2),'Color','red');     axis tight;
    end
    hold off
end

load([handles.d_path 'mirone_pref.mat']);
if (~isempty(home_dir)),            handles.home_dir = home_dir;
else                                handles.home_dir = pwd;     end
try
    if (iscell(directory_list)),    handles.last_dir = directory_list{1};
    else                            handles.last_dir = directory_list;  end
catch
    handles.last_dir = pwd;
end

handles.hFig = hObject;         handles.hAxes = gca;    % sometimes needed in 1 st or 2 nd derivatives
%age_start = dir(x_fname_out);

% Choose default command line output for ecran_export
handles.output = hObject;
guidata(hObject, handles);

set(hObject,'Visible','on');

% UIWAIT makes ecran_export wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% NOTE: If you make uiwait active you have also to uncomment the next three lines
% handles = guidata(hObject);
% out = ecran_OutputFcn(hObject, [], handles);
% varargout{1} = out;

% --- Outputs from this function are returned to the command line.
function varargout = ecran_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% Get default command line output from handles structure
varargout{1} = handles.output;
delete(handles.figure1);

function xy = read_xy(handles,file,n_col)
%global nx_head
% build the format string to read the data n_columns
format = [];    fid = fopen(file,'r');
for (i=1:n_col),    format = [format '%f '];    end

% Jump header lines
for (i = 1:handles.n_head),    tline = fgetl(fid);  end

todos = fread(fid,'*char');
xy = sscanf(todos,format,[n_col inf])';    % After hours strugling agains this FILHO DA PUTA, I might have found
fclose(fid);
%if (n_col > 2)    xy = xy(:,1:2);   end

% --- Executes on button press in checkbox_geog.
function checkbox_geog_Callback(hObject, eventdata, handles)
if get(hObject,'Value')
    set(handles.popup_selectPlot,'String',{'Distance along profile (data units)';
    'Distance along profile (km)';'Distance along profile (NM)'});
    set(handles.popup_selectSave,'String',{'Save Profile on disk';'Distance,Z (data units -> ascii)';
        'Distance,Z (data units -> binary)';'Distance,Z (km -> ascii)';'Distance,Z (km -> binary)';
        'Distance,Z (NM -> ascii)';'Distance,Z (NM -> binary)';
        'X,Y,Z (data units -> ascii)';'X,Y,Z (data units -> binary)';
        'Distance,Z (NM -> mat file)';'Distance,Z (km -> mat file)';'Distance,Z (data units -> mat file)'});
else
    set(handles.popup_selectPlot,'String','Distance along profile (data units)','Value',1);
    set(handles.popup_selectSave,'String',{'Save Profile on disk';'Distance,Z (data units -> ascii)';
        'Distance,Z (data units -> binary)';'X,Y,Z (data units -> ascii)';'X,Y,Z (data units -> binary)';
        'Distance,Z (data units -> mat file)'},'Value',1);
end

% ---------------------------------------------------------------------------------
function popup_selectPlot_Callback(hObject, eventdata, handles)
val = get(hObject,'Value');     str = get(hObject, 'String');
D2R = pi/180;
h = findobj(handles.axes1,'Type','line');
switch str{val};
    case 'Distance along profile (data units)'  % Compute the accumulated distance along profile in data units
        xd = diff(handles.data(:,1));               yd = diff(handles.data(:,2));
        tmp = sqrt(xd.*xd + yd.*yd);                rd = zeros(size(handles.data,1),1);
        for i=2:length(handles.data(:,1));          rd(i) = rd(i-1) + tmp(i-1);        end
        %plot(rd,handles.data(:,3));                 axis tight;
        set(h,'XData',rd);                          axis tight;
    case 'Distance along profile (km)'           % Compute the accumulated distance along profile in km
        deg2km = 111.1949;
        xd = diff(handles.data(:,1)*deg2km).*cos(handles.data(2:end,2)*D2R);     yd = diff(handles.data(:,2))*deg2km;
        tmp = sqrt(xd.*xd + yd.*yd);                rd = zeros(size(handles.data,1),1);
        for i=2:length(handles.data(:,1));          rd(i) = rd(i-1) + tmp(i-1);        end
        %plot(rd,handles.data(:,3));                 axis tight;
        set(h,'XData',rd);                          axis tight;
    case 'Distance along profile (NM)'            % Compute the accumulated distance along profile in Nmiles
        deg2nm = 60.04;     %deg2km = 111194.9;
        xd = diff(handles.data(:,1)*deg2nm).*cos(handles.data(2:end,2)*D2R);     yd = diff(handles.data(:,2))*deg2nm;
        tmp = sqrt(xd.*xd + yd.*yd);                rd = zeros(size(handles.data,1),1);
        for i=2:length(handles.data(:,1)),          rd(i) = rd(i-1) + tmp(i-1);        end
        %plot(rd,handles.data(:,3));                 axis tight;
        set(h,'XData',rd);                          axis tight;
end
guidata(hObject, handles);

% ---------------------------------------------------------------------------------
function popup_selectSave_Callback(hObject, eventdata, handles)
val = get(hObject,'Value');     str = get(hObject, 'String');
D2R = pi/180;
switch str{val};
    case 'Save Profile on disk'                    %
    case 'Distance,Z (data units -> ascii)'                    % Save profile in ascii data units
        cd(handles.last_dir);
        [FileName,PathName] = uiputfile({'*.dat', 'Dist Z (*.dat)';'*.*', 'All Files (*.*)'}, 'Distance,Z (ascii)');
        pause(0.01);    cd(handles.home_dir);
        if isequal(FileName,0);     set(hObject,'Value',1);     return;      end     % User gave up
        double2ascii([PathName FileName],[handles.dist' handles.data(:,3)'],'%f\t%f');
        set(hObject,'Value',1);
    case 'Distance,Z (data units -> binary)'                    % Save profile in binary data units
        cd(handles.last_dir);
        [FileName,PathName] = uiputfile({'*.dat', 'Dist Z (*.dat)';'*.*', 'All Files (*.*)'}, 'Distance,Z (binary float)');
        pause(0.01);    cd(handles.home_dir);
        if isequal(FileName,0);     set(hObject,'Value',1);     return;      end     % User gave up
        fid = fopen([PathName FileName],'wb');
        fwrite(fid,[handles.dist' handles.data(:,3)']','float');  fclose(fid);
        set(hObject,'Value',1);
    case 'Distance,Z (km -> ascii)'                    % Save profile in ascii (km Z) 
        cd(handles.last_dir);
        [FileName,PathName] = uiputfile({'*.dat', 'Dist Z (*.dat)';'*.*', 'All Files (*.*)'}, 'Distance (km),Z (ascii)');
        pause(0.01);    cd(handles.home_dir);
        if isequal(FileName,0);     set(hObject,'Value',1);     return;      end     % User gave up
        deg2km = 111.1949;
        xd = diff(handles.data(:,1)*deg2km).*cos(handles.data(2:end,2)*D2R);     yd = diff(handles.data(:,2))*deg2km;
        tmp = sqrt(xd.*xd + yd.*yd);                rd = zeros(size(handles.data,1),1);
        for i=2:length(handles.data(:,1));          rd(i) = rd(i-1) + tmp(i-1);        end
        double2ascii([PathName FileName],[rd handles.data(:,3)],'%f\t%f')
        set(hObject,'Value',1);
    case 'Distance,Z (km -> binary)'                    % Save profile in binary (km Z) 
        cd(handles.last_dir);
        [FileName,PathName] = uiputfile({'*.dat', 'Dist Z (*.dat)';'*.*', 'All Files (*.*)'}, 'Distance (km),Z (binary float)');
        pause(0.01);    cd(handles.home_dir);
        if isequal(FileName,0);     set(hObject,'Value',1);     return;      end     % User gave up
        deg2km = 111.1949;
        xd = diff(handles.data(:,1)*deg2km).*cos(handles.data(2:end,2)*D2R);     yd = diff(handles.data(:,2))*deg2km;
        tmp = sqrt(xd.*xd + yd.*yd);                rd = zeros(1,size(handles.data,1));
        for i=2:length(handles.data(:,1));          rd(i) = rd(i-1) + tmp(i-1);        end
        fid = fopen([PathName FileName],'wb');
        fwrite(fid,[rd' handles.data(:,3)']','float');  fclose(fid);
        set(hObject,'Value',1);
    case 'Distance,Z (NM -> ascii)'                    % Save profile in ascii (NM Z) 
        cd(handles.last_dir);
        [FileName,PathName] = uiputfile({'*.dat', 'Dist Z (*.dat)';'*.*', 'All Files (*.*)'}, 'Distance (m),Z (ascii)');
        pause(0.01);    cd(handles.home_dir);
        if isequal(FileName,0);     set(hObject,'Value',1);     return;      end     % User gave up
        deg2nm = 60.04;
        xd = diff(handles.data(:,1)*deg2nm).*cos(handles.data(2:end,2)*D2R);     yd = diff(handles.data(:,2))*deg2nm;
        tmp = sqrt(xd.*xd + yd.*yd);                rd = zeros(size(handles.data,1),1);
        for i=2:length(handles.data(:,1));          rd(i) = rd(i-1) + tmp(i-1);        end
        double2ascii([PathName FileName],[rd handles.data(:,3)],'%f\t%f')
        set(hObject,'Value',1);
    case 'Distance,Z (NM -> binary)'                    % Save profile in binary (NM Z) 
        cd(handles.last_dir);
        [FileName,PathName] = uiputfile({'*.dat', 'Dist Z (*.dat)';'*.*', 'All Files (*.*)'}, 'Distance (m),Z (binary float)');
        pause(0.01);    cd(handles.home_dir);
        if isequal(FileName,0);     set(hObject,'Value',1);     return;      end     % User gave up
        deg2nm = 60.04;
        xd = diff(handles.data(:,1)*deg2nm).*cos(handles.data(2:end,2)*D2R);     yd = diff(handles.data(:,2))*deg2nm;
        tmp = sqrt(xd.*xd + yd.*yd);                rd = zeros(1,size(handles.data,1));
        for i=2:length(handles.data(:,1));          rd(i) = rd(i-1) + tmp(i-1);        end
        fid = fopen([PathName FileName],'wb');
        fwrite(fid,[rd' handles.data(:,3)']','float');  fclose(fid);
        set(hObject,'Value',1);
    case 'X,Y,Z (data units -> ascii)'                    % Save profile in ascii (km Z) 
        cd(handles.last_dir);
        [FileName,PathName] = uiputfile({'*.dat', 'x,y,z (*.dat)';'*.*', 'All Files (*.*)'}, 'X,Y,Z (ascii)');
        pause(0.01);    cd(handles.home_dir);
        if isequal(FileName,0);     set(hObject,'Value',1);     return;      end     % User gave up
        double2ascii([PathName FileName],[handles.data(:,1) handles.data(:,2) handles.data(:,3)],'%f\t%f\t%f')
        set(hObject,'Value',1);
    case 'X,Y,Z (data units -> binary)'                    % Save profile in binary (km Z) 
        cd(handles.last_dir);
        [FileName,PathName] = uiputfile({'*.dat', 'x,y,z (*.dat)';'*.*', 'All Files (*.*)'}, 'X,Y,Z (binary float)');
        pause(0.01);    cd(handles.home_dir);
        if isequal(FileName,0);     set(hObject,'Value',1);     return;      end     % User gave up
        fid = fopen([PathName FileName],'wb');
        fwrite(fid,[handles.data(:,1)' handles.data(:,2)' handles.data(:,3)']','float');  fclose(fid);
        set(hObject,'Value',1);
    case 'Distance,Z (NM -> mat file)'                    % Save profile in mat file (m Z) 
        cd(handles.last_dir);
        [FileName,PathName] = uiputfile({'*.mat', 'Dist Z (*.mat)';'*.*', 'All Files (*.*)'}, 'Distance (m),Z (Matlab mat file)');
        pause(0.01);    cd(handles.home_dir);
        if isequal(FileName,0);     set(hObject,'Value',1);     return;      end     % User gave up
        deg2nm = 60.04;
        xd = diff(handles.data(:,1)*deg2nm).*cos(handles.data(2:end,2)*D2R);     yd = diff(handles.data(:,2))*deg2nm;
        tmp = sqrt(xd.*xd + yd.*yd);                rd = zeros(1,size(handles.data,1));
        for i=2:length(handles.data(:,1));          rd(i) = rd(i-1) + tmp(i-1);        end
        R = rd';   Z = handles.data(:,3)';          % More one BUG, handles.data(:,3) canot be saved
        save([PathName FileName],'R','Z')
        set(hObject,'Value',1);
    case 'Distance,Z (km -> mat file)'                    % Save profile in mat file (km Z)
        cd(handles.last_dir);
        [FileName,PathName] = uiputfile({'*.mat', 'Dist Z (*.mat)';'*.*', 'All Files (*.*)'}, 'Distance (km),Z (Matlab mat file))');
        pause(0.01);    cd(handles.home_dir);
        if isequal(FileName,0);     set(hObject,'Value',1);     return;      end     % User gave up
        deg2km = 111.1949;
        xd = diff(handles.data(:,1)*deg2km).*cos(handles.data(2:end,2)*D2R);     yd = diff(handles.data(:,2))*deg2km;
        tmp = sqrt(xd.*xd + yd.*yd);                rd = zeros(1,size(handles.data,1));
        for i=2:length(handles.data(:,1));          rd(i) = rd(i-1) + tmp(i-1);        end
        R = rd';   Z = handles.data(:,3)';          % More one BUG, handles.data(:,3) canot be saved
        save([PathName FileName],'R','Z')
        set(hObject,'Value',1);
    case 'Distance,Z (data units -> mat file)'                    % Save profile in binary data units
        cd(handles.last_dir);
        [FileName,PathName] = uiputfile({'*.mat', 'Dist Z (*.mat)';'*.*', 'All Files (*.*)'}, 'Distance,Z (Matlab mat file))');
        pause(0.01);    cd(handles.home_dir);
        if isequal(FileName,0);     set(hObject,'Value',1);     return;      end     % User gave up
        R = handles.dist';   Z = handles.data(:,3)';      % More one BUG, handles.data(:,3) canot be saved
        save([PathName FileName],'R','Z')
        set(hObject,'Value',1);
end

% --------------------------------------------------------------------
function menuFile_Callback(hObject, eventdata, handles)
% --------------------------------------------------------------------

% --------------------------------------------------------------------
function FileExport_Callback(hObject, eventdata, handles)
FILEMENUFCN(gcf,'FileExport')

% --------------------------------------------------------------------
function FilePrintSetup_Callback(hObject, eventdata, handles)
print -dsetup

% --------------------------------------------------------------------
function FilePrint_Callback(hObject, eventdata, handles)
if (ispc);      print -v
else            print;  end

% --------------------------------------------------------------------
function FileExit_Callback(hObject, eventdata, handles)
uiresume(handles.figure1);

% --------------------------------------------------------------------
function FileOpen_Callback(hObject, eventdata, handles)
% Read the file and select what columns to plot
cd(handles.last_dir);
[FileName,PathName] = uigetfile({'*.dat;*.DAT;', 'Data files (*.dat,*.DAT)';'*.*', 'All Files (*.*)'},'Select input data');
pause(0.01);    cd(handles.home_dir);
if isequal(FileName,0),   return;     end
fname = [PathName FileName];
handles.last_dir = PathName;

data = text_read(fname,NaN);
if (isempty(data))
    errordlg('File doesn''t have any recognized nymeric data (Quiting).','Error');    return
end

% If msgbox exist we have to move it from behind the main window. So get it's handle
hMsgFig = gcf;
if (handles.hFig ~= hMsgFig)
    uistack(hMsgFig,'top');    % If error msgbox exists, bring it forward
    % Here we have a stupid problem. If don't kill the message window before the
    % select_cols is called this later wont work. FDS I have no more patiente for this.
    pause(1)
    try    delete(hMsgFig);    end
end

out = select_cols(data,'xy',fname,1000);
if (isempty(out)) ,   return;    end

%lineHand = plot(data(:,out(1)),data(:,out(2)));    axis tight;     zoom_j on;
lineHand = line('Parent',handles.axes1,'XData',data(:,out(1)),'YData',data(:,out(2)));    axis tight;     zoom_j on;
handles.n_plot = handles.n_plot + 1;
if (handles.n_plot > 1)
    c_order = get(handles.axes1,'ColorOrder');
    if (handles.n_plot <= 7)
        nc = handles.n_plot;
    else
        nc = rem(handles.n_plot,7);     % recycle through the default colors
    end
    cor = c_order(nc,:);
    set(lineHand,'Color',cor)
end
handles.data = [data(:,out(1)) data(:,out(2))];     % NOTE, if handles.n_plot > 1 only last data is saved
guidata(hObject, handles);

% --------------------------------------------------------------------
function FileSave_Callback(hObject, eventdata, handles)
cd(handles.last_dir);
[FileName,PathName] = uiputfile({'*.dat', 'X,Y (*.dat)';'*.*', 'All Files (*.*)'}, 'X,Y (ascii)');
pause(0.01);    cd(handles.home_dir);
if isequal(FileName,0);     return;      end     % User gave up
h_lin=findobj(get(handles.hAxes,'Children'),'LineStyle','-');
xx = get(h_lin,'XData');            yy = get(h_lin,'YData');
double2ascii([PathName FileName],[xx' yy'],'%f\t%f');

% --------------------------------------------------------------------
function VoidMenuAnalysis_Callback(hObject, eventdata, handles)
% --------------------------------------------------------------------

% --------------------------------------------------------------------
function VoidAnalysisFFT_Callback(hObject, eventdata, handles)
% --------------------------------------------------------------------

% --------------------------------------------------------------------
function AnalysisFFT_AmpSpectrum_Callback(hObject, eventdata, handles)
h_lin = findobj(get(handles.hAxes,'Children'),'LineStyle','-');
x = get(h_lin,'XData');
Fs = 1 / (x(2) - x(1));         % Sampling frequency
Fn=Fs/2;                        % Nyquist frequency
x = get(h_lin,'YData');
NFFT=2.^(ceil(log(length(x))/log(2)));      % Next highest power of 2 greater than or equal to length(x)
FFTX=fft(x,NFFT);                           % Take fft, padding with zeros, length(FFTX)==NFFT
NumUniquePts = ceil((NFFT+1)/2);
FFTX=FFTX(1:NumUniquePts);                  % fft is symmetric, throw away second half
MX=abs(FFTX);                               % Take magnitude of X
% Multiply by 2 to take into account the fact that we threw out second half of FFTX above
MX=MX*2;                MX(1)=MX(1)/2;      % Account for endpoint uniqueness
MX(length(MX))=MX(length(MX))/2;            % We know NFFT is even
% Scale the FFT so that it is not a function of the length of x.
MX=MX/length(x);        f=(0:NumUniquePts-1)*2*Fn/NFFT;
ecran('reuse',f,MX,[],'Amplitude Spectrum','frequency',[],'Amplitude Spectrum')

% --------------------------------------------------------------------
function AnalysisFFT_PSD_Callback(hObject, eventdata, handles)
h_lin = findobj(get(handles.hAxes,'Children'),'LineStyle','-');
x = get(h_lin,'XData');         Fs = 1 / (x(2) - x(1));         % Sampling frequency
x = get(h_lin,'YData');
[Pxx,w] = psd(x,Fs);
% We want to guarantee that the result is an integer if X is a negative power of 10.
% To do so, we force some rounding of precision by adding 300-300.
Pxx = (10.*log10(Pxx)+300)-300;    % Compute db
ecran('reuse',w,Pxx,[],'Power Spectrum','Frequency (Hz)','Power Spectral Density (dB/Hz)','Periodogram PSD Estimate')

% --------------------------------------------------------------------
function AnalysisAutocorrelation_Callback(hObject, eventdata, handles)
h_lin = findobj(get(handles.hAxes,'Children'),'LineStyle','-');
xx = get(h_lin,'XData');        yy = get(h_lin,'YData');
c = autocorr(yy);               n = length(yy);
ecran('reuse',xx,c(n:end),[],'Normalized Autocorrelation','Lag in user X units')

% --------------------------------------------------------------------
function AnalysisRemoveMean_Callback(hObject, eventdata, handles)
h_lin = findobj(get(handles.hAxes,'Children'),'LineStyle','-');
xx = get(h_lin,'XData');        yy = get(h_lin,'YData');
ecran('reuse',xx,yy-mean(yy),[],'Mean Removed')

% --------------------------------------------------------------------
function AnalysisRemoveTrend_Callback(hObject, eventdata, handles)
h_lin = findobj(get(handles.hAxes,'Children'),'LineStyle','-');
xx = get(h_lin,'XData');        yy = get(h_lin,'YData');
p = polyfit(xx,yy,1);           y = polyval(p,xx);
ecran('reuse',xx,yy-y,[],'Trend Removed')

% --------------------------------------------------------------------
function AnalysisSmoothSpline_Callback(hObject, eventdata, handles)
h = get(gca,'Children');      xx = get(h,'XData');        yy = get(h,'YData');
[pp,p] = spl_fun('csaps',xx,yy);      % This is just to get csaps's p estimate
y = spl_fun('csaps',xx,yy,p,xx);
delete(findobj(get(handles.hAxes,'Children')));
hold on;    plot(xx,yy,'r.');   plot(xx,y);     hold off;   axis tight;     zoom_j on;
smoothing_param(p,[xx(1) xx(2)-xx(1) xx(end)],handles.hFig,handles.hAxes);
guidata(hObject, handles);

% --------------------------------------------------------------------
function Analysis1derivative_Callback(hObject, eventdata, handles)
h_lin = findobj(get(handles.hAxes,'Children'),'LineStyle','-');       % this is the one to be replaced
xx = get(h_lin,'XData');            yy = get(h_lin,'YData');
pp = spl_fun('csaps',xx,yy,1);       % Use 1 for not smoothing, just interpolate
v = spl_fun('ppual',pp,xx,'l','first');
ecran('reuse',xx,v,[],'First derivative')

% --------------------------------------------------------------------
function Analysis2derivative_Callback(hObject, eventdata, handles)
h_lin = findobj(get(handles.hAxes,'Children'),'LineStyle','-');       % this is the one to be replaced
xx = get(h_lin,'XData');            yy = get(h_lin,'YData');
pp = spl_fun('csaps',xx,yy,1);       % Use 1 for not smoothing, just interpolate
v = spl_fun('ppual',pp,xx,'l','second');
ecran('reuse',xx,v,[],'Second derivative')

% --------------------------------------------------------------------
function ac = autocorr(x)
%AUTOCORR Computes normalized auto-correlation of vector X.
[x,nshift] = shiftdim(x);
maxlag = size(x,1) - 1;
x = x(:);   m = size(x,1);
% Compute Autocorrelation via FFT
X = fft(x,2^nextpow2(2*m-1));
ac = ifft(abs(X).^2);

ac = real(ac);            % We want only the real part
% Move negative lags before positive
ac = [ac(end-maxlag+1:end,:);ac(1:maxlag+1,:)];
ac = ac./ac(maxlag+1);     % Normalize by ac[0]

% If first vector is a row, return a row
ac = shiftdim(ac,-nshift);

% --------------------------------------------------------------------
function [Pxx,w] = psd(xw,Fs)
%Power Spectral Density estimate via periodogram method.
N = length(xw);     xw = xw(:);
nfft  = max(256, 2^nextpow2(N));

nx = size(xw,2);
xw = [xw; zeros(nfft-N,1)];     % pad with zeros (I REALY don't like this)
if (nx~=1),  xw = xw.';  end;    clear nx;

% Compute the periodogram power spectrum [Power] estimate
Sxx =(abs(fft(xw)).^2)./N; 

% Generate the frequency vector in [rad/sample] at which Sxx was computed
w = 2.*pi.*(0 : 1./nfft : 1-1./nfft);

% Compute the Power/freq (PSD), the Power and the frequency at which it is computed
w = w(:);

% Generate the spectrum
if rem(nfft,2),   select = 1:(nfft+1)/2;          % odd
else                select = 1:nfft/2+1;    end     % even
Sxx_unscaled = Sxx(select);
w = w(select);
Sxx = [Sxx_unscaled(1); 2*Sxx_unscaled(2:end-1); Sxx_unscaled(end)];

Pxx = Sxx./Fs;      % Scale by the sampling frequency to obtain the psd
w = w.*Fs./(2.*pi); % Scale the frequency vector from rad/sample to Hz   

% --------------------------------------------------------------------
% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
delete(handles.figure1);

% --- Executes on key press over figure1 with no controls selected.
function figure1_KeyPressFcn(hObject, eventdata, handles)
delete(handles.figure1);

% --- Creates and returns a handle to the GUI figure. 
function ecran_LayoutFcn(h1,handles)

set(h1,'Units','centimeters',...
'PaperUnits',get(0,'defaultfigurePaperUnits'),...
'CloseRequestFcn',{@figure1_CloseRequestFcn,handles},...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'KeyPressFcn',{@figure1_KeyPressFcn,handles},...
'MenuBar','none',...
'Name','ecran',...
'NumberTitle','off',...
'PaperPosition',[1 2.5 15.625 6.953],...
'PaperSize',[20.98404194812 29.67743169791],...
'PaperType',get(0,'defaultfigurePaperType'),...
'Position',[13.7214409375 12.0822707291667 21.4414038541667 9.041874375],...
'Renderer',get(0,'defaultfigureRenderer'),...
'RendererMode','manual',...
'Tag','figure1',...
'UserData',[]);

axes('Parent',h1,...
'CameraPosition',[0.5 0.5 9.16025403784439],...
'CameraPositionMode',get(0,'defaultaxesCameraPositionMode'),...
'NextPlot','Add',...
'Color',get(0,'defaultaxesColor'),...
'ColorOrder',get(0,'defaultaxesColorOrder'),...
'Position',[0.0493218249075216 0.184210526315789 0.939580764488286 0.730994152046784],...
'XColor',get(0,'defaultaxesXColor'),...
'YColor',get(0,'defaultaxesYColor'),...
'ZColor',get(0,'defaultaxesZColor'),...
'Tag','axes1');

uicontrol('Parent',h1,...
'Units','normalized',...
'Callback',{@ecran_uicallback,h1,'checkbox_geog_Callback'},...
'Position',[0.0480887792848335 0.0360655737704918 0.23921085080148 0.0524590163934426],...
'String','Geographical coordinates',...
'Style','checkbox',...
'Tag','checkbox_geog');

uicontrol('Parent',h1,...
'Units','normalized',...
'BackgroundColor',[1 1 1],...
'Callback',{@ecran_uicallback,h1,'popup_selectPlot_Callback'},...
'Position',[0.335099337748344 0.0192926045016077 0.327152317880795 0.0803858520900321],...
'String',{  'Distance along profile (data units)'; 'Distance along profile (km)'; 'Distance along profile (NM)' },...
'Style','popupmenu',...
'Value',1,...
'Tag','popup_selectPlot');

uicontrol('Parent',h1,...
'Units','normalized',...
'BackgroundColor',[1 1 1],...
'Callback',{@ecran_uicallback,h1,'popup_selectSave_Callback'},...
'Position',[0.691390728476821 0.0192926045016077 0.288741721854305 0.0803858520900321],...
'String',{  'Save Profile on disk'; 'distance Z (data units -> ascii)'; 'distance Z (data units -> binary)'; 'distance Z (km -> ascii)'; 'distance Z (km -> binary)'; 'distance Z (NM -> ascii)'; 'distance Z (NM -> binary)'; 'X Y Z (data units -> ascii)'; 'X Y Z (data units -> binary)' },...
'Style','popupmenu',...
'Value',1,...
'Tag','popup_selectSave');

h10 = uimenu('Parent',h1,...
'Callback',{@ecran_uicallback,h1,'menuFile_Callback'},...
'Label','File',...
'Tag','menuFile');

uimenu('Parent',h10,...
'Callback',{@ecran_uicallback,h1,'FileOpen_Callback'},...
'Label','Open',...
'Tag','FileOpen');

uimenu('Parent',h10,...
'Callback',{@ecran_uicallback,h1,'FileSave_Callback'},...
'Label','Save',...
'Tag','FileSave');

uimenu('Parent',h10,...
'Callback',{@ecran_uicallback,h1,'FileExport_Callback'},...
'Label','Export...',...
'Separator','on',...
'Tag','FileExport');

uimenu('Parent',h10,...
'Callback',{@ecran_uicallback,h1,'FilePrintSetup_Callback'},...
'Label','Print Setup',...
'Separator','on',...
'Tag','FilePrintSetup');

uimenu('Parent',h10,...
'Callback',{@ecran_uicallback,h1,'FilePrint_Callback'},...
'Label','Print...',...
'Tag','FilePrint');

uimenu('Parent',h10,...
'Callback',{@ecran_uicallback,h1,'FileExit_Callback'},...
'Label','Exit',...
'Separator','on',...
'Tag','FileExit');

h17 = uimenu('Parent',h1,...
'Callback',{@ecran_uicallback,h1,'VoidMenuAnalysis_Callback'},...
'Label','Analysis',...
'Tag','VoidMenuAnalysis');

uimenu('Parent',h17,...
'Callback',{@ecran_uicallback,h1,'AnalysisRemoveMean_Callback'},...
'Label','Remove Mean',...
'Tag','AnalysisRemoveMean');

uimenu('Parent',h17,...
'Callback',{@ecran_uicallback,h1,'AnalysisRemoveTrend_Callback'},...
'Label','Remove Trend',...
'Tag','AnalysisRemoveTrend');

h20 = uimenu('Parent',h17,...
'Callback',{@ecran_uicallback,h1,'VoidAnalysisFFT_Callback'},...
'Label','FFT',...
'Separator','on',...
'Tag','VoidAnalysisFFT');

uimenu('Parent',h20,...
'Callback',{@ecran_uicallback,h1,'AnalysisFFT_AmpSpectrum_Callback'},...
'Label','Amplitude Spectrum',...
'Tag','AnalysisFFT_AmpSpectrum');

uimenu('Parent',h20,...
'Callback',{@ecran_uicallback,h1,'AnalysisFFT_PSD_Callback'},...
'Label','Power Spectrum Density',...
'Tag','AnalysisFFT_PSD');

uimenu('Parent',h17,...
'Callback',{@ecran_uicallback,h1,'AnalysisAutocorrelation_Callback'},...
'Label','Autocorrelation',...
'Tag','AnalysisAutocorrelation');

uimenu('Parent',h17,...
'Callback',{@ecran_uicallback,h1,'AnalysisSmoothSpline_Callback'},...
'Label','Smoothing Spline',...
'Separator','on',...
'Tag','AnalysisSmoothSpline');

uimenu('Parent',h17,...
'Callback',{@ecran_uicallback,h1,'Analysis1derivative_Callback'},...
'Label','1 st derivative',...
'Tag','Analysis1derivative');

uimenu('Parent',h17,...
'Callback',{@ecran_uicallback,h1,'Analysis2derivative_Callback'},...
'Label','2 nd derivative',...
'Tag','Analysis2derivative');

function ecran_uicallback(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
feval(callback_name,hObject,[],guidata(h1));

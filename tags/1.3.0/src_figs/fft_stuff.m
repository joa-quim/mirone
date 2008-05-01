function varargout = fft_stuff(varargin)
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
fft_stuff_LayoutFcn(hObject,handles);
handles = guihandles(hObject);

global home_dir
movegui(hObject,'center')
% Case when this function was called directly
if isempty(home_dir),   home_dir = pwd;     end

if isempty(home_dir)        % Case when this function was called directly
    handles.path_data = ['data' filesep];
else
    handles.path_data = [home_dir filesep 'data' filesep];
end

handles.h_calling_fig = [];     % Handles to the calling figure
handles.geog = 0;               % Set this as default
handles.Z1 = [];
handles.Z2 = [];
last_dir = [];
grid_in_continue = 0;

if (~isempty(varargin))         % When called from a Mirone window
    handles.h_calling_fig = varargin{1};
    handles.Z1 = varargin{2};
    handles.head_Z1 = varargin{3};
    handles.geog = varargin{4};
    mode = varargin{5};
    if (strcmp(mode,'Allopts'))     % "Slow mode" show all options in this figure
        grid_in_continue = 1;
    end
    [handles.orig_nrows,handles.orig_ncols] = size(handles.Z1);
    [w,nlist] = mboard([],handles.orig_ncols,handles.orig_nrows,0,0);
    handles.new_nx = w(1);
    handles.new_ny = w(2);
    rlat = (handles.head_Z1(4) + handles.head_Z1(3)) / 2;
	if (handles.geog)
        [sclat,sclon] = scltln(rlat);
        dx = handles.head_Z1(8) * sclon;
        dy = handles.head_Z1(9) * sclat;
        handles.scaled_dx = dx;     handles.scaled_dy = dy;
        handles.is_meters = 0;      handles.is_km = 0;
	else        % Guess if grid units are meters or km
        dx = handles.head_Z1(2) - handles.head_Z1(1);
        dy = handles.head_Z1(4) - handles.head_Z1(3);
        len = sqrt(dx.*dx + dy.*dy);         % Distance in user unites
        if (len > 1e5)      % If grid's diagonal > 1e5 consider we have meters
            handles.is_meters = 1;     handles.is_km = 0;   handles.geog = 0;
        else                % km
            handles.is_meters = 0;     handles.is_km = 1;   handles.geog = 0;
        end
        handles.scaled_dx = handles.head_Z1(8);
        handles.scaled_dy = handles.head_Z1(9);
        if (handles.is_km)
            handles.scaled_dx = handles.scaled_dx * 1000;
            handles.scaled_dy = handles.scaled_dy * 1000;
        end
	end
    if (~grid_in_continue)       % Called in the "quick mode" 
        if (any(strcmp(mode,{'Power' 'Autocorr' 'Amplitude'})))
            sectrumFun(handles, handles.Z1, handles.head_Z1, mode)
        end
        delete(hObject)
        return
    end
end

if (grid_in_continue)       % Grid recieved in argument. Fill the listboxes
	% The easeast way of not leting the user screw things by selecting a nnx and/or nny
	% lower than nx or ny is to delete the forbiden numbers from the listboxes
	ind = find(cat(1,nlist{:}) > handles.orig_ncols);
	nlist_t = [{handles.orig_ncols}; nlist(ind)];
	ind = find(cat(1,nlist_t{:}) == handles.new_nx);        % Find index of new_nx
	set(handles.listbox_nnx,'String',nlist_t,'Value',ind)
	set(handles.edit_Ncols,'string',num2str(handles.new_nx))
	
	ind = find(cat(1,nlist{:}) > handles.orig_nrows);
	nlist_t = [{handles.orig_nrows}; nlist(ind)];
	ind = find(cat(1,nlist_t{:}) == handles.new_ny);        % Find index of new_ny
	set(handles.listbox_nny,'String',nlist_t,'Value',ind)
	set(handles.edit_Nrows,'string',num2str(handles.new_ny))
    
    handles.X = linspace(handles.head_Z1(1),handles.head_Z1(2),handles.orig_ncols);
    handles.Y = linspace(handles.head_Z1(3),handles.head_Z1(4),handles.orig_nrows);
    set(handles.edit_Grid1,'String','In Memory array')
end

% Import icons
load([handles.path_data 'mirone_icons.mat'],'Mfopen_ico');
set(handles.pushbutton_Grid2,'CData',Mfopen_ico)
set(handles.pushbutton_Grid1,'CData',Mfopen_ico)
clear Mfopen_ico;

% Set upt some useful tooltips
str = sprintf(['The default value is the number of rows in the grid\n',...
    'However, for reducing border effects you may want to apply\n',...
    'a skirt to the grid. For that, select a value from the side\n',...
    'listbox. Extra points will be padded by mirrowiong the west side.']);
set(handles.edit_Nrows,'TooltipString',str)
str = sprintf(['The default value is the number of cols in the grid\n',...
    'However, for reducing border effects you may want to apply\n',...
    'a skirt to the grid. For that, select a value from the side\n',...
    'listbox. Extra points will be padded by mirrowiong the south side.']);
set(handles.edit_Ncols,'TooltipString',str)

str = sprintf('Good FFT numbers for padding the grid');
set(handles.listbox_nnx,'TooltipString',str)
set(handles.listbox_nny,'TooltipString',str)

% Give a Pro look (3D) to the frame boxes 
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

% Choose default command line output for fft_stuff_export
handles.output = hObject;
guidata(hObject, handles);
% UIWAIT makes fft_stuff_export wait for user response (see UIRESUME)
% uiwait(handles.figure1);

set(hObject,'Visible','on');
% NOTE: If you make uiwait active you have also to uncomment the next three lines
% handles = guidata(hObject);
% out = fft_stuff_OutputFcn(hObject, [], handles);
% varargout{1} = out;

% --- Outputs from this function are returned to the command line.
function varargout = fft_stuff_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure

% Get default command line output from handles structure
varargout{1} = handles.output;

% -------------------------------------------------------------------------------------------------
function edit_Grid1_Callback(hObject, eventdata, handles)
fname = get(hObject,'String');
if isempty(fname),   handles.Z1 = [];    return;     end
% Let the pushbutton_Grid1_Callback do all the work
fft_stuff('pushbutton_Grid1_Callback',gcbo,[],guidata(gcbo),fname)

% -------------------------------------------------------------------------------------------------
function pushbutton_Grid1_Callback(hObject, eventdata, handles,opt)
if (nargin == 3),    opt = [];    end
if (nargin == 4),    fname = opt;    end

if (isempty(opt))       % Otherwise 'opt' already transmited the file name.
    if (~isempty(handles.h_calling_fig))                    % If we know the handle to the calling fig
        cfig_handles = guidata(handles.h_calling_fig);      % get handles of the calling fig
        last_dir = cfig_handles.last_dir;
        home = cfig_handles.home_dir;
    else
        last_dir = [];
    end

    if (~isempty(last_dir)),    cd(last_dir);   end
    [FileName,PathName] = uigetfile({'*.grd;*.GRD', 'Grid files (*.grd,*.GRD)';'*.*', 'All Files (*.*)'},'Select GMT grid');
    pause(0.01);
    if (~isempty(last_dir)),    cd(home);   end
    if isequal(FileName,0);     return;     end
    fname = [PathName FileName];
end

% Because GMT and Surfer share the .grd extension, find out which kind grid we are dealing with
[fid, msg] = fopen(fname, 'r');
if fid < 0
    errordlg([PathName FileName ': ' msg],'ERROR'); return
end
ID = fread(fid,4,'*char');
ID = strread(ID,'%s');
if (strcmp(ID,'DSBB') || strcmp(ID,'DSRB'))
    fname = [fname '=6'];
elseif strcmp(ID,'DSAA')
    warndlg('I don''t know and do not intend to learn how to read ASCII Surfer grids.','Warning')
    return
else        % It must (we hope) be a gmt grid
end
[handles.X,handles.Y,handles.Z1,handles.head_Z1] = grdread_m(fname);
if (grdutils(handles.Z1,'-N'))
    errordlg('This grid has NaNs. That is not allowed in FFTs','Error');    return;
end

% See if Grid2 is already loaded and, if yes, if both grids are compatible
if (~isempty(get(handles.edit_Grid2,'String')))
    if ( abs(handles.head_Z1(1) - handles.head_Z1(1)) > 1e-4 || abs(handles.head_Z1(2) - handles.head_Z1(2)) > 1e-4 ||...
            abs(handles.head_Z1(3) - handles.head_Z1(3)) > 1e-4 || abs(handles.head_Z1(4) - handles.head_Z1(4)) > 1e-4)
        errordlg('Error: Grid1 & Grid2 do not cover the same region','Error');  return
    elseif(abs(handles.head_Z1(8) - handles.head_Z1(8)) > 1e-6 || abs(handles.head_Z1(9) - handles.head_Z1(9)) > 1e-6)
        errordlg('Error: Grid1 & Grid2 do not have the same size.','Error');     return
    end
end
[handles.orig_nrows,handles.orig_ncols] = size(handles.Z1);
set(handles.edit_Grid1,'String',fname)

[ns,nlist] = mboard([],handles.orig_ncols,handles.orig_nrows,0,0);
handles.new_nx = ns(1);
handles.new_ny = ns(2);

% The easeast way of not leting the user screw things by selecting a nnx and/or nny lower
% than nx or ny is to delete the forbiden numbers from the listboxes
ind = find(cat(1,nlist{:}) > handles.orig_ncols);
nlist_t = [{handles.orig_ncols}; nlist(ind)];
ind = find(cat(1,nlist_t{:}) == handles.new_nx);        % Find index of new_nx
set(handles.listbox_nnx,'String',nlist_t,'Value',ind)
set(handles.edit_Ncols,'string',num2str(handles.new_nx))

ind = find(cat(1,nlist{:}) > handles.orig_nrows);
nlist_t = [{handles.orig_nrows}; nlist(ind)];
ind = find(cat(1,nlist_t{:}) == handles.new_ny);        % Find index of new_ny
set(handles.listbox_nny,'String',nlist_t,'Value',ind)
set(handles.edit_Nrows,'string',num2str(handles.new_ny))

% Try to guess if grid is in geogs
if (abs(handles.head_Z1(2)-handles.head_Z1(1)) < 180 || abs(handles.head_Z1(4)-handles.head_Z1(3)) < 170)
    handles.geog = 1;   % We probably have a geog grid
    handles.is_meters = 0;     handles.is_km = 0;
    set(handles.popup_GridCoords,'Value',1)
else
    dx = handles.head_Z1(2) - handles.head_Z1(1);
    dy = handles.head_Z1(4) - handles.head_Z1(3);
    len = sqrt(dx.*dx + dy.*dy);         % Distance in user unites
    if (len > 1e5)      % If grid's diagonal > 1e5 consider we have meters
        handles.is_meters = 1;     handles.is_km = 0;   handles.geog = 0;
        set(handles.popup_GridCoords,'Value',2)
    else                % km
        handles.is_meters = 0;     handles.is_km = 1;   handles.geog = 0;
        set(handles.popup_GridCoords,'Value',3)
    end
end

rlat = (handles.head_Z1(4) + handles.head_Z1(3)) / 2;

if (handles.geog)
    [sclat,sclon] = scltln(rlat);
    dx = handles.head_Z1(8) * sclon;
    dy = handles.head_Z1(9) * sclat;
    handles.scaled_dx = dx;     handles.scaled_dy = dy;
else
    handles.scaled_dx = handles.head_Z1(8);
    handles.scaled_dy = handles.head_Z1(9);
    if (handles.is_km)
        handles.scaled_dx = handles.scaled_dx * 1000;
        handles.scaled_dy = handles.scaled_dy * 1000;
    end
end
guidata(hObject,handles)

% -------------------------------------------------------------------------------------------------
function edit_Grid2_Callback(hObject, eventdata, handles)
fname = get(hObject,'String');
if isempty(fname),    handles.Z2 = [];    return;     end
% Let the pushbutton_Grid2_Callback do all the work
fft_stuff('pushbutton_Grid2_Callback',gcbo,[],guidata(gcbo),fname)

% -------------------------------------------------------------------------------------------------
function pushbutton_Grid2_Callback(hObject, eventdata, handles, opt)
if (nargin == 3),    opt = [];       end
if (nargin == 4),    fname = opt;    end

if (isempty(opt))       % Otherwise 'opt' already transmited the file name.
    if (~isempty(handles.h_calling_fig))                    % If we know the handle to the calling fig
        cfig_handles = guidata(handles.h_calling_fig);      % get handles of the calling fig
        last_dir = cfig_handles.last_dir;
        home = cfig_handles.home_dir;
    else
        last_dir = [];
    end

    if (~isempty(last_dir)),    cd(last_dir);   end
    [FileName,PathName] = uigetfile({'*.grd;*.GRD', 'Grid files (*.grd,*.GRD)';'*.*', 'All Files (*.*)'},'Select GMT grid');
    pause(0.01);
    if (~isempty(last_dir)),    cd(home);   end
    if isequal(FileName,0);     return;     end
    fname = [PathName FileName];
end

% Because GMT and Surfer share the .grd extension, find out which kind grid we are dealing with
[fid, msg] = fopen(fname, 'r');
if (fid < 0),    errordlg([fname ': ' msg],'ERROR');     return;     end
ID = fread(fid,4,'*char');
ID = strread(ID,'%s');      fclose(fid);
if (strcmp(ID,'DSBB') || strcmp(ID,'DSRB'))
    fname = [fname '=6'];
elseif strcmp(ID,'DSAA')
    warndlg('I don''t know and do not intend to learn how to read ASCII Surfer grids.','Warning')
    return
else        % It must (we hope) be a gmt grid
end
[X,Y,handles.Z2,handles.head_Z2] = grdread_m(fname);
if (grdutils(handles.Z2,'-N'))
    errordlg('This grid has NaNs. That is not allowed in FFTs','Error');    return;
end
% See if Grid1 grid is already loaded and, if yes, if they are compatible
if (~isempty(get(handles.edit_Grid1,'String')))
    if ( abs(handles.head_Z2(1) - handles.head_Z2(1)) > 1e-4 || abs(handles.head_Z2(2) - handles.head_Z2(2)) > 1e-4 ||...
            abs(handles.head_Z2(3) - handles.head_Z2(3)) > 1e-4 || abs(handles.head_Z2(4) - handles.head_Z2(4)) > 1e-4)
        errordlg('Error: Grid1 & Grid2 do not cover the same region','Error');  return
    elseif(abs(handles.head_Z2(8) - handles.head_Z2(8)) > 1e-6 || abs(handles.head_Z2(9) - handles.head_Z2(9)) > 1e-6)
        errordlg('Error: Grid1 & Grid2 do not have the same size.','Error');     return
    end
end
set(handles.edit_Grid2,'String',fname)
guidata(hObject,handles)

% -------------------------------------------------------------------------------------------------
function checkbox_leaveTrend_Callback(hObject, eventdata, handles)

% -------------------------------------------------------------------------------------------------
function listbox_nnx_Callback(hObject, eventdata, handles)
contents = get(hObject,'String');
nnx = str2double(contents{get(hObject,'Value')});
set(handles.edit_Ncols,'String',num2str(nnx))
handles.new_nx = nnx;     guidata(hObject,handles)

% -------------------------------------------------------------------------------------------------
function listbox_nny_Callback(hObject, eventdata, handles)
contents = get(hObject,'String');
nny = str2double(contents{get(hObject,'Value')});
set(handles.edit_Nrows,'String',num2str(nny))
handles.new_ny = nny;     guidata(hObject,handles)

% -------------------------------------------------------------------------------------------------
function edit_Ncols_Callback(hObject, eventdata, handles)
xx = str2double(get(hObject,'String'));
if (isempty(get(hObject,'String')))
    try,    set(hObject,'String',num2str(handles.ncols));   return;     end
end
if (xx < handles.cols)
    set(hObject,'String',num2str(handles.ncols));    return;
end
handles.ncols = xx;     guidata(hObject,handles)

% -------------------------------------------------------------------------------------------------
function edit_Nrows_Callback(hObject, eventdata, handles)
xx = str2double(get(hObject,'String'));
if (isnan(xx))
    try    set(hObject,'String',num2str(handles.nrows));   return;     end
end
if (xx < handles.nrows)
    set(hObject,'String',num2str(handles.nrows));   return;
end
handles.nrows = xx;     guidata(hObject,handles)

% -------------------------------------------------------------------------------------------------
function edit_UDcont_Callback(hObject, eventdata, handles)
% Do nothing here. If junk in -> junk out

% -------------------------------------------------------------------------------------------------
function edit_dirDerivative_Callback(hObject, eventdata, handles)
azim = str2double(get(hObject,'String'));
if (isnan(azim)),   set(hObject,'String','0');  end

% -------------------------------------------------------------------------------------------------
function edit_derivative_Callback(hObject, eventdata, handles)
n_der = str2double(get(hObject,'String'));
if (isnan(n_der)),   set(hObject,'String','1');  end

% -------------------------------------------------------------------------------------------------
function pushbutton_powerSpectrum_Callback(hObject, eventdata, handles)
if (isempty(handles.Z1))
    errordlg('No grid loaded. Rebooting ...','Error');  return
end
sectrumFun(handles, handles.Z1, handles.head_Z1, 'Power')

% --------------------------------------------------------------------
function pushbutton_crossSpectra_Callback(hObject, eventdata, handles)
if (isempty(handles.Z1) || isempty(handles.Z2))
    errordlg('"Cross" means one grid against the other. You got the idea?','Error');
    return
end
sectrumFun(handles, handles.Z1, handles.head_Z1, 'CrossPower', handles.Z2)

% -------------------------------------------------------------------------------------------------
function pushbutton_autoCorr_Callback(hObject, eventdata, handles)
sectrumFun(handles, handles.Z1, handles.head_Z1, 'Autocorr')

% -------------------------------------------------------------------------------------------------
function pushbutton_crossCorr_Callback(hObject, eventdata, handles)
if (isempty(handles.Z1) || isempty(handles.Z2))
    errordlg('"Cross" means one grid against the other. You got the idea?','Error');
    return
end
sectrumFun(handles, handles.Z1, handles.head_Z1, 'CrossCorrel', handles.Z2)

% -------------------------------------------------------------------------------------------------
function pushbutton_radialPowerAverage_Callback(hObject, eventdata, handles)
% This is modeled on the 1-D case, using the following ideas:
% In 1-D, we ensemble average over samples of length L = n * dt.
% This gives n/2 power spectral estimates, at frequencies i/L,
% where i = 1, n/2.  If we have a total data set of ndata,
% we can make nest=ndata/n independent estimates of the power.
% Our standard error is then 1/sqrt(nest). In making 1-D estimates
% from 2-D data, we may choose n and L from nx2 or ny2 and delta_kx, delta_ky 
% as appropriate.  In this routine, we are giving the sum over all frequencies
% in the other dimension; that is, an approximation of the integral.
%   ORIGINAL TEXT FROM WALTER SMITH IN GRDFFT

if (isempty(handles.Z1)),    errordlg('No grid loaded yet.','Error');    return; end

nx2 = handles.new_nx;   ny2 = handles.new_ny;
delta_kx = 2*pi / (nx2 * handles.scaled_dx);
delta_ky = 2*pi / (ny2 * handles.scaled_dy);
if (delta_kx < delta_ky)
    delta_k = delta_kx;    nk = fix(nx2/2);
else
    delta_k = delta_ky;    nk = fix(ny2/2);
end
r_delta_k = 1 / delta_k;
set(handles.figure1,'pointer','watch')
if (get(handles.checkbox_leaveTrend,'Value'))       % Remove trend
    handles.Z1 = double(grdtrend_m(single(handles.Z1),handles.head_Z1,'-D','-N3'));
end
[Z,band,modk] = wavenumber_and_mboard(handles);
ifreq = round(modk*r_delta_k) + 1;     clear modk;
Z = fft2(Z);    Z(1,1) = 0;
Z = Z.* conj(Z);

ifreq(ifreq < 1) = 1;           % Might happen when doing r spectrum
nk1 = nk + 1;
power = zeros(nk1,1);
n_used = 0;
for (m=1:ny2)
    for(n=1:nx2)
	    if (ifreq(m,n) > nk1),  continue;	end % Might happen when doing r spectrum
        power(ifreq(m,n)) = power(ifreq(m,n)) + Z(m,n);
        n_used = n_used + 1;
    end
end
power(1) = [];              % Remove DC component

eps_pow = 1 / sqrt(n_used/nk);
delta_k = delta_k / (2*pi);         % Write out frequency, not wavenumber
powfactor = 1 / (nx2*ny2)^2;
power = power * powfactor;
eps_pow = eps_pow * power;          % Wrong for admitance
freq = (1:nk) * delta_k;
%freq = 1 ./ freq;               % Wavelength
if (handles.geog),   freq = freq * 1000;     end     % Report frequency in 1/km

h=figure('NumberTitle','off','Name','Radial average power spectrum',...
    'Color',get(0,'factoryUicontrolBackgroundColor'));
ud.h_mir_fig = handles.h_calling_fig;
ud.data = [freq; power'; eps_pow'];     % ready for the stupid fwrite
set(h,'UserData',ud)
uimenu('Label','Save data','callback',{@calSave,h});
semilogy(freq,power)
if (handles.geog || handles.is_km)
    xlabel('Frequency (1/km)','FontSize',12)
else
    xlabel('Frequency (1/m)','FontSize',12)
end
ylabel('Log(Power)','FontSize',12)
set(handles.figure1,'pointer','arrow')

% -------------------------------------------------------------------------------
function calSave(obj,eventdata,h_fig)
% Save data in file
ud = get(h_fig,'UserData');          % Get userdata
if (~isempty(ud.h_mir_fig))
	handles_mir = guidata(ud.h_mir_fig);    % Get the Mirone fig handles
	work_dir = handles_mir.work_dir;
	home_dir = handles_mir.home_dir;
	cd(work_dir)
end
[FileName,PathName] = uiputfile({ ...
    '*.dat;*.DAT', 'radial spectrum file (*.dat,*.DAT)'; '*.*', 'All Files (*.*)'}, 'Select File name');
try     % If fft_stuff was called directly home_dir does not exist
    cd(home_dir);       % allways go home
end
if isequal(FileName,0);   return;     end
pause(0.01)
fname = [PathName FileName];
[PATH,FNAME,EXT] = fileparts(fname);
if (isempty(EXT)),    fname = [fname '.dat'];    end

fid = fopen(fname, 'w');
if (fid < 0),    errordlg(['Can''t open file:  ' fname],'Error');    return;     end
fprintf(fid,'%f\t%f\t%f\n',[ud.data(1,:); ud.data(2,:); ud.data(3,:)]);
fclose(fid);

% -------------------------------------------------------------------------------------------------
function pushbutton_integrate_Callback(hObject, eventdata, handles, opt)
if (isempty(handles.Z1)),    errordlg('No grid loaded yet.','Error');    return; end
if (isempty(opt))
    scale = 1;
else
    scale = 980619.9203;    % Moritz's 1980 IGF value for gravity in mGal at 45 degrees latit
end
set(handles.figure1,'pointer','watch')
if (get(handles.checkbox_leaveTrend,'Value'))       % Remove trend
    handles.Z1 = double(grdtrend_m(single(handles.Z1),handles.head_Z1,'-D','-N3'));
end
[Z,band,k] = wavenumber_and_mboard(handles);
k(1,1) = eps;       % Avoid a devided by zero warning
Z = fft2(Z) ./ (k*scale);    Z(1,1) = 0;    clear k;
Z = real(ifft2(Z));
[Z,tmp] = unband(handles,Z,band);
if (scale == 1), tmp.name = 'Integrated grid';
else                tmp.name = 'Geoid height';
end
set(handles.figure1,'pointer','arrow')
mirone(single(Z),tmp);

% --------------------------------------------------------------------
function sectrumFun(handles, Z, head, opt1, Z2)
% OPT1 = 'Amplitude'   -> compute amplitude spectrum
% OPT1 = 'Power'       -> compute power spectrum
% OPT1 = 'Autocorr'    -> compute autocorrelation
% OPT1 = 'CrossPower'  -> compute cross power spectra between Z & Z2
% OPT1 = 'CrossCorrel' -> compute cross correlation between Z & Z2
% Z2 -> needed when OPT1 = 'CrossPower' OR OPT1 = 'CrossCorrel'

two_grids = 0;
if (nargin == 5),    two_grids = 1;  end
set(handles.figure1,'pointer','watch')
if (get(handles.checkbox_leaveTrend,'Value'))       % Remove trend
    Z = double(grdtrend_m(single(Z),handles.head_Z1,'-D','-N3'));
    if (two_grids)
        Z2 = double(grdtrend_m(single(Z2),handles.head_Z1,'-D','-N3'));
    end
end
nx = handles.orig_ncols;        ny = handles.orig_nrows;
[Z,band] = mboard(double(Z),nx,ny,handles.new_nx,handles.new_ny);
if (two_grids)
    Z2 = mboard(double(Z2),nx,ny,handles.new_nx,handles.new_ny);
end
if (any(strcmp(opt1,{'Amplitude' 'Power' 'CrossPower'})))   % For these this is more efficient
    Z = fftshift(fft2(Z));
    m1 = band(1)+1;     m2 = m1 + ny - 1;
    n1 = band(3)+1;     n2 = n1 + nx - 1;
    %Z = ifft2(ifftshift(Z));
    Z = Z(m1:m2,n1:n2);                 % Remove the padding skirt
    if (two_grids)
        Z2 = fftshift(fft2(Z2));        Z2 = Z2(m1:m2,n1:n2);
    end
end

if (strcmp(opt1,'Power'))
	Z = log10( (Z .* conj(Z) / (nx*ny)) + 1);       % An offset of 1 is added to avoid log(0)
    tmp.name = 'Power spectrum';
elseif (strcmp(opt1,'CrossPower'))                  % Cross spectra
	Z = log10( (real(Z).*real(Z2) + imag(Z).*imag(Z2)) / (nx*ny) + 1);  % An offset of 1 is added to avoid log(0)
    tmp.name = 'Cross spectra';    clear Z2;
elseif (strcmp(opt1,'Amplitude'))                   % Amplitude spectrum
    Z = log10(abs(Z) / (nx*ny) + 1);
    tmp.name = 'Amplitude spectrum';
elseif (strcmp(opt1,'Autocorr'))                    % Autocorrelation
    Z = fftshift(real(ifft2(abs(fft2(Z)).^2)));
    m1 = band(1)+1;     m2 = m1 + ny - 1;
    n1 = band(3)+1;     n2 = n1 + nx - 1;
    Z = Z(m1:m2,n1:n2);                             % Remove the padding skirt
    tmp.name = 'Autocorrelation';
elseif (strcmp(opt1,'CrossCorrel'))                 % Cross Correlation
    Z = ifft2(abs(fft2(Z) .* fft2(Z2)));
    Z = fftshift(real(Z));
    m1 = band(1)+1;     m2 = m1 + ny - 1;
    n1 = band(3)+1;     n2 = n1 + nx - 1;
    Z = Z(m1:m2,n1:n2);                             % Remove the padding skirt
    tmp.name = 'Cross Correlation';
end

if ( (strcmp(opt1,'Autocorr') || strcmp(opt1,'CrossCorrel')) )
    delta_kx = 1;   delta_ky = 1;
else            % make wave number array
    delta_kx = 2*pi / (handles.new_nx * handles.scaled_dx);
    delta_ky = 2*pi / (handles.new_ny * handles.scaled_dy);
end
if (handles.is_km)      % In km case scaled_dx|dy where in meters
    delta_kx = delta_kx * 1000;
    delta_ky = delta_ky * 1000;
end

nx2 = fix(nx/2);    ny2 = fix(ny/2);
if (rem(nx,2) == 0),   sft_x = 1;
else                   sft_x = 0;     end
if (rem(ny,2) == 0),   sft_y = 1;
else                   sft_y = 0;     end
tmp.X = (-nx2:nx2-sft_x).*delta_kx;     tmp.Y = (-ny2:ny2-sft_y).*delta_ky;
tmp.head = [tmp.X(1) tmp.X(end) tmp.Y(1) tmp.Y(end) min(Z(:)) max(Z(:)) 0 delta_kx delta_ky];
set(handles.figure1,'pointer','arrow')
mirone(single(Z),tmp);

% --------------------------------------------------------------------
function pushbutton_goUDcont_Callback(hObject, eventdata, handles)
if (isempty(handles.Z1)),    errordlg('No grid loaded yet.','Error');    return; end
zup = str2double(get(handles.edit_UDcont,'String'));
if (isnan(zup)),     return;     end
set(handles.figure1,'pointer','watch')
if (get(handles.checkbox_leaveTrend,'Value'))       % Remove trend
    handles.Z1 = double(grdtrend_m(single(handles.Z1),handles.head_Z1,'-D','-N3'));
end
[Z,band,k] = wavenumber_and_mboard(handles);
Z = fft2(Z) .* exp(-k.*zup);    clear k;
Z = real(ifft2((Z)));
[Z,tmp] = unband(handles,Z,band);
tmp.name = 'U/D Continuation';
set(handles.figure1,'pointer','arrow')
mirone(single(Z),tmp);

% --------------------------------------------------------------------
function pushbutton_goDerivative_Callback(hObject, eventdata, handles, opt)
% Compute N-vertical derivative
if (isempty(handles.Z1)),   errordlg('No grid loaded yet.','Error');    return; end
if (isempty(opt)),   scale = 1;
else                scale = 980619.9203;
end                 % Moritz's 1980 IGF value for gravity in mGal at lat = 45 degrees
n_der = fix(str2double(get(handles.edit_derivative,'String')));
set(handles.figure1,'pointer','watch')
if (get(handles.checkbox_leaveTrend,'Value'))       % Remove trend
    handles.Z1 = double(grdtrend_m(single(handles.Z1),handles.head_Z1,'-D','-N3'));
end
[Z,band,k] = wavenumber_and_mboard(handles);
if (scale == 980619.9203),   n_der = 1;  end
if (n_der > 1),  k = k.^n_der;   end
Z = fft2(Z) .* k*scale;    Z(1,1) = 0;    clear k;
Z = real(ifft2((Z)));
[Z,tmp] = unband(handles,Z,band);
if (scale == 1),  tmp.name = [num2str(n_der) 'th Vertical Derivative'];
else                 tmp.name = 'Gravity anomaly (mGal)';
end
set(handles.figure1,'pointer','arrow')
mirone(single(Z),tmp);

% --------------------------------------------------------------------
function pushbutton_goDirDerivative_Callback(hObject, eventdata, handles)
% Compute a directional derivative
if (isempty(handles.Z1)),    errordlg('No grid loaded yet.','Error');    return; end
azim = str2double(get(handles.edit_dirDerivative,'String'));
set(handles.figure1,'pointer','watch')
if (get(handles.checkbox_leaveTrend,'Value'))       % Remove trend
    handles.Z1 = double(grdtrend_m(single(handles.Z1),handles.head_Z1,'-D','-N3'));
end
[Z,band,k] = wavenumber_and_mboard(handles,1);
fact = (sin(azim*pi/180) * k.x + cos(azim*pi/180) * k.y);   clear k;
Z = fft2(Z);    Z(1,1) = 0;
Z = complex(-imag(Z) .* fact, real(Z) .* fact);
Z = real(ifft2((Z)));
[Z,tmp] = unband(handles,Z,band);
tmp.name = [num2str(azim) ' Azimuthal Derivative'];
set(handles.figure1,'pointer','arrow')
mirone(single(Z),tmp);

% --------------------------------------------------------------------
function [Z,band,k] = wavenumber_and_mboard(handles,opt)
% Compute the wavenumber array and call mboard
% Z is the expanded array (to handles.new_nx & handles.new_ny)
% BAND is the second output od mboard
% K is the wavenumber array
% If opt = 1 K will be instead a structure with KX & KY fields

if (nargin == 1),    opt = [];   end
% calculate wavenumber array
nx2 = fix(handles.new_nx/2);     ny2 = fix(handles.new_ny/2);

if (rem(handles.new_nx,2) == 0), sft_x = 1;
else                            sft_x = 0;     end
if (rem(handles.new_ny,2) == 0), sft_y = 1;
else                            sft_y = 0;     end

%nx2 = handles.new_nx/2;     ny2 = handles.new_ny/2;    % Tivey way -> but and if n is odd?
dkx = 2*pi / (handles.new_nx * handles.scaled_dx);
dky = 2*pi / (handles.new_ny * handles.scaled_dy);
%kx = (-nx2:nx2-1).*dkx;     ky = (-ny2:ny2-1).*dky;    % Tivey way
kx = (-nx2:nx2-sft_x).*dkx;  ky = (-ny2:ny2-sft_y).*dky;
%X = ones(size(ky))'*kx;     Y = ky'*ones(size(kx));
X = repmat(kx,length(ky),1);    Y = repmat(ky',1,length(kx));
if (opt)
    k.x = ifftshift(X);     clear X;
    k.y = ifftshift(Y);     clear Y;
else
    k = ifftshift(sqrt(X.^2+Y.^2));      % wavenumber array   (Tivey used fftshift)
    clear X Y;
end
[Z,band] = mboard(double(handles.Z1),handles.orig_ncols,handles.orig_nrows,handles.new_nx,handles.new_ny);

% --------------------------------------------------------------------
function [Z,hdr] = unband(handles,Z,band)
% Take the BAND vector and remove the padding skirt previously applyed
% (by mboard) to the Z matrix
% Return also the HDR structure used the call a new Mirone window
m1 = band(1)+1;     m2 = m1 + handles.orig_nrows - 1;
n1 = band(3)+1;     n2 = n1 + handles.orig_ncols - 1;
Z = Z(m1:m2,n1:n2);                 % Remove the padding skirt
hdr.X = handles.X;      hdr.Y = handles.Y;
hdr.head = handles.head_Z1;
hdr.head(5) = min(min(Z));      hdr.head(6) = max(max(Z));

% --------------------------------------------------------------------
function popup_GridCoords_Callback(hObject, eventdata, handles)
xx = get(hObject,'Value');
if (xx == 1),        handles.geog = 1;       handles.is_meters = 0;  handles.is_km = 0;
elseif (xx == 2),   handles.is_meters = 1;  handles.is_geog = 0;    handles.is_km = 0;
elseif (xx == 3)
    handles.is_km = 1;      handles.is_geog = 0;    handles.is_meters = 0;
    handles.scaled_dx = handles.scaled_dx * 1000;
    handles.scaled_dy = handles.scaled_dy * 1000;
end
guidata(hObject,handles)

% -------------------------------------------------------------------------------------------------
function [sclat,sclon] = scltln(orlat)
% Routine to determine lat-lon scales, km/deg, for ellipsoids
% of revolution,  using equations of:
%       Snyder, J.P., 1987, Map Projections -- A Working Manual,
%       USGS Professional Paper 1395, Washington DC, 383p. cf. pp 24-25.
%
% Currently, this is hard-wired for the WGS-84 ellipsoid.
%
% The basic equations are:
% 	sclat = a * (1-e*e)    /  (1 - e*e * sin(orlat)*sin(orlat))**(1.5)
%	sclon = a * cos(orlat) /  (1 - e*e * sin(orlat)*sin(orlat))**(0.5)
%
% where:    a  is the equatorial radius
%           b  is the polar radius
%           e  is the eccentricity
%           f  is the flattening
% Also:
%	e*e = 1. - b*b/a*a
%	f   = 1. - b/a
%
% Dan Scheirer, 21 May 1991

% These constants belong to the: WGS, 1984 ellipsoid (gmt_defaults.h)
%a = 6378.137;   b = 6356.7521;
a = 6378137;    b = 6356752.1;      % Use meters

% Now, do the calculations...
e2 = 1 - (b*b)/(a*a);
sinlat = sin(orlat*pi/180);
denom  = sqrt(1 - e2 * sinlat * sinlat);
sclat = (pi/180) * a * (1 - e2)  / denom / denom / denom;
sclon = (pi/180) * a * cos(orlat*pi/180) / denom;


% --- Creates and returns a handle to the GUI figure. 
function fft_stuff_LayoutFcn(h1,handles)

set(h1,'PaperUnits',get(0,'defaultfigurePaperUnits'),...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'MenuBar','none',...
'Name','FFT stuff',...
'NumberTitle','off',...
'Position',[520 579.000000000001 570 221],...
'Renderer',get(0,'defaultfigureRenderer'),...
'RendererMode','manual',...
'Resize','off',...
'Tag','figure1',...
'UserData',[]);

uicontrol('Parent',h1,'Position',[330 139 231 79],'String',{''},'Style','frame','Tag','frame2');

uicontrol('Parent',h1,'Position',[10 139 301 79],'String',{''},'Style','frame','Tag','frame1');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@fft_stuff_uicallback,h1,'edit_Grid1_Callback'},...
'HorizontalAlignment','left',...
'Position',[50 187 231 22],...
'Style','edit',...
'TooltipString','Enter grid 1',...
'Tag','edit_Grid1');

uicontrol('Parent',h1,...
'Callback',{@fft_stuff_uicallback4,h1,[],'pushbutton_Grid1_Callback'},...
'Position',[278 185 23 23],...
'Tag','pushbutton_Grid1');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@fft_stuff_uicallback,h1,'edit_Grid2_Callback'},...
'HorizontalAlignment','left',...
'Position',[50 157 231 22],...
'Style','edit',...
'TooltipString','Enter grid 2 (ONLY USED IN "CROSS" OPERATIONS)',...
'Tag','edit_Grid2');

uicontrol('Parent',h1,...
'Callback',{@fft_stuff_uicallback4,h1,[],'pushbutton_Grid2_Callback'},...
'Position',[278 155 23 23],...
'Tag','pushbutton_Grid2');

uicontrol('Parent',h1,'HorizontalAlignment','left',...
'Position',[17 161 31 15],'String','Grid2','Style','text','Tag','text1');

uicontrol('Parent',h1,'HorizontalAlignment','left',...
'Position',[18 190 26 16],'String','Grid1','Style','text','Tag','text_FieldMag');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@fft_stuff_uicallback,h1,'listbox_nny_Callback'},...
'Position',[390 147 51 60],'Style','listbox','Value',1,'Tag','listbox_nny');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@fft_stuff_uicallback,h1,'listbox_nnx_Callback'},...
'Position',[498 150 51 60],'Style','listbox','Value',1,'Tag','listbox_nnx');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@fft_stuff_uicallback,h1,'edit_Nrows_Callback'},...
'Position',[350 168 41 21],...
'Style','edit',...
'TooltipString','Number of grid rows',...
'Tag','edit_Nrows');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@fft_stuff_uicallback,h1,'edit_Ncols_Callback'},...
'Position',[456 168 41 21],...
'Style','edit',...
'TooltipString','Number of grid columns',...
'Tag','edit_Ncols');

uicontrol('Parent',h1,'Position',[351 192 39 15],'String','# Rows',...
'Style','text','Tag','text11');

uicontrol('Parent',h1,'Position',[456 192 39 15],'String','# Cols',...
'Style','text','Tag','text12');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@fft_stuff_uicallback,h1,'edit_UDcont_Callback'},...
'Position',[127 73 41 21],...
'Style','edit',...
'TooltipString','Height of continuation (in meters)',...
'Tag','edit_UDcont');

uicontrol('Parent',h1,'HorizontalAlignment','right','Position',[11 75 112 17],...
'String','Up/Down Continuation','Style','text','Tag','text13');

uicontrol('Parent',h1,...
'Callback',{@fft_stuff_uicallback,h1,'pushbutton_powerSpectrum_Callback'},...
'Position',[240 109 101 23],...
'String','Power Spectrum',...
'Tag','pushbutton_powerSpectrum');

uicontrol('Parent',h1,...
'Callback',{@fft_stuff_uicallback,h1,'pushbutton_autoCorr_Callback'},...
'Position',[240 76 101 23],...
'String','Auto Correlation',...
'Tag','pushbutton_autoCorr');

uicontrol('Parent',h1,...
'Callback',{@fft_stuff_uicallback,h1,'pushbutton_crossCorr_Callback'},...
'Position',[360 76 121 23],...
'String','Cross Correlation',...
'Tag','pushbutton_crossCorr');

uicontrol('Parent',h1,...
'Callback',{@fft_stuff_uicallback,h1,'pushbutton_radialPowerAverage_Callback'},...
'Position',[360 43 121 23],...
'String','Radial power average',...
'Tag','pushbutton_radialPowerAverage');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@fft_stuff_uicallback,h1,'edit_derivative_Callback'},...
'Position',[127 43 41 21],...
'String','1',...
'Style','edit',...
'TooltipString','Order of vertical differentiation',...
'Tag','edit_derivative',...
'UserData',[]);

uicontrol('Parent',h1,'HorizontalAlignment','right','Position',[21 43 96 19],...
'String','Derivative (N-order)','Style','text','Tag','text14','UserData',[]);

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@fft_stuff_uicallback,h1,'edit_dirDerivative_Callback'},...
'Position',[127 13 41 21],...
'String','0',...
'Style','edit',...
'TooltipString','Angle of directional derivative',...
'Tag','edit_dirDerivative');

uicontrol('Parent',h1,'HorizontalAlignment','right','Position',[20 13 100 19],...
'String','Directional derivative','Style','text','Tag','text16');

uicontrol('Parent',h1,...
'Callback',{@fft_stuff_uicallback4,h1,[],'pushbutton_integrate_Callback'},...
'Position',[240 43 101 23],...
'String','Integrate',...
'Tag','pushbutton_integrate');

uicontrol('Parent',h1,...
'Callback',{@fft_stuff_uicallback4,h1,'FAA2Geoid','pushbutton_integrate_Callback'},...
'Position',[240 10 101 23],...
'String','FAA2Geoid',...
'Tag','pushbutton_faa2geoid');

uicontrol('Parent',h1,...
'Callback',{@fft_stuff_uicallback4,h1,'Geoid2FAA','pushbutton_goDerivative_Callback'},...
'Position',[360 10 121 23],...
'String','Geoid2FAA',...
'Tag','pushbutton_geoid2faa');

uicontrol('Parent',h1,...
'Callback',{@fft_stuff_uicallback,h1,'pushbutton_goUDcont_Callback'},...
'Position',[168 73 31 21],...
'String','Go',...
'Tag','pushbutton_goUDcont');

uicontrol('Parent',h1,...
'Callback',{@fft_stuff_uicallback4,h1,[],'pushbutton_goDerivative_Callback'},...
'Position',[168 43 31 21],...
'String','Go',...
'Tag','pushbutton_goDerivative');

uicontrol('Parent',h1,...
'Callback',{@fft_stuff_uicallback,h1,'pushbutton_goDirDerivative_Callback'},...
'Position',[168 13 31 21],...
'String','Go',...
'TooltipString','Azimuth of directional derivative CCW from north',...
'Tag','pushbutton_goDirDerivative');

uicontrol('Parent',h1,...
'Callback',{@fft_stuff_uicallback,h1,'pushbutton_crossSpectra_Callback'},...
'Position',[360 109 121 23],...
'String','Cross Spectra',...
'Tag','pushbutton_crossSpectra');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@fft_stuff_uicallback,h1,'popup_GridCoords_Callback'},...
'Position',[86 106 82 22],...
'String',{'Geogs'; 'Meters'; 'Kilometers'},...
'Style','popupmenu',...
'TooltipString','GRID COORDINATES: IT IS YOUR RESPONSABILITY THAT THIS IS CORRECT',...
'Value',1,...
'Tag','popup_GridCoords');

uicontrol('Parent',h1,'ForegroundColor',[1 0 0],'Position',[31 109 51 15],...
'String','CONFIRM','Style','text','Tag','text17');

uicontrol('Parent',h1,...
'Callback',{@fft_stuff_uicallback,h1,'checkbox_leaveTrend_Callback'},...
'Position',[50 140 105 15],...
'String','Leave trend alone',...
'Style','checkbox',...
'TooltipString','If checked remove a plane before transformations',...
'Tag','checkbox_leaveTrend');

function fft_stuff_uicallback(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
feval(callback_name,hObject,[],guidata(h1));

function fft_stuff_uicallback4(hObject, eventdata, h1, opt, callback_name)
% This function is executed by the callback and than the handles is allways updated.
feval(callback_name,hObject,[],guidata(h1),opt);

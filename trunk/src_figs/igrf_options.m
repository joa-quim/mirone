function varargout = igrf_options(varargin)
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
	igrf_options_LayoutFcn(hObject,handles);
	handles = guihandles(hObject);
	movegui(hObject,'east')

    if (numel(varargin) > 0)
        handMir = varargin{1};
        d_path = handMir.path_data;
	    handles.h_calling_fig = handMir.figure1;
	else
        d_path = [pwd filesep 'data' filesep];
	    handles.h_calling_fig = [];
    end

	% Import icons
	load([d_path 'mirone_icons.mat'],'Mfopen_ico');
	set(handles.pushbutton_InputFile,'CData',Mfopen_ico)
	set(handles.pushbutton_OutputFile,'CData',Mfopen_ico)
	clear Mfopen_ico;

	handles.input_file = [];        % May contain input total field measurements
	handles.output_file = [];       %
	handles.dms_xinc = 0;
	handles.dms_yinc = 0;
	handles.x_min = -180;
	handles.x_max = 180;
	handles.y_min = -89;
	handles.y_max = 89;

	time = clock;
	set(handles.edit_DateYY,'String',sprintf('%d',time(1)))
	set(handles.edit_DateMM,'String',sprintf('%d',time(2)))
	set(handles.edit_DateDD,'String',sprintf('%d',time(3)))

	yd = dec_year;
	set(handles.edit_DateDec,'String',sprintf('%.3f',yd))

	% Compute the igrf field using the default values and the clock time
	lat = str2double(get(handles.edit_LatDec,'String'));
	lon = str2double(get(handles.edit_LonDec,'String'));

	start = 1900; stop = 2010;

	handles.start_stop_epoch = [start stop];
	out = igrf_m(lon, lat, 0, yd);

	% Fill in the region boxes with default values (a global grid)
	set(handles.edit_Ymin,'String','-89')
	set(handles.edit_Ymax,'String','89')
	set(handles.edit_Xinc,'String','1')
	set(handles.edit_Yinc,'String','1')
	set(handles.edit_Ncols,'String','361')
	set(handles.edit_Nrows,'String','179')

	% Give a Pro look (3D) to the frame boxes 
	bgcolor = get(0,'DefaultUicontrolBackgroundColor');
	framecolor = max(min(0.65*bgcolor,[1 1 1]),[0 0 0]);
	%set(0,'Units','pixels');    % Pixels are easier to reason with
    h_f = [handles.frame1 handles.frame2 handles.frame3 handles.frame4];
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

	% Recopy the text fields on top of previously created frames (uistack is too slow)
    h_t = [handles.text31 handles.text39 handles.text40 handles.text41];
	for i=1:numel(h_t)
        usr_d = get(h_t(i),'UserData');
        t_size = get(h_t(i),'Position');   t_str = get(h_t(i),'String');    fw = get(h_t(i),'FontWeight');
        fs = get(h_t(i),'FontSize');
        bgc = get (h_t(i),'BackgroundColor');   fgc = get (h_t(i),'ForegroundColor');tag=get(h_t(i),'Tag');
        uicontrol('Parent',hObject, 'Style','text', 'Position',t_size,'String',t_str, 'FontSize', fs, ...
            'BackgroundColor',bgc,'ForegroundColor',fgc,'FontWeight',fw,'UserData',usr_d,'Tag',tag);
	end
	delete(h_t)

	% Import the world map
	handles.w_map = flipdim(imread([d_path 'etopo2.jpg']),1);
	image([-180 180],[-90 90],handles.w_map);   set(handles.axes1,'YDir','normal');

	% Fill the various mag field texts
	set(handles.text_FnT,'String',[sprintf('%.0f',out(1)) '  nT'],'FontSize',10)
	set(handles.text_HnT,'String',[sprintf('%.0f',out(2)) '  nT'],'FontSize',10)
	set(handles.text_XnT,'String',[sprintf('%.0f',out(3)) '  nT'],'FontSize',10)
	set(handles.text_YnT,'String',[sprintf('%.0f',out(4)) '  nT'],'FontSize',10)
	set(handles.text_ZnT,'String',[sprintf('%.0f',out(5)) '  nT'],'FontSize',10)
	set(handles.text_Ddeg,'String',[sprintf('%.1f',out(6)) ' º'],'FontSize',10)
	set(handles.text_Ideg,'String',[sprintf('%.1f',out(7)) ' º'],'FontSize',10)

	set(hObject,'WindowButtonDownFcn',{@wbd,handles});

	% Choose default command line output for igrf_options
	handles.output = hObject;
	guidata(hObject, handles);

    set(hObject,'Visible','on');
    if (nargout),   varargout{1} = handles.output;     end

% -------------------------------------------------------------------------------------------------
function wbd(obj,eventdata,handles)
	Pt = get(handles.axes1,'CurrentPoint');
	xlim = get(handles.axes1,'xlim');
	ylim = get(handles.axes1,'ylim');
	if (xlim(1) <= Pt(1,1) && Pt(1,1) <= xlim(2) && ylim(1) <= Pt(1,2) && Pt(1,2) <= ylim(2))
        set(handles.edit_LonDec,'String',sprintf('%.3f',Pt(1,1)))
        set(handles.edit_LatDec,'String',sprintf('%.3f',Pt(1,2)))
        dms = dec2deg(Pt(1,1),'opt');       % Convert longitude decimal degrees in deg,min,sec
        set(handles.edit_LonDeg,'String',sprintf('%.0f',dms(1)))
        set(handles.edit_LonMin,'String',sprintf('%.0f',dms(2)))
        set(handles.edit_LonSec,'String',sprintf('%.2f',dms(3)))
        dms = dec2deg(Pt(1,2),'opt');       % Convert latitude decimal degrees in deg,min,sec
        set(handles.edit_LatDeg,'String',sprintf('%.0f',dms(1)))
        set(handles.edit_LatMin,'String',sprintf('%.0f',dms(2)))
        set(handles.edit_LatSec,'String',sprintf('%.2f',dms(3)))
        elev = str2double(get(handles.edit_Elev,'String'));
        date = str2double(get(handles.edit_DateDec,'String'));
        out = igrf_m(Pt(1,1), Pt(1,2), elev, date);
        set_field_boxes(out, handles)   % Fill the various mag field texts
	end

% -------------------------------------------------------------------------------------------------
function edit_LatDeg_Callback(hObject, eventdata, handles)
	if (isnan(get(hObject,'String'))),    set(hObject,'String','37'); return; end
	lat_deg = str2double(get(hObject,'String'));
	lat_deg = fix(lat_deg);     % Ensure that the value is an integer
	lat_min = str2double(get(handles.edit_LatMin,'String'));
	lat_sec = str2double(get(handles.edit_LatSec,'String'));
	if (lat_deg < 0),   lat = lat_deg - lat_min/60 - lat_sec/3600;
	else                lat = lat_deg + lat_min/60 + lat_sec/3600;      end
	set(handles.edit_LatDec,'String',sprintf('%.4f',lat))
	lon = str2double(get(handles.edit_LonDec,'String'));
	elev = str2double(get(handles.edit_Elev,'String'));
	date = str2double(get(handles.edit_DateDec,'String'));
	out = igrf_m(lon, lat, elev, date);
	set_field_boxes(out, handles)   % Fill the various mag field boxes

% -------------------------------------------------------------------------------------------------
function edit_LatMin_Callback(hObject, eventdata, handles)
	if (isnan(get(hObject,'String'))),    set(hObject,'String','0'); return; end
	lat_min = str2double(get(hObject,'String'));
	if (lat_min < 0 || lat_min >= 60), set(hObject,'String','0'); return; end
	lat_min = fix(lat_min);     % Ensure that the value is an integer
	lat_deg = str2double(get(handles.edit_LatDeg,'String'));
	lat_sec = str2double(get(handles.edit_LatSec,'String'));
	if (lat_deg < 0),   lat = lat_deg - lat_min/60 - lat_sec/3600;
	else                lat = lat_deg + lat_min/60 + lat_sec/3600;      end
	set(handles.edit_LatDec,'String',sprintf('%.4f',lat))
	lon = str2double(get(handles.edit_LonDec,'String'));
	elev = str2double(get(handles.edit_Elev,'String'));
	date = str2double(get(handles.edit_DateDec,'String'));
	out = igrf_m(lon, lat, elev, date);
	set_field_boxes(out, handles)   % Fill the various mag field boxes

% -------------------------------------------------------------------------------------------------
function edit_LatSec_Callback(hObject, eventdata, handles)
	if (isnan(get(hObject,'String'))),    set(hObject,'String','0'); return; end
	lat_sec = str2double(get(hObject,'String'));
	if (lat_sec < 0 || lat_sec >= 60), set(hObject,'String','0'); return; end
	lat_deg = str2double(get(handles.edit_LatDeg,'String'));
	lat_min = str2double(get(handles.edit_LatMin,'String'));
	if (lat_deg < 0),   lat = lat_deg - lat_min/60 - lat_sec/3600;
	else                lat = lat_deg + lat_min/60 + lat_sec/3600;      end
	set(handles.edit_LatDec,'String',sprintf('%.4f',lat))
	lon = str2double(get(handles.edit_LonDec,'String'));
	elev = str2double(get(handles.edit_Elev,'String'));
	date = str2double(get(handles.edit_DateDec,'String'));
	out = igrf_m(lon, lat, elev, date);
	set_field_boxes(out, handles)   % Fill the various mag field boxes

% -------------------------------------------------------------------------------------------------
function edit_LatDec_Callback(hObject, eventdata, handles)
	if (isnan(get(hObject,'String'))),    set(hObject,'String','37'); return; end
	lat = str2double(get(hObject,'String'));
	lon = str2double(get(handles.edit_LonDec,'String'));
	elev = str2double(get(handles.edit_Elev,'String'));
	date = str2double(get(handles.edit_DateDec,'String'));
	dms = dec2deg(lat,'opt');       % Convert decimal degrees in deg,min,sec
	set(handles.edit_LatDeg,'String',sprintf('%.0f',dms(1)))
	set(handles.edit_LatMin,'String',sprintf('%.0f',dms(2)))
	set(handles.edit_LatSec,'String',sprintf('%.2f',dms(3)))
	out = igrf_m(lon, lat, elev, date);
	set_field_boxes(out, handles)   % Fill the various mag field boxes

% -------------------------------------------------------------------------------------------------
function edit_LonDeg_Callback(hObject, eventdata, handles)
	if (isnan(get(hObject,'String'))),    set(hObject,'String','-8'); return; end
	lon_deg = str2double(get(hObject,'String'));
	lon_deg = fix(lon_deg);     % Ensure that the value is an integer
	lon_min = str2double(get(handles.edit_LonMin,'String'));
	lon_sec = str2double(get(handles.edit_LonSec,'String'));
	if (lon_deg < 0),   lon = lon_deg - lon_min/60 - lon_sec/3600;
	else                lon = lon_deg + lon_min/60 + lon_sec/3600;      end
	set(handles.edit_LonDec,'String',sprintf('%.4f',lon))
	lat = str2double(get(handles.edit_LatDec,'String'));
	elev = str2double(get(handles.edit_Elev,'String'));
	date = str2double(get(handles.edit_DateDec,'String'));
	out = igrf_m(lon, lat, elev, date);
	set_field_boxes(out, handles)   % Fill the various mag field boxes

% -------------------------------------------------------------------------------------------------
function edit_LonMin_Callback(hObject, eventdata, handles)
	if (isnan(get(hObject,'String'))),    set(hObject,'String','0'); return; end
	lon_min = str2double(get(hObject,'String'));
	if (lon_min < 0 || lon_min >= 60), set(hObject,'String','0'); return; end
	lon_min = fix(lon_min);     % Ensure that the value is an integer
	lon_deg = str2double(get(handles.edit_LonDeg,'String'));
	lon_sec = str2double(get(handles.edit_LonSec,'String'));
	if (lon_deg < 0),   lon = lon_deg - lon_min/60 - lon_sec/3600;
	else                lon = lon_deg + lon_min/60 + lon_sec/3600;      end
	set(handles.edit_LonDec,'String',sprintf('%.4f',lon))
	lat = str2double(get(handles.edit_LatDec,'String'));
	elev = str2double(get(handles.edit_Elev,'String'));
	date = str2double(get(handles.edit_DateDec,'String'));
	dms = dec2deg(lon,'opt');       % Convert decimal degrees in deg,min,sec
	set(handles.edit_LonDeg,'String',sprintf('%.0f',dms(1)))
	set(handles.edit_LonMin,'String',sprintf('%.0f',dms(2)))
	set(handles.edit_LonSec,'String',sprintf('%.2f',dms(3)))
	out = igrf_m(lon, lat, elev, date);
	set_field_boxes(out, handles)   % Fill the various mag field boxes

% -------------------------------------------------------------------------------------------------
function edit_LonSec_Callback(hObject, eventdata, handles)
	if (isnan(get(hObject,'String'))),    set(hObject,'String','0'); return; end
	lon_sec = str2double(get(hObject,'String'));
	if (lon_sec < 0 || lon_sec >= 60), set(hObject,'String','0'); return; end
	lon_deg = str2double(get(handles.edit_LonDeg,'String'));
	lon_min = str2double(get(handles.edit_LonMin,'String'));
	if (lon_deg < 0),   lon = lon_deg - lon_min/60 - lon_sec/3600;
	else                lon = lon_deg + lon_min/60 + lon_sec/3600;      end
	set(handles.edit_LatDec,'String',sprintf('%.4f',lon))
	lat = str2double(get(handles.edit_LatDec,'String'));
	elev = str2double(get(handles.edit_Elev,'String'));
	date = str2double(get(handles.edit_DateDec,'String'));
	out = igrf_m(lon, lat, elev, date);
	set_field_boxes(out, handles)   % Fill the various mag field boxes

% -------------------------------------------------------------------------------------------------
function edit_LonDec_Callback(hObject, eventdata, handles)
	if (isnan(get(hObject,'String'))),    set(hObject,'String','-8'); return; end
	lon = str2double(get(hObject,'String'));
	lat = str2double(get(handles.edit_LatDec,'String'));
	elev = str2double(get(handles.edit_Elev,'String'));
	date = str2double(get(handles.edit_DateDec,'String'));
	out = igrf_m(lon, lat, elev, date);
	set_field_boxes(out, handles)   % Fill the various mag field boxes

% -------------------------------------------------------------------------------------------------
function edit_Elev_Callback(hObject, eventdata, handles)
	xx = get(hObject,'String');
	if (isnan(xx)),    set(hObject,'String','0'); return; end
	lon = str2double(get(handles.edit_LonDec,'String'));
	lat = str2double(get(handles.edit_LatDec,'String'));
	elev = str2double(xx);
	date = str2double(get(handles.edit_DateDec,'String'));
	dms = dec2deg(lon,'opt');       % Convert decimal degrees in deg,min,sec
	set(handles.edit_LonDeg,'String',sprintf('%.0f',dms(1)))
	set(handles.edit_LonMin,'String',sprintf('%.0f',dms(2)))
	set(handles.edit_LonSec,'String',sprintf('%.2f',dms(3)))
	out = igrf_m(lon, lat, elev, date);
	set_field_boxes(out, handles)   % Fill the various mag field boxes

% -------------------------------------------------------------------------------------------------
function edit_DateDD_Callback(hObject, eventdata, handles)
	if (isnan(get(hObject,'String'))),    set(hObject,'String','1'); return; end
	day = fix(str2double(get(hObject,'String')));     % Make sure the value is integer
	month = str2double(get(handles.edit_DateMM,'String'));
	year  = str2double(get(handles.edit_DateYY,'String'));
	if (day < 1 || day > 31),    day = 1;    set(hObject,'String','1');  end
	if (day > 28 && month == 2 && ~(~rem(year, 4) && rem(year, 100) ) || ~rem(year, 400))
        day = 28;   set(hObject,'String','28');
	elseif (day > 29 && month == 2 && (~rem(year, 4) && rem(year, 100) ) || ~rem(year, 400))
        day = 29;   set(hObject,'String','29');
	elseif (day == 31 && (month == 4 || month == 6 || month == 9 || month == 11))
        day = 30;   set(hObject,'String','30');
	end
	yd = dec_year(year,month,day);
	set(handles.edit_DateDec,'String',sprintf('%.4f',yd))
	lat = str2double(get(handles.edit_LatDec,'String'));
	lon = str2double(get(handles.edit_LonDec,'String'));
	elev = str2double(get(handles.edit_Elev,'String'));
	out = igrf_m(lon, lat, elev, yd);
	set_field_boxes(out, handles)   % Fill the various mag field boxes

% -------------------------------------------------------------------------------------------------
function edit_DateMM_Callback(hObject, eventdata, handles)
	if (isnan(get(hObject,'String'))),    set(hObject,'String','1'); return; end
	month = fix(str2double(get(hObject,'String')));     % Make sure the value is integer
	if (month < 1 || month > 12),    month = 1;  set(hObject,'String','1');  end
	year  = str2double(get(handles.edit_DateYY,'String'));
	day = str2double(get(handles.edit_DateDD,'String'));
	yd = dec_year(year,month,day);
	set(handles.edit_DateDec,'String',sprintf('%.4f',yd))
	lat = str2double(get(handles.edit_LatDec,'String'));
	lon = str2double(get(handles.edit_LonDec,'String'));
	elev = str2double(get(handles.edit_Elev,'String'));
	out = igrf_m(lon, lat, elev, yd);
	set_field_boxes(out, handles)   % Fill the various mag field boxes

% -------------------------------------------------------------------------------------------------
function edit_DateYY_Callback(hObject, eventdata, handles)
	if (isnan(get(hObject,'String'))),    set(hObject,'String','2004'); return; end
	year = fix(str2double(get(hObject,'String')));
	if (year < handles.start_stop_epoch(1) || year > handles.start_stop_epoch(2))
        errordlg('Date is outside the coefficients period','Error');
        set(hObject,'String','2004');    return;
	end
	day  = str2double(get(handles.edit_DateDD,'String'));
	month = str2double(get(handles.edit_DateMM,'String'));
	yd = dec_year(year,month,day);
	set(handles.edit_DateDec,'String',sprintf('%.4f',yd))
	lat = str2double(get(handles.edit_LatDec,'String'));
	lon = str2double(get(handles.edit_LonDec,'String'));
	elev = str2double(get(handles.edit_Elev,'String'));
	out = igrf_m(lon, lat, elev, yd);
	set_field_boxes(out, handles)   % Fill the various mag field boxes

% -------------------------------------------------------------------------------------------------
function edit_DateDec_Callback(hObject, eventdata, handles)
	if (isnan(get(hObject,'String'))),    set(hObject,'String','2005'); return; end
	date = str2double(get(hObject,'String'));
	if (date < handles.start_stop_epoch(1) || date > handles.start_stop_epoch(2))
        errordlg('Date is outside the coefficients period','Error');
        set(hObject,'String','2005');    return;
	end
	lat = str2double(get(handles.edit_LatDec,'String'));
	lon = str2double(get(handles.edit_LonDec,'String'));
	elev = str2double(get(handles.edit_Elev,'String'));
	out = igrf_m(lon, lat, elev, date);
	set_field_boxes(out, handles)   % Fill the various mag field boxes

% -------------------------------------------------------------------------------------------------
function checkbox_Option_H_Callback(hObject, eventdata, handles)
	if (get(hObject,'Value'))
        set(handles.edit_nHeaders,'Enable','on')
	else
        set(handles.edit_nHeaders,'String','', 'Enable','inactive')
	end

% -------------------------------------------------------------------------------------------------
function edit_nHeaders_Callback(hObject, eventdata, handles)
	xx = str2double(get(hObject,'String'));
	if (xx <= 0 || isnan(xx))
        set(hObject,'String','')
	else
        set(hObject,'String',sprintf('%.0f',xx))      % Prevent any "jokes"
	end

% -------------------------------------------------------------------------------------------------
function pushbutton_HelpInputFile_Callback(hObject, eventdata, handles)
	%message_win('create',sprintf('\t\t\tHelp on ''Input File'' '));
	message = sprintf(['These option is used when you have a file with postions where you want to ',...
        'compute the IGRF. The minimum necessary information in the file is two columns. One for longitudes ',...
        'and the other for latitudes. With these you can only ask for computing the IGRF values at that ',...
        'positions. However, for computing the IGRF I also need to know a date and elevation. Given that ',...
        'your file doesn''t have that info I''ll get it from the "Elevation" and "Date" fields above in ',...
        'this widow.\n If, on the other hand, the file has that information (that is, it has more than ',...
        '2 columns) I''ll use that information instead. But things are not so simple. There are a lot of ',...
        'different ways by which information may be stored in a file. Being impossible to guess them, you ',...
        'are presented with a window column selecter. It is than up to you to guive me the correct info.\n',...
        'Furthermore, dates exist in a multitude of formats, including text formats. Before calling the ',...
        'column selector I already tryied to guess if dates in text format were present. Recognized ',...
        'formats (that is, formats that will be converted to decimal years) are:\n',...
        '\tdd-mmm-yyyy\tmmm-dd-yyyy\tdd-mm-yyyy\tmm-dd-yyyy\n',...
        '\tdd/mmm/yyyy\tmmm/dd/yyyy\tdd/mm/yyyy\tmm/dd/yyyy\n',...
        '\tdd:mmm:yyyy\tmmm:dd:yyyy\tdd:mm:yyyy\tmm:dd:yyyy\n',...
        'Where yyyy stands for year (e.g. 2004)\t mm means month in the 1-12 range (the same for days) ',...
        'and mmm is a string with month names (e.g. Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec)\n',...
        'Two checkboxes in the column selector window let you chose what to add to the output file. ',...
        'For each of the selected option, an extra column will be appended to the output file.']);
	message_win('create',message);

% -------------------------------------------------------------------------------------------------
function edit_InputFile_Callback(hObject, eventdata, handles)
	fname = get(hObject,'String');
	if isempty(fname)               % Reset the file computing section
        set(handles.edit_OutputFile,'String','')
        set(handles.checkbox_Option_H,'Value',0)
        set(handles.edit_nHeaders,'String','','Enable','inactive')
        handles.input_file = [];    handles.ind_igrf = [];
        return;
	end
	% Let the pushbutton_InputFile_Callback do all the work
	pushbutton_InputFile_Callback(hObject,[],guidata(gcbo),fname)

% -------------------------------------------------------------------------------------------------
function pushbutton_InputFile_Callback(hObject, eventdata, handles,opt)
	if (nargin == 4),   fname = opt;     end

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

	hFig = gcf;

	% See if number of headers are requested
	n_headers = [];
	if (get(handles.checkbox_Option_H,'Value'))
        xx = get(handles.edit_nHeaders,'String');
        if (~isempty(xx))
            n_headers = str2double(xx);
        end
	end

	set(hFig,'Name',['Reading file:' fname])
	if (isempty(n_headers))
        [numeric_data,date,headerlines,str_col] = text_read(fname,NaN);
	else
        [numeric_data,date,headerlines,str_col] = text_read(fname,NaN,n_headers);
	end
	set(hFig,'Name','IGRF Calculator')

	if (isempty(numeric_data))      % Even if the rest is not empty we have to quit and return
        errordlg('File doesn''t have any recognized nymeric data (Quiting).','Error');
        set(handles.checkbox_Option_H,'Value',0)
        set(handles.checkbox_nHeaders,'String','','Enable','inactive')
        return
	end

	% If msgbox exist we have to move it from behind the main window. So get it's handle
	hMsgFig = gcf;
	if (hFig ~= hMsgFig)
        uistack(hMsgFig,'top');    % If error msgbox exists, bring it forward
        % Here we have a stupid problem. If don't kill the message window before the
        % select_cols is called this later wont work. FDS I have no more patiente for this.
        pause(1)
        try     delete(hMsgFig);    end
	end

	if (~isempty(str_col))      % Keep the ralative position of the date column
        numeric_data = [numeric_data(:,1:str_col-1) date numeric_data(:,str_col:end)];
        out = select_cols(numeric_data,'igrf',fname,1000);
	else
        out = select_cols(numeric_data,'igrf',fname,1000);
	end

    if (isempty(out)),   return;    end

	if (headerlines)        % Update headers info in checkbox and editbox
        set(handles.checkbox_Option_H,'Value',1)
        set(handles.edit_nHeaders,'String',sprintf('%d',headerlines),'Enable','on')
	end

	[m,n] = size(numeric_data);
	ind_igrf = out(1:5);    ind_opt = out(6:7);     % Split indices to more convinient vars
	if (ind_igrf(3) == 0)   % Total Field is missing. OK it only means that we cannot compute the anomaly
        ind_opt(2) = 0;
	end
	if (ind_igrf(4) == 0)   % Date is missing (get it from the year edit box)
        date = str2double(get(handles.edit_DateDec,'String'));
        date = repmat(date,m,1);
	else
        date = numeric_data(:,ind_igrf(4));
	end
	if (ind_igrf(5) == 0)   % Elevation is missing (get it from the Elevation edit box)
        elev = str2double(get(handles.edit_Elev,'String'));
        elev = repmat(elev,m,1);
	else
        elev = numeric_data(:,ind_igrf(5));
	end

	% Because it's very hard to deal with all possibilities of column positions in the
	% array and at the same time allways keep track on what they contain, I'll add the
	% elev & date columns to the end of numeric_data. Like this I know exactly where they
	% are and can remove them before writing the result in a file in the pushbutton_ComputeFile
	numeric_data = [numeric_data elev date];

	% Save numeric_data in handles
	handles.input_file = numeric_data;
	handles.ind_igrf = out;     % We'll need to know the indices inside pushbutton_ComputeFile

	set(handles.edit_InputFile,'String',fname);
	[pathstr,name] = fileparts(fname);
	fname = [pathstr filesep name '_igrf.dat'];     % Proposed output filename
	set(handles.edit_OutputFile,'String',fname)
	guidata(hObject,handles)

% ------------------------
function mon = monstr2monnum(monstr)
	M = ['jan'; 'feb'; 'mar'; 'apr'; 'may'; 'jun'; 'jul'; 'aug'; 'sep'; 'oct'; 'nov'; 'dec'];
	m = size(monstr,1);
	mon = zeros(m,1);
	for i=1:m
        mon(i) = find(all((M == monstr(ones(12,1),1:3))'));
	end

% -------------------------------------------------------------------------------------------------
function pushbutton_OutputFile_Callback(hObject, eventdata, handles)
	if (~isempty(handles.h_calling_fig))                    % If we know the handle to the calling fig
        cfig_handles = guidata(handles.h_calling_fig);      % get handles of the calling fig
        last_dir = cfig_handles.last_dir;
        home = cfig_handles.home_dir;
	else
        last_dir = [];
	end
	if (~isempty(last_dir)),    cd(last_dir);   end
	[FileName,PathName] = uiputfile({'*.dat;*.DAT', 'Mag file (*.dat,*.DAT)';'*.*', 'All Files (*.*)'},'Select file');
	pause(0.01);
	if (~isempty(last_dir)),    cd(home);   end
	if isequal(FileName,0);     return;     end
	fname = [PathName FileName];
	set(handles.edit_OutputFile,'String',fname)

% -------------------------------------------------------------------------------------------------
function pushbutton_ComputeFile_Callback(hObject, eventdata, handles)
	% Compute the IGRF and output the results into a file

	fname = get(handles.edit_OutputFile,'String');
	if (~isempty(fname) && ~isempty(handles.input_file))
        [m,n] = size(handles.input_file);       % Remember that n is in excess of 2
        ind_igrf = handles.ind_igrf;            % ind_igrf(6) tells if write Total Field and
                                                % ind_igrf(7) if write anomaly.    
        f = igrf_m(handles.input_file(:,ind_igrf(1)), handles.input_file(:,ind_igrf(2)),...
            handles.input_file(:,end-1), handles.input_file(:,end),'-Ft');
        
        % See what was selected for writing, but first rip the last 2 columns that are repeated
        handles.input_file = handles.input_file(:,1:end-2);
        [m,n_old] = size(handles.input_file);   n = n_old;
        if (ind_igrf(6))
            handles.input_file = [handles.input_file f];
            n = n + 1;
        end
        if (ind_igrf(7))    % Compute anomaly. Note, this relyies in that a column with T exists
            handles.input_file = [handles.input_file handles.input_file(:,ind_igrf(3))-f];
            n = n + 1;
        end
        % Now write the results in the file fname
        format = repmat('%f\t',1,n_old);
        format = [format repmat('%.0f\t',1,n-n_old)];   % At least those I know they don't need decimals
        format = format(1:end-2);               % We don't want the last '\t'
        try
            double2ascii(fname,handles.input_file,format);
            msgbox('File successefully writen on disk','Sorte');
        catch
            errordlg('There was an error (origin unknown) while writing file on disk','Error');
        end
        % Reset the file computing section
        set(handles.edit_OutputFile,'String','')
        set(handles.edit_InputFile,'String','')
        set(handles.checkbox_Option_H,'Value',0)
        set(handles.edit_nHeaders,'String','','Enable','inactive')
        handles.input_file = [];    handles.ind_igrf = [];
	end
	guidata(hObject,handles)

% -------------------------------------------------------------------------------------------------
function edit_Xmin_Callback(hObject, eventdata, handles)
	xx = get(hObject,'String');     val = test_dms(xx);
	if ~isempty(val)            % when dd:mm or dd:mm:ss was given
        x_min = 0;
        if str2double(val{1}) > 0
            for i = 1:length(val),  x_min = x_min + str2double(val{i}) / (60^(i-1));    end
        else
            for i = 1:length(val),  x_min = x_min - abs(str2double(val{i})) / (60^(i-1));   end
        end
        handles.x_min = x_min;
        if ~isempty(handles.x_max) && x_min >= handles.x_max
            errordlg('West Longitude >= East Longitude ','Error in Longitude limits')
            set(hObject,'String','');   guidata(hObject, handles);  return
        end
        nc = get(handles.edit_Ncols,'String');
        if ~isempty(handles.x_max) && ~isempty(nc)       % x_max and ncols boxes are filled
            % Compute Ncols, but first must recompute x_inc
            x_inc = ivan_the_terrible((handles.x_max - x_min),round(abs(str2double(nc))),1);
            xx = floor((handles.x_max - str2double(xx)) / (str2double(get(handles.edit_Xinc,'String')))+0.5) + 1;
            set(handles.edit_Xinc,'String',sprintf('%.12f',x_inc))
        elseif ~isempty(handles.x_max)      % x_max box is filled but ncol is not, so put to the default (100)
            x_inc = ivan_the_terrible((handles.x_max - x_min),100,1);
            set(handles.edit_Xinc,'String',sprintf('%.12f',x_inc))
            set(handles.edit_Ncols,'String','100')
        end
	else                % box is empty, so clear also x_inc and ncols
        set(handles.edit_Xinc,'String','');     set(handles.edit_Ncols,'String','');
        set(hObject,'String','');
	end
	guidata(hObject, handles);

% -------------------------------------------------------------------------------------------------
function edit_Xmax_Callback(hObject, eventdata, handles)
	xx = get(hObject,'String');     val = test_dms(xx);
	if ~isempty(val)
        x_max = 0;
        if str2double(val{1}) > 0
            for i = 1:length(val),  x_max = x_max + str2double(val{i}) / (60^(i-1));    end
        else
            for i = 1:length(val),  x_max = x_max - abs(str2double(val{i})) / (60^(i-1));   end
        end
        handles.x_max = x_max;
        if ~isempty(handles.x_min) && x_max <= handles.x_min 
            errordlg('East Longitude <= West Longitude','Error in Longitude limits')
            set(hObject,'String','');   guidata(hObject, handles);  return
        end
        nc = get(handles.edit_Ncols,'String');
        if ~isempty(handles.x_min) && ~isempty(nc)       % x_max and ncols boxes are filled
            % Compute Ncols, but first must recompute x_inc
            x_inc = ivan_the_terrible((x_max - handles.x_min),round(abs(str2double(nc))),1);
            xx = floor((handles.x_min - str2double(xx)) / (str2double(get(handles.edit_Xinc,'String')))+0.5) + 1;
            set(handles.edit_Xinc,'String',sprintf('%.12f',x_inc))
        elseif ~isempty(handles.x_min)      % x_min box is filled but ncol is not, so put to the default (100)
            x_inc = ivan_the_terrible((x_max - handles.x_min),100,1);
            set(handles.edit_Xinc,'String',sprintf('%.12f',x_inc))
            set(handles.edit_Ncols,'String','100')
        end
	else                % box is empty, so clear also x_inc and ncols
        set(handles.edit_Xinc,'String','');     set(handles.edit_Ncols,'String','');
        set(hObject,'String','');
	end
	guidata(hObject, handles);

% -------------------------------------------------------------------------------------------------
function edit_Xinc_Callback(hObject, eventdata, handles)
	dms = 0;
	xx = get(hObject,'String');     val = test_dms(xx);
	if isempty(val),    return;     end
	% If it survived then ...
	if (length(val) > 1),    dms = 1;      end         % inc given in dd:mm or dd:mm:ss format
	x_inc = 0;
	for (i = 1:length(val)),   x_inc = x_inc + str2double(val{i}) / (60^(i-1));    end
	if ~isempty(handles.x_min) && ~isempty(handles.x_max)
        % Make whatever x_inc given compatible with GMT_grd_RI_verify
        x_inc = ivan_the_terrible((handles.x_max - handles.x_min), x_inc,2);
        if ~dms         % case of decimal unities
            set(hObject,'String',sprintf('%.12f',x_inc))
            ncol = floor((handles.x_max - handles.x_min) / x_inc + 0.5) + 1;
        else            % inc was in dd:mm or dd:mm:ss format
            ncol = floor((handles.x_max - handles.x_min) / x_inc + 0.5) + 1;
            ddmm = dec2deg(x_inc);
            set(hObject,'String',ddmm)
        end
        set(handles.edit_Ncols,'String',sprintf('%d',ncol))
	end
	handles.dms_xinc = dms;     handles.x_inc = str2double(xx);
	guidata(hObject, handles);

% -------------------------------------------------------------------------------------------------
function edit_Ncols_Callback(hObject, eventdata, handles)
	xx = get(hObject,'String');
	if isnan(xx)          % Idiot user attempt. Reset ncols.
        set(hObject,'String',handles.ncols);    return;
	end
	if ~isempty(get(handles.edit_Xmin,'String')) && ~isempty(get(handles.edit_Xmax,'String')) && ...
            ~isempty(get(handles.edit_Xinc,'String')) && ~isempty(xx)
        x_inc = ivan_the_terrible((handles.x_max - handles.x_min),round(abs(str2double(xx))),1);
        if handles.dms_xinc        % x_inc was given in dd:mm:ss format
            ddmm = dec2deg(x_inc);
            set(handles.edit_Xinc,'String',ddmm)
        else                    % x_inc was given in decimal format
            set(handles.edit_Xinc,'String',sprintf('%.12f',x_inc));
        end
        handles.ncols = str2double(xx);
	end
	guidata(hObject, handles);

% -------------------------------------------------------------------------------------------------
function edit_Ymin_Callback(hObject, eventdata, handles)
	% Read value either in decimal or in the dd:mm or dd_mm:ss formats and do some tests
	xx = get(hObject,'String');     val = test_dms(xx);
	if ~isempty(val)
        y_min = 0;
        if str2double(val{1}) > 0
            for i = 1:length(val),  y_min = y_min + str2double(val{i}) / (60^(i-1));    end
        else
            for i = 1:length(val),  y_min = y_min - abs(str2double(val{i})) / (60^(i-1));   end
        end
        handles.y_min = y_min;
        if ~isempty(handles.y_max) && y_min >= handles.y_max
            errordlg('South Latitude >= North Latitude','Error in Latitude limits')
            set(hObject,'String','');   guidata(hObject, handles);  return
        end
        nr = get(handles.edit_Nrows,'String');
        if ~isempty(handles.y_max) && ~isempty(nr)       % y_max and nrows boxes are filled
            % Compute Nrows, but first must recompute y_inc
            y_inc = ivan_the_terrible((handles.y_max - y_min),round(abs(str2double(nr))),1);
            xx = floor((handles.y_max - str2double(xx)) / (str2double(get(handles.edit_Yinc,'String')))+0.5) + 1;
            set(handles.edit_Yinc,'String',sprintf('%.12f',y_inc))
        elseif ~isempty(handles.y_max)      % y_max box is filled but nrows is not, so put to the default (100)
            y_inc = ivan_the_terrible((handles.y_max - y_min),100,1);
            set(handles.edit_Yinc,'String',sprintf('%.12f',y_inc))
            set(handles.edit_Nrows,'String','100')
        end
	else                % box is empty, so clear also y_inc and nrows
        set(handles.edit_Yinc,'String','');     set(handles.edit_Nrows,'String','');
        set(hObject,'String','');
	end
	guidata(hObject, handles);

% -------------------------------------------------------------------------------------------------
function edit_Ymax_Callback(hObject, eventdata, handles)
	% Read value either in decimal or in the dd:mm or dd_mm:ss formats and do some tests
	xx = get(hObject,'String');     val = test_dms(xx);
	if ~isempty(val)
        y_max = 0;
        if str2double(val{1}) > 0
            for i = 1:length(val),  y_max = y_max + str2double(val{i}) / (60^(i-1));    end
        else
            for i = 1:length(val),  y_max = y_max - abs(str2double(val{i})) / (60^(i-1));   end
        end
        handles.y_max = y_max;
        if ~isempty(handles.y_min) && y_max <= handles.y_min 
            errordlg('North Latitude <= South Latitude','Error in Latitude limits')
            set(hObject,'String','');   guidata(hObject, handles);  return
        end
        nr = get(handles.edit_Nrows,'String');
        if ~isempty(handles.y_min) && ~isempty(nr)       % y_min and nrows boxes are filled
            % Compute Nrows, but first must recompute y_inc
            y_inc = ivan_the_terrible((y_max - handles.y_min),round(abs(str2double(nr))),1);
            xx = floor((handles.y_min - str2double(xx)) / (str2double(get(handles.edit_Yinc,'String')))+0.5) + 1;
            set(handles.edit_Yinc,'String',sprintf('%.12f',y_inc))
        elseif ~isempty(handles.y_min)      % y_min box is filled but nrows is not, so put to the default (100)
            y_inc = ivan_the_terrible((y_max - handles.y_min),100,1);
            set(handles.edit_Yinc,'String',sprintf('%.12f',y_inc))
            set(handles.edit_Nrows,'String','100')
        end
	else                % This box is empty, so clear also y_inc and nrows
        set(handles.edit_Yinc,'String','');     set(handles.edit_Nrows,'String','');
        set(hObject,'String','');
	end
	guidata(hObject, handles);

% -------------------------------------------------------------------------------------------------
function edit_Yinc_Callback(hObject, eventdata, handles)
	dms = 0;
	xx = get(hObject,'String');     val = test_dms(xx);
	if isempty(val)
        set(hObject, 'String', '');    return
	end
	% If it survived then ...
	if (length(val) > 1),    dms = 1;      end         % inc given in dd:mm or dd:mm:ss format
	y_inc = 0;
	for i = 1:length(val),   y_inc = y_inc + str2double(val{i}) / (60^(i-1));    end
	if ~isempty(handles.y_min) && ~isempty(handles.y_max)
        % Make whatever y_inc given compatible with GMT_grd_RI_verify
        y_inc = ivan_the_terrible((handles.y_max - handles.y_min), y_inc,2);
        if ~dms         % case of decimal unities
            set(hObject,'String',sprintf('%.12f',y_inc))
            nrow = floor((handles.y_max - handles.y_min) / y_inc + 0.5) + 1;
        else            % inc was in dd:mm or dd:mm:ss format
            nrow = floor((handles.y_max - handles.y_min) / y_inc + 0.5) + 1;
            ddmm = dec2deg(y_inc);
            set(hObject,'String',ddmm)
        end
        set(handles.edit_Nrows,'String',sprintf('%d',nrow))
	end
	handles.dms_yinc = dms;     handles.y_inc = str2double(xx);
	guidata(hObject, handles);

% -------------------------------------------------------------------------------------------------
function edit_Nrows_Callback(hObject, eventdata, handles)
	xx = get(hObject,'String');
	if isnan(xx)          % Idiot user attempt. Reset nrows.
        set(hObject,'String',handles.nrows);    return;
	end
	if ~isempty(get(handles.edit_Ymin,'String')) && ~isempty(get(handles.edit_Ymax,'String')) && ...
            ~isempty(get(handles.edit_Yinc,'String'))
        y_inc = ivan_the_terrible((handles.y_max - handles.y_min),round(abs(str2double(xx))),1);
        if handles.dms_yinc        % y_inc was given in dd:mm:ss format
            ddmm = dec2deg(y_inc);
            set(handles.edit_Yinc,'String',ddmm)
        else                    % y_inc was given in decimal format
            set(handles.edit_Yinc,'String',sprintf('%.12f',y_inc));
        end
        handles.nrows = str2double(xx);
	end
	guidata(hObject, handles);

% -------------------------------------------------------------------------------------------------
function pushbutton_ComputeGrid_Callback(hObject, eventdata, handles)
	% Compute a grid of the field specified on popup_FieldComponent, but first do some tests
	xmin = get(handles.edit_Xmin,'String');     xmax = get(handles.edit_Xmax,'String');
	xinc = get(handles.edit_Xinc,'String');
	ymin = get(handles.edit_Ymin,'String');     ymax = get(handles.edit_Ymax,'String');
	yinc = get(handles.edit_Yinc,'String');
	if (isempty(xmin))
        errordlg('West grid limit is empty. Give a plausible value','Error');    return
	end
	if (isempty(xmax))
        errordlg('East grid limit is empty. Give a plausible value','Error');    return
	end
	if (isempty(xinc))
        errordlg('Xinc box is empty. Give a reasonable value','Error');    return
	end
	if (isempty(ymin))
        errordlg('South grid limit is empty. Give a plausible value','Error');    return
	end
	if (isempty(ymax))
        errordlg('North grid limit is empty. Give a plausible value','Error');    return
	end
	if (isempty(yinc))
        errordlg('Yinc box is empty. Give a reasonable value','Error');    return
	end
	elev = str2double(get(handles.edit_Elev,'String'));
	date = str2double(get(handles.edit_DateDec,'String'));
	xmin = str2double(xmin);    xmax = str2double(xmax);    xinc = str2double(xinc);
	ymin = str2double(ymin);    ymax = str2double(ymax);    yinc = str2double(yinc);

	% OK, now find out what to compute
	contents = get(handles.popup_FieldComponent,'String');
	switch contents{get(handles.popup_FieldComponent,'Value')}
        case 'Total field'
            opt_F = '-Ft';  name = 'IGRF Total Field';
        case 'H component'
            opt_F = '-Fh';  name = 'IGRF Horizontal component';
        case 'X component'
            opt_F = '-Fx';  name = 'IGRF X component';
        case 'Y component'
            opt_F = '-Fy';  name = 'IGRF Y component';
        case 'Y component'
            opt_F = '-Fz';  name = 'IGRF Z component';
        case 'Declination'
            opt_F = '-Fd';  name = 'IGRF Declination';
        case 'Inclination'
            opt_F = '-Fi';  name = 'IGRF Inclination';
	end

	X = xmin:xinc:xmax;     Y = ymin:yinc:ymax;
	m = length(Y);          n = length(X);
	XX = reshape(repmat(X,m,1),1,m*n);
	YY = repmat(Y,1,n);
	f = igrf_m(XX,YY,elev, date,opt_F);         clear XX YY;
	Zmin = min(f(:));    Zmax = max(f(:));
	f = single(reshape(f,m,n));

	tmp.head = [X(1) X(end) Y(1) Y(end) Zmin Zmax 0 xinc yinc];
	tmp.X = X;    tmp.Y = Y;    tmp.name = name;
	mirone(f,tmp);

% -------------------------------------------------------------------------------------------------
function pushbutton_HelpGrid_Callback(hObject, eventdata, handles)

%-----------------------------------------------------------------------------
function set_field_boxes(out, handles)
	% Fill the various mag field texts
	set(handles.text_FnT,'String',[sprintf('%.0f',out(1)) '  nT'])
	set(handles.text_HnT,'String',[sprintf('%.0f',out(2)) '  nT'])
	set(handles.text_XnT,'String',[sprintf('%.0f',out(3)) '  nT'])
	set(handles.text_YnT,'String',[sprintf('%.0f',out(4)) '  nT'])
	set(handles.text_ZnT,'String',[sprintf('%.0f',out(5)) '  nT'])
	set(handles.text_Ddeg,'String',[sprintf('%.1f',out(6)) ' º'])
	set(handles.text_Ideg,'String',[sprintf('%.1f',out(7)) ' º'])


% --- Creates and returns a handle to the GUI figure. 
function igrf_options_LayoutFcn(h1,handles)

set(h1,'PaperUnits',get(0,'defaultfigurePaperUnits'),...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'units','pixels',...
'MenuBar','none',...
'Name','IGRF Calculator',...
'NumberTitle','off',...
'Position',[520 129 480 671],...
'Renderer',get(0,'defaultfigureRenderer'),...
'RendererMode','manual',...
'Resize','off',...
'Tag','figure1',...
'Visible','off');

uicontrol('Parent',h1,'Position',[10 224 461 75],'FontSize',9,'Style','frame','Tag','frame4');
uicontrol('Parent',h1,'Position',[10 319 461 111],'Style','frame','Tag','frame3');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@igrf_options_uicallback,h1,'edit_LatDeg_Callback'},...
'Position',[70 394 41 21],...
'String','37',...
'Style','edit',...
'Tag','edit_LatDeg');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@igrf_options_uicallback,h1,'edit_LatMin_Callback'},...
'Position',[159 394 41 21],...
'String','0',...
'Style','edit',...
'Tag','edit_LatMin');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@igrf_options_uicallback,h1,'edit_LatSec_Callback'},...
'Position',[239 394 51 21],...
'String','0',...
'Style','edit',...
'Tag','edit_LatSec');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@igrf_options_uicallback,h1,'edit_LatDec_Callback'},...
'Position',[359 394 61 21],...
'String','37.0',...
'Style','edit',...
'Tag','edit_LatDec');

uicontrol('Parent',h1,'Position',[20 398 50 15],...
'String','Latitude',...
'Style','text');

uicontrol('Parent',h1,'Position',[20 366 50 15],...
'String','Longitude',...
'Style','text');

uicontrol('Parent',h1,'Position',[117 398 21 15],...
'String','deg',...
'Style','text');

uicontrol('Parent',h1,'Position',[206 398 21 15],...
'String','min','Style','text');

uicontrol('Parent',h1,'Position',[297 398 21 15],...
'String','sec',...
'Style','text');

uicontrol('Parent',h1,'Position',[426 397 41 15],...
'String','Decimal','Style','text');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@igrf_options_uicallback,h1,'edit_LonDeg_Callback'},...
'Position',[70 363 41 21],...
'String','-8',...
'Style','edit',...
'Tag','edit_LonDeg');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@igrf_options_uicallback,h1,'edit_LonMin_Callback'},...
'Position',[159 363 41 21],...
'String','0',...
'Style','edit',...
'Tag','edit_LonMin');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@igrf_options_uicallback,h1,'edit_LonSec_Callback'},...
'Position',[239 363 51 21],...
'String','0',...
'Style','edit',...
'Tag','edit_LonSec');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@igrf_options_uicallback,h1,'edit_LonDec_Callback'},...
'Position',[359 363 61 21],...
'String','-8.0',...
'Style','edit',...
'Tag','edit_LonDec');

uicontrol('Parent',h1,'Position',[117 367 21 15],...
'String','deg',...
'Style','text');

uicontrol('Parent',h1,'Position',[206 367 21 15],...
'String','min',...
'Style','text');

uicontrol('Parent',h1,'Position',[297 367 21 15],...
'String','sec',...
'Style','text');

uicontrol('Parent',h1,'Position',[426 366 41 15],...
'String','Decimal',...
'Style','text');

uicontrol('Parent',h1,'Position',[20 333 50 15],...
'String','Elevation',...
'Style','text');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@igrf_options_uicallback,h1,'edit_Elev_Callback'},...
'Position',[70 330 41 21],...
'String','0',...
'Style','edit',...
'Tag','edit_Elev');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@igrf_options_uicallback,h1,'edit_DateDD_Callback'},...
'Position',[227 330 25 21],...
'Style','edit',...
'TooltipString','day',...
'Tag','edit_DateDD');

uicontrol('Parent',h1,...
'HorizontalAlignment','left',...
'Position',[142 334 81 15],...
'String','Date (d/m/yyyy)',...
'Style','text');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@igrf_options_uicallback,h1,'edit_DateMM_Callback'},...
'Position',[268 330 25 21],...
'Style','edit',...
'TooltipString','month',...
'Tag','edit_DateMM');

uicontrol('Parent',h1,'Position',[255 333 11 15],...
'String','/','Style','text');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@igrf_options_uicallback,h1,'edit_DateYY_Callback'},...
'Position',[308 330 35 21],...
'Style','edit',...
'TooltipString','year',...
'Tag','edit_DateYY');

uicontrol('Parent',h1,'Position',[295 333 11 15],...
'String','/','Style','text');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@igrf_options_uicallback,h1,'edit_DateDec_Callback'},...
'Position',[360 330 61 21],...
'Style','edit',...
'Tag','edit_DateDec');

uicontrol('Parent',h1,'Position',[426 333 41 15],...
'HorizontalAlignment','left',...
'String','Decimal',...
'Style','text');

uicontrol('Parent',h1,...
'FontSize',9,...
'FontWeight','bold',...
'HorizontalAlignment','left',...
'Position',[17 267 70 16],...
'String','Total Field',...
'Style','text',...
'Tag','text_F');

uicontrol('Parent',h1,...
'FontSize',9,...
'HorizontalAlignment','left',...
'Position',[86 269 60 15],...
'String','nT',...
'Style','text',...
'Tag','text_FnT');

uicontrol('Parent',h1,...
'FontSize',10,...
'FontWeight','bold',...
'HorizontalAlignment','left',...
'Position',[158 269 19 16],...
'String','Inc',...
'Style','text',...
'Tag','text_I');

uicontrol('Parent',h1,...
'HorizontalAlignment','left',...
'Position',[188 270 60 15],...
'String','º',...
'Style','text',...
'Tag','text_Ideg');

uicontrol('Parent',h1,'Position',[259 268 30 17],...
'FontSize',10,...
'FontWeight','bold',...
'HorizontalAlignment','left',...
'String','Dec',...
'Style','text',...
'Tag','text_D');

uicontrol('Parent',h1,...
'FontSize',9,...
'HorizontalAlignment','left',...
'Position',[300 270 60 15],...
'String','º',...
'Style','text',...
'Tag','text_Ddeg');

uicontrol('Parent',h1,...
'FontSize',10,...
'FontWeight','bold',...
'HorizontalAlignment','right',...
'Position',[55 238 15 16],...
'String','X',...
'Style','text',...
'Tag','text_X');

uicontrol('Parent',h1,...
'FontSize',9,...
'HorizontalAlignment','left',...
'Position',[86 239 60 15],...
'String','nT',...
'Style','text',...
'Tag','text_XnT');

uicontrol('Parent',h1,...
'FontSize',10,...
'FontWeight','bold',...
'HorizontalAlignment','left',...
'Position',[158 239 15 15],...
'String','Y',...
'Style','text',...
'Tag','text_Y');

uicontrol('Parent',h1,...
'FontSize',9,...
'HorizontalAlignment','left',...
'Position',[182 239 60 15],...
'String','nT',...
'Style','text',...
'Tag','text_YnT');

uicontrol('Parent',h1,...
'FontSize',10,...
'FontWeight','bold',...
'HorizontalAlignment','left',...
'Position',[260 239 15 15],...
'String','Z',...
'Style','text',...
'Tag','text_Z');

uicontrol('Parent',h1,...
'FontSize',9,...
'HorizontalAlignment','left',...
'Position',[286 239 60 15],...
'String','nT',...
'Style','text',...
'Tag','text_ZnT');

uicontrol('Parent',h1,...
'FontSize',10,...
'FontWeight','bold',...
'HorizontalAlignment','left',...
'Position',[359 239 15 15],...
'String','H',...
'Style','text',...
'Tag','text_H');

uicontrol('Parent',h1,...
'FontSize',9,...
'HorizontalAlignment','left',...
'Position',[380 239 60 15],...
'String','nT',...
'Style','text',...
'Tag','text_HnT');

uicontrol('Parent',h1,'Position',[114 334 12 15],...
'HorizontalAlignment','left',...
'String','m',...
'Style','text');

uicontrol('Parent',h1,'Position',[10 13 461 93],...
'Enable','inactive','Style','frame','Tag','frame1');

uicontrol('Parent',h1,...
'Callback',{@igrf_options_uicallback,h1,'checkbox_Option_H_Callback'},...
'Position',[20 79 70 15],...
'String','Headers?',...
'Style','checkbox',...
'TooltipString','Are there any header lines in the input file?',...
'Tag','checkbox_Option_H');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@igrf_options_uicallback,h1,'edit_nHeaders_Callback'},...
'Enable','inactive',...
'HorizontalAlignment','left',...
'Position',[181 77 31 20],...
'Style','edit',...
'TooltipString','How many?',...
'Tag','edit_nHeaders');

uicontrol('Parent',h1,...
'BackgroundColor',[0.831372559070587 0.815686285495758 0.7843137383461],...
'Callback',{@igrf_options_uicallback,h1,'pushbutton_HelpInputFile_Callback'},...
'FontWeight','bold',...
'ForegroundColor',[0 0 1],...
'Position',[427 54 31 21],...
'String','?',...
'Tag','pushbutton_HelpInputFile');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@igrf_options_uicallback,h1,'edit_InputFile_Callback'},...
'HorizontalAlignment','left',...
'Position',[20 51 341 22],...
'Style','edit',...
'TooltipString','Enter a file with total field measurements',...
'Tag','edit_InputFile');

uicontrol('Parent',h1,...
'Callback',{@igrf_options_uicallback4,h1,[],'pushbutton_InputFile_Callback'},...
'Position',[359 50 23 23],...
'Tag','pushbutton_InputFile');

uicontrol('Parent',h1,'Position',[31 98 75 15],...
'Enable','inactive',...
'String','Input Mag File',...
'Style','text',...
'Tag','text31');

uicontrol('Parent',h1,...
'HorizontalAlignment','left',...
'Position',[110 79 67 15],...
'String','Nº of headers',...
'Style','text',...
'TooltipString','How many?',...
'Tag','text32');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'HorizontalAlignment','left',...
'Position',[20 23 341 22],...
'Style','edit',...
'TooltipString','Output File Name for the cumputed field',...
'Tag','edit_OutputFile');

uicontrol('Parent',h1,...
'Callback',{@igrf_options_uicallback,h1,'pushbutton_OutputFile_Callback'},...
'Position',[359 22 23 23],...
'Tag','pushbutton_OutputFile');

uicontrol('Parent',h1,'Position',[10 125 461 81],...
'Enable','inactive','Style','frame','Tag','frame2');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@igrf_options_uicallback,h1,'edit_Xmin_Callback'},...
'Position',[76 164 71 21],...
'Style','edit',...
'String', '-180',...
'TooltipString','X min value',...
'Tag','edit_Xmin');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Callback',{@igrf_options_uicallback,h1,'edit_Xmax_Callback'},...
'Position',[152 164 71 21],...
'Style','edit',...
'String', '180',...
'TooltipString','X max value',...
'Tag','edit_Xmax');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@igrf_options_uicallback,h1,'edit_Xinc_Callback'},...
'Position',[228 164 71 21],...
'Style','edit',...
'TooltipString','DX grid spacing',...
'Tag','edit_Xinc');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@igrf_options_uicallback,h1,'edit_Ncols_Callback'},...
'Position',[304 164 45 21],...
'Style','edit',...
'TooltipString','Number of columns in the grid',...
'Tag','edit_Ncols');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@igrf_options_uicallback,h1,'edit_Ymin_Callback'},...
'Position',[76 138 71 21],...
'Style','edit',...
'TooltipString','Y min value',...
'Tag','edit_Ymin');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@igrf_options_uicallback,h1,'edit_Ymax_Callback'},...
'Position',[152 138 71 21],...
'Style','edit',...
'TooltipString','Y max value',...
'Tag','edit_Ymax');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@igrf_options_uicallback,h1,'edit_Yinc_Callback'},...
'Position',[228 138 71 21],...
'Style','edit',...
'TooltipString','DY grid spacing',...
'Tag','edit_Yinc');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@igrf_options_uicallback,h1,'edit_Nrows_Callback'},...
'Position',[304 138 45 21],...
'Style','edit',...
'TooltipString','Number of columns in the grid',...
'Tag','edit_Nrows');

uicontrol('Parent',h1,...
'BackgroundColor',[0.831372559070587 0.815686285495758 0.7843137383461],...
'Callback',{@igrf_options_uicallback,h1,'pushbutton_HelpGrid_Callback'},...
'FontWeight','bold',...
'ForegroundColor',[0 0 1],...
'Position',[369 135 31 23],...
'String','?',...
'Tag','pushbutton_HelpGrid');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[18 169 55 15],...
'String','X Direction',...
'Style','text');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[17 143 55 15],...
'String','Y Direction',...
'Style','text');

uicontrol('Parent',h1,'Position',[169 186 41 13],...
'Enable','inactive',...
'String','Max',...
'Style','text');

uicontrol('Parent',h1,'Position',[91 187 41 13],...
'Enable','inactive',...
'String','Min',...
'Style','text');

uicontrol('Parent',h1,...
'Enable','inactive',...
'Position',[246 187 41 13],...
'String','Spacing',...
'Style','text');

uicontrol('Parent',h1,'Position',[302 187 51 13],...
'Enable','inactive',...
'String','# of lines',...
'Style','text');

uicontrol('Parent',h1,'Position',[27 198 110 15],...
'Enable','inactive',...
'String','Griding Line Geometry',...
'Style','text',...
'Tag','text39');

uicontrol('Parent',h1,...
'Callback',{@igrf_options_uicallback,h1,'pushbutton_ComputeFile_Callback'},...
'FontWeight','bold',...
'Position',[404 22 56 23],...
'String','Compute',...
'Tag','pushbutton_ComputeFile');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Position',[370 163 91 22],...
'String',{'Total field'; 'H component'; 'X component'; 'Y component'; 'Z component'; 'Declination'; 'Inclination' },...
'Style','popupmenu',...
'Value',1,...
'Tag','popup_FieldComponent');

uicontrol('Parent',h1,...
'Callback',{@igrf_options_uicallback,h1,'pushbutton_ComputeGrid_Callback'},...
'FontWeight','bold',...
'Position',[405 135 56 23],...
'String','Compute',...
'Tag','pushbutton_ComputeGrid');

uicontrol('Parent',h1,'Position',[31 290 65 15],...
'Enable','inactive',...
'String','Point values',...
'Style','text',...
'Tag','text40');

uicontrol('Parent',h1,'Position',[30 422 51 15],...
'Enable','inactive',...
'String','Location',...
'Style','text',...
'Tag','text41');

axes('Parent',h1,'Units','pixels',...
'FontSize',9,...
'Position',[27 457 441 211],...
'Tag','axes1');

function igrf_options_uicallback(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
feval(callback_name,hObject,[],guidata(h1));

function igrf_options_uicallback4(hObject, eventdata, h1, opt, callback_name)
% This function is executed by the callback and than the handles is allways updated.
feval(callback_name,hObject,[],guidata(h1),opt);

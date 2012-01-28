function varargout = read_FlatFile(varargin)
% Helper window to read a flat naked true raw file

%	Copyright (c) 2004-2012 by J. Luis
%
% 	This program is part of Mirone and is free software; you can redistribute
% 	it and/or modify it under the terms of the GNU Lesser General Public
% 	License as published by the Free Software Foundation; either
% 	version 2.1 of the License, or any later version.
% 
% 	This program is distributed in the hope that it will be useful,
% 	but WITHOUT ANY WARRANTY; without even the implied warranty of
% 	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% 	Lesser General Public License for more details.
%
%	Contact info: w3.ualg.pt/~jluis/mirone
% --------------------------------------------------------------------
 
	hObject = figure('Vis','off');
	read_FlatFile_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hObject,'center')

	handles.interleave = '';

	if (~isempty(varargin))
		handles.fname = varargin{1}{1};
		[pathname,fname,ext] = fileparts(handles.fname);
		set(hObject,'Name',['Read Flat File (' fname '.' ext ')'])
	end

	tip = sprintf('%s\n%s\n%s','Specifies that the bands within the file are stored',...
			'in Band Sequential (BSQ) format. Landsat files',...
			'are normally supplied in this format.');
	set(handles.radio_BSQ,'TooltipString',tip)
	tip = sprintf('%s\n%s\n%s\n%s','Specifies that the bands within the file are stored',...
			'in Band Interleaved by Line (BIL) format. This means',...
			'that each line of pixels is stored one band after the',...
			'next within the file.');
	set(handles.radio_BIL,'TooltipString',tip)
	tip = sprintf('%s\n%s\n%s','Specifies that the bands within the file are stored',...
			'in Band Interleaved by Pixel (BIP) format. Landsat',...
			'MSS files are occasionally supplied in this format.');
	set(handles.radio_BIP,'TooltipString',tip)

	%------------ Give a Pro look (3D) to the frame boxes  -------------------------------
	new_frame3D(hObject, [handles.text_pixelForm handles.text_interleaveForm handles.text_optional])
	%------------- END Pro look (3D) -----------------------------------------------------

	% Choose default command line output for read_FlatFile_export
	handles.output = cell(3,1);
	guidata(hObject, handles);

	set(hObject,'Vis','on');
	% UIWAIT makes read_FlatFile_export wait for user response (see UIRESUME)
	uiwait(handles.figure1);

	handles = guidata(hObject);
	varargout = handles.output;
	delete(handles.figure1)

% ------------------------------------------------------------------------------------------
function edit_nHeadeBytes_CB(hObject, handles)
	hdr = str2double(get(hObject,'String'));
	if (isnan(hdr) || hdr < 0),    set(hObject,'String','0');   end

% ------------------------------------------------------------------------------------------
function edit_pixelsPerLine_CB(hObject, handles)
	n = str2double(get(hObject,'String'));
	if (isnan(n) || n < 0),		set(hObject,'String','')
	else						set(handles.edit_Xlast,'String',get(hObject,'String'))
	end

% ------------------------------------------------------------------------------------------
function edit_nLines_CB(hObject, handles)
	n = str2double(get(hObject,'String'));
	if (isnan(n) || n < 0),		set(hObject,'String','')
	else						set(handles.edit_Ylast,'String',get(hObject,'String'))
	end

% ------------------------------------------------------------------------------------------
function edit_nBands_CB(hObject, handles)
	n = str2double(get(hObject,'String'));
	if (isnan(n) || n < 0),		set(hObject,'String','1')
	else						set(handles.edit_lastBand,'String',get(hObject,'String'))
	end

% ------------------------------------------------------------------------------------------
function edit_firstBand_CB(hObject, handles)
	n = str2double(get(hObject,'String'));
	if (isnan(n) || n < 0 || n > str2double(get(handles.edit_lastBand,'String')))
		set(hObject,'String','1')
	end

% ------------------------------------------------------------------------------------------
function edit_lastBand_CB(hObject, handles)
	n = str2double(get(hObject,'String'));
	if (isnan(n) || n < 0 || n > str2double(get(handles.edit_nBands,'String')))
		set(hObject,'String',get(handles.edit_nBands,'String'))
	end

% ------------------------------------------------------------------------------------------
function edit_Xfirst_CB(hObject, handles)
	n = str2double(get(hObject,'String'));
	if (isnan(n) || n < 0 || n >= str2double(get(handles.edit_pixelsPerLine,'String')))
		set(hObject,'String','1')
	end

% ------------------------------------------------------------------------------------------
function edit_Xlast_CB(hObject, handles)
	n = str2double(get(hObject,'String'));
	if (isnan(n) || n < 0 || n >= str2double(get(handles.edit_pixelsPerLine,'String')))
		set(hObject,'String',get(handles.edit_pixelsPerLine,'String'))
	end

% ------------------------------------------------------------------------------------------
function edit_Xsample_CB(hObject, handles)
	x = str2double(get(hObject,'String'));
	if (isnan(x) || x < 0),    set(hObject,'String','1');   end

% ------------------------------------------------------------------------------------------
function edit_Yfirst_CB(hObject, handles)
	n = str2double(get(hObject,'String'));
	if (isnan(n) || n < 0 || n >= str2double(get(handles.edit_nLines,'String')))
		set(hObject,'String','1')
	end

% ------------------------------------------------------------------------------------------
function edit_Ylast_CB(hObject, handles)
	n = str2double(get(hObject,'String'));
	if (isnan(n) || n < 0 || n >= str2double(get(handles.edit_nLines,'String')))
		set(hObject,'String',get(handles.edit_nLines,'String'))
	end

% ------------------------------------------------------------------------------------------
function edit_Ysample_CB(hObject, handles)
	x = str2double(get(hObject,'String'));
	if (isnan(x) || x < 0),    set(hObject,'String','1');   end

% ------------------------------------------------------------------------------------------
function radio_BSQ_CB(hObject, handles)
% Hint: get(hObject,'Value') returns toggle state of radio_BSQ
	if (get(hObject,'Value'))
		set(handles.radio_BIL,'Value',0)
		set(handles.radio_BIP,'Value',0)
		handles.interleave = 'bsq';
		guidata(handles.figure1,handles)
	else
		set(hObject,'Value',1)
	end

% ------------------------------------------------------------------------------------------
function radio_BIL_CB(hObject, handles)
	if (get(hObject,'Value'))
		set(handles.radio_BSQ,'Value',0)
		set(handles.radio_BIP,'Value',0)
		handles.interleave = 'bil';
		guidata(handles.figure1,handles)
	else
		set(hObject,'Value',1)
	end

% ------------------------------------------------------------------------------------------
function radio_BIP_CB(hObject, handles)
	if (get(hObject,'Value'))
		set(handles.radio_BSQ,'Value',0)
		set(handles.radio_BIL,'Value',0)
		handles.interleave = 'bip';
		guidata(handles.figure1,handles)
	else
		set(hObject,'Value',1)
	end

% ------------------------------------------------------------------------------------------
function push_OK_CB(hObject, handles)

% Get the skip header bytes
hdr = str2double(get(handles.edit_nHeadeBytes,'String'));

% Get the Pixels per line
n_col = str2double(get(handles.edit_pixelsPerLine,'String'));
if (isnan(n_col))
    errordlg('ERROR: Must inform me on the number of "Pixels per Line".','ERROR');  return
end

% Get the number of lines
n_row = str2double(get(handles.edit_nLines,'String'));
if (isnan(n_row))
	errordlg('ERROR: Must inform me on the number of "Number of Lines".','ERROR');  return
end

% Get number of bands
n_band = str2double(get(handles.edit_nBands,'String'));
if (isnan(n_band))
	errordlg('ERROR: Must inform me on the number of "Number of Bands".','ERROR');  return
end

% Initialize those
got_subset_row = 0;         got_subset_column = 0;          got_subset_band = 0;

% Test for First and Last bands
kf = str2double(get(handles.edit_firstBand,'String'));
kl = str2double(get(handles.edit_lastBand,'String'));
if (kf > kl)
    errordlg('IDIOT COICE OF First and Last bands parameters.','Chico Clever'); return
end
if (kf ~= 1 || kl ~= n_band)    % We have band subset request
    subset_band = {'Band','Direct',kf:kl};
    got_subset_band = 1;
end

% Test for a row subset request
krf = str2double(get(handles.edit_Yfirst,'String'));        % First row
krl = str2double(get(handles.edit_Ylast,'String'));         % Last row
krs = str2double(get(handles.edit_Ysample,'String'));       % Increment row
if ( (~isnan(krl) && krl ~= n_row) || krs ~= 1)
    subset_row = {'Row','Range',[krf krs krl]};
    got_subset_row = 1;
elseif (isnan(krl) && krs ~= 1)
    errordlg('ERROR: To select a ROW step, you must also give the last row.','ERROR')
    return
elseif (krl-krf < 10*krs)
    errordlg('ERROR: Improper choice of First and Last rows.','ERROR');     return
end

% Test for a column subset request
kcf = str2double(get(handles.edit_Xfirst,'String'));        % First column
kcl = str2double(get(handles.edit_Xlast,'String'));         % Last column
kcs = str2double(get(handles.edit_Xsample,'String'));       % Increment column
if ( (~isnan(kcl) && kcl ~= n_col) || kcs ~= 1)
    subset_column = {'Column','Range',[kcf kcs kcl]};
    got_subset_column = 1;
elseif (isnan(kcl) && kcs ~= 1)
    errordlg('ERROR: To select a COLUMN step, you must also give the last column.','ERROR')
    return
elseif (kcl-kcf < 10*kcs)
    errordlg('ERROR: Improper choice of First and Last columns.','ERROR');     return
end

% Test for the Interleave format
if (isempty(handles.interleave))
    errordlg('Must select the Interleave format','Error');    return
end

% Get the Pixel resolution
switch (get(handles.popup_dataType,'Value'))
    case 1,        fmt = '*uint8';
    case 2,        fmt = '*uint16';
    case 3,        fmt = '*int16';
    case 4,        fmt = '*uint32';
    case 5,        fmt = '*int32';
    case 6,        fmt = '*float';
end

% Get the endianess
endian = 'ieee-le';         % Default ot little-endian
if (get(handles.check_swapBytes,'Value'))
    endian = 'ieee-be';
end

% OK, do the reading. But shit we still need to test for 8 possibilities
fname = handles.fname;
dims = [n_row n_col n_band];
interl = handles.interleave;
if (got_subset_row && got_subset_column && got_subset_band)         % All three
    raw = multibandread_j(fname,dims,fmt,hdr,interl,endian,subset_row,subset_column,subset_band);
    third_out = [subset_row; subset_column; subset_band];
elseif (got_subset_row && got_subset_column && ~got_subset_band)    % Row & Column
    raw = multibandread_j(fname,dims,fmt,hdr,interl,endian,subset_row,subset_column);
    third_out = [subset_row; subset_column];
elseif (got_subset_row && got_subset_band && ~got_subset_column)    % Row & Band
    raw = multibandread_j(fname,dims,fmt,hdr,interl,endian,subset_row,subset_band);
    third_out = [subset_row; subset_band];
elseif (got_subset_column && got_subset_band && ~got_subset_row)    % Column & Band
    raw = multibandread_j(fname,dims,fmt,hdr,interl,endian,subset_column,subset_band);
    third_out = [subset_column; subset_band];
elseif (got_subset_row && ~got_subset_column && ~got_subset_band)   % Row only
    raw = multibandread_j(fname,dims,fmt,hdr,interl,endian,subset_row);
    third_out = subset_row;
elseif (got_subset_column && ~got_subset_row && ~got_subset_band)   % Column only
    raw = multibandread_j(fname,dims,fmt,hdr,interl,endian,subset_column);
    third_out = subset_column;
elseif (got_subset_band && ~got_subset_row && ~got_subset_column)   % Band only
    raw = multibandread_j(fname,dims,fmt,hdr,interl,endian,subset_band);
    third_out = subset_band;
else
    raw = multibandread_j(fname,dims,fmt,hdr,interl,endian);
    third_out = [];
end

if (~isa(raw,'uint8'))
    img = alloc_mex(size(raw),'uint8');     % Pre allocation
	for i=1:size(raw,3)
        img(:,:,i) = scaleto8(raw(:,:,i));
	end
    handles.output{1} = img;
else
    handles.output{1} = raw;
end

handles.output{2} = {fmt, hdr, handles.interleave, endian};
handles.output{3} = third_out;
guidata(hObject, handles);
uiresume(handles.figure1);

% ------------------------------------------------------------------------------------------
function push_cancel_CB(hObject, handles)
	if (exist('OCTAVE_VERSION','builtin'))		% To know if we are running under Octave
		do_uiresume = ( isprop(hObject, '__uiwait_state__') && strcmp(get(hObject, '__uiwait_state__'), 'active') );
	else
		do_uiresume = strcmp(get(hObject, 'waitstatus'), 'waiting');
	end
	if (do_uiresume)		% The GUI is still in UIWAIT, us UIRESUME
		uiresume(hObject);
	else					% The GUI is no longer waiting, just close it
		delete(handles.figure1)
	end

% ------------------------------------------------------------------------------------------
function figure1_CloseRequestFcn(hObject, evt)
	handles = guidata(hObject);
	push_cancel_CB(hObject, handles)

% ----------- Creates and returns a handle to the GUI figure. 
function read_FlatFile_LayoutFcn(h1)

set(h1,...
'CloseRequestFcn',@figure1_CloseRequestFcn,...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'MenuBar','none',...
'Name','Read Flat File',...
'NumberTitle','off',...
'Position',[520 544 504 206],...
'Resize','off',...
'Tag','figure1');

uicontrol('Parent',h1,'Position',[225 5 191 41],'Style','frame');
uicontrol('Parent',h1,'Position',[341 65 157 111],'Style','frame');
uicontrol('Parent',h1,'Position',[5 65 330 111],'Style','frame');
uicontrol('Parent',h1,'Position',[5 5 211 41],'Style','frame');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@flatfile_uiCB,...
'Position',[124 134 51 21],'String','0','Style','edit',...
'TooltipString','Skip this number of bytes (that is skip the header)',...
'Tag','edit_nHeadeBytes');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@flatfile_uiCB,...
'Position',[124 104 51 21],'Style','edit',...
'TooltipString','Number of columns in the dataset',...
'Tag','edit_pixelsPerLine');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@flatfile_uiCB,...
'Position',[124 74 51 21],'Style','edit',...
'TooltipString','Number of rows in the dataset',...
'Tag','edit_nLines');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@flatfile_uiCB,...
'Position',[292 134 31 21],'String','1','Style','edit',...
'TooltipString','Number of bands in the dataset',...
'Tag','edit_nBands');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@flatfile_uiCB,...
'Position',[292 104 31 21],'String','1','Style','edit',...
'TooltipString','First band to load from the dataset',...
'Tag','edit_firstBand');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@flatfile_uiCB,...
'Position',[292 74 31 21],'String','1','Style','edit',...
'TooltipString','Last band to load from the dataset',...
'Tag','edit_lastBand');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'FontName','Helvetica','Position',[11 12 101 22],...
'String',{'8-bit unsigned'; '16-bit unsigned'; '16-bit signed'; '32-bit unsigned'; '32-bit signed'; '32-bit float'},...
'Style','popupmenu',...
'TooltipString','If you don''t know what this is, you better study a little',...
'Value',1,'Tag','popup_dataType');

uicontrol('Parent',h1,...
'FontName','Helvetica','Position',[127 16 85 15],'String','Swap bytes',...
'Style','checkbox',...
'TooltipString','Swap byte order (default is Little endian. eg. Intel order)',...
'Tag','check_swapBytes');

uicontrol('Parent',h1,...
'Call',@flatfile_uiCB,...
'FontName','Helvetica','Position',[234 17 41 15],'String','BSQ',...
'Style','radiobutton','Tag','radio_BSQ');

uicontrol('Parent',h1,...
'Call',@flatfile_uiCB,...
'FontName','Helvetica','Position',[307 17 41 15],'String','BIL',...
'Style','radiobutton','Tag','radio_BIL');

uicontrol('Parent',h1,...
'Call',{@flatfile_uiCB,h1,'radio_BIP_CB'},...
'FontName','Helvetica','Position',[372 17 41 15],'String','BIP',...
'Style','radiobutton','Tag','radio_BIP');

uicontrol('Parent',h1,...
'Call',@flatfile_uiCB,...
'FontName','Helvetica','FontSize',9,'Position',[434 6 65 23],...
'String','OK','Tag','push_OK');

uicontrol('Parent',h1,...
'Call',@flatfile_uiCB,...
'FontName','Helvetica','FontSize',9,'Position',[434 35 65 23],...
'String','Cancel','Tag','push_cancel');

uicontrol('Parent',h1,'FontName','Helvetica','HorizontalAlignment','left',...
'Position',[17 137 105 15],'String','Header length (bytes)','Style','text');

uicontrol('Parent',h1,'FontName','Helvetica','HorizontalAlignment','left',...
'Position',[16 109 105 15],'String','Pixels per Line','Style','text');

uicontrol('Parent',h1,'FontName','Helvetica','HorizontalAlignment','left',...
'Position',[16 79 105 15],'String','Number of Lines','Style','text');

uicontrol('Parent',h1,'FontName','Helvetica','HorizontalAlignment','left',...
'Position',[205 138 85 15],'String','Number of Bands','Style','text');

uicontrol('Parent',h1,'FontName','Helvetica','HorizontalAlignment','left',...
'Position',[205 109 85 15],'String','First Band','Style','text');

uicontrol('Parent',h1,'FontName','Helvetica','HorizontalAlignment','left',...
'Position',[205 79 85 15],'String','Last Band','Style','text');

uicontrol('Parent',h1,'FontName','Helvetica','FontSize',9,...
'Position',[16 38 80 15],'String','Pixel Format','Style','text','Tag','text_pixelForm');

uicontrol('Parent',h1,'FontName','Helvetica','FontSize',9,...
'Position',[236 38 100 15],'String','Interleave Format',...
'Style','text','Tag','text_interleaveForm');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@flatfile_uiCB,...
'Position',[387 134 41 21],'String','1','Style','edit',...
'TooltipString','Start column index','Tag','edit_Xfirst');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',{@flatfile_uiCB,h1,'edit_Xlast_CB'},...
'Position',[387 104 41 21],'Style','edit',...
'TooltipString','Last column index','Tag','edit_Xlast');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@flatfile_uiCB,...
'Position',[387 74 41 21],'String','1','Style','edit',...
'TooltipString','Increment in columns (e.g. 2 = read every other column)',...
'Tag','edit_Xsample');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@flatfile_uiCB,...
'Position',[447 134 41 21],'String','1','Style','edit',...
'TooltipString','Start row index','Tag','edit_Yfirst');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@flatfile_uiCB,...
'Position',[447 104 41 21],'Style','edit',...
'TooltipString','Last row index','Tag','edit_Ylast');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@flatfile_uiCB,...
'Position',[447 74 41 21],'String','1','Style','edit',...
'TooltipString','Increment in rows (e.g. 2 = read every other row)',...
'Tag','edit_Ysample');

uicontrol('Parent',h1,'FontName','Helvetica',...
'HorizontalAlignment','left','Position',[351 138 35 15],...
'String','First','Style','text');

uicontrol('Parent',h1,'FontName','Helvetica',...
'HorizontalAlignment','left','Position',[352 108 35 15],...
'String','Last','Style','text');

uicontrol('Parent',h1,'FontName','Helvetica','HorizontalAlignment','left',...
'Position',[352 78 35 15],'String','Sample',...
'Style','text');

uicontrol('Parent',h1,'FontName','Helvetica','FontSize',9,...
'Position',[348 167 105 18],'String','Window (optional)',...
'Style','text','Tag','text_optional');

uicontrol('Parent',h1,'Position',[189 65 3 101], 'Style','frame');

uicontrol('Parent',h1,'FontName','Helvetica',...
'Position',[404 155 10 15],'String','X',...
'Style','text');

uicontrol('Parent',h1,'FontName','Helvetica',...
'Position',[463 155 10 15],'String','Y',...
'Style','text');

uicontrol('Parent',h1,'FontName','Helvetica',...
'FontSize',11,'Position',[91 186 254 17],...
'String','Please enter the following information',...
'Style','text','Tag','text_header');

function flatfile_uiCB(hObject, eventdata)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));

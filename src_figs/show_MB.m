function varargout = show_MB(varargin)
% Helper window to decide what to do when a multi-beam file, or datalist.mb-1, is loaded

%	Copyright (c) 2004-2017 by J. Luis
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

% $Id$

	hObject = figure('Vis','off');
	show_MB_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hObject,'center');

	handMir = varargin{1};
	handles.fnameMB = varargin{2};
	handles.whichFleder = handMir.whichFleder;
	handles.TDRver   = handMir.TDRver;
	handles.path_tmp = handMir.path_tmp;
	handles.no_file  = handMir.no_file;
	handles.hMirFig  = handMir.figure1;

	%------------ Give a Pro look (3D) to the frame boxes  ---------------------------------------------
	new_frame3D(hObject, NaN)
	%------------- END Pro look (3D) -------------------------------------------------------------------

	% Make this Fig be always at the top
	if (strncmp(computer,'PC',2) && handMir.version7 < 8.4)
		WindowAPI(hObject, 'TopMost')
	end

	guidata(hObject, handles);
	set(hObject,'Visible','on');
	if (nargout),   varargout{1} = hObject;     end

% ----------------------------------------------------------------------
function radio_pointCloud_CB(hObject, handles)
	if (~get(hObject,'Val')),	set(hObject,'Val',1),	return,		end
	set(handles.radio_image,'Value',0)
	set([handles.popup_symb handles.edit_symbSize],'Enable','on')
	set([handles.radio_imgSimple handles.radio_imgShaded handles.radio_imgAmp ...
		handles.radio_imgSS],'Enable','off')

% ----------------------------------------------------------------------
function radio_image_CB(hObject, handles)
	if (~get(hObject,'Val')),	set(hObject,'Val',1),	return,		end
	set(handles.radio_pointCloud,'Value',0)
	set([handles.radio_imgSimple handles.radio_imgShaded handles.radio_imgAmp ...
		handles.radio_imgSS],'Enable','on', 'Val',0)
	set(handles.radio_imgSimple, 'Val', 1)
	set([handles.popup_symb handles.edit_symbSize],'Enable','off')

% ----------------------------------------------------------------------
function radio_imgSimple_CB(hObject, handles)
	if (~get(hObject,'Val')),	set(hObject,'Val',1),	return,		end
	set([handles.radio_imgShaded handles.radio_imgAmp handles.radio_imgSS],'Value',0)

% ----------------------------------------------------------------------
function radio_imgShaded_CB(hObject, handles)
	if (~get(hObject,'Val')),	set(hObject,'Val',1),	return,		end
	set([handles.radio_imgSimple handles.radio_imgAmp handles.radio_imgSS],'Value',0)

% ----------------------------------------------------------------------
function radio_imgAmp_CB(hObject, handles)
	if (~get(hObject,'Val')),	set(hObject,'Val',1),	return,		end
	set([handles.radio_imgSimple handles.radio_imgShaded handles.radio_imgSS],'Value',0)

% ----------------------------------------------------------------------
function radio_imgSS_CB(hObject, handles)
	if (~get(hObject,'Val')),	set(hObject,'Val',1),	return,		end
	set([handles.radio_imgSimple handles.radio_imgShaded handles.radio_imgAmp],'Value',0)

% ----------------------------------------------------------------------
function edit_symbSize_CB(hObject, handles)
% Just check that impossible figures are not acepted
	str = get(hObject,'String');        s = str2double(str);
	if (isnan(s) || s <= 0),	set(hObject,'String','0.01'),	return,		end

% ----------------------------------------------------------------------
function push_help_CB(hObject, handles)


% ----------------------------------------------------------------------
function push_OK_CB(hObject, handles)
% ...
	if (get(handles.radio_pointCloud, 'Val'))
		contents = cellstr(get(handles.popup_symb,'String'));
		val = get(handles.popup_symb,'Value');
		switch contents{val}
			case 'Point',    symb = 7;
			case 'Cube',     symb = 3;
			case 'Sphere',   symb = 6;
			case 'Dyamond',  symb = 4;
			case 'Cylinder', symb = 5;
		end
		siz = str2double(get(handles.edit_symbSize,'String'));
		PTparams.Symbol = symb;	PTparams.PointRad = siz;
		par.TDRver = handles.TDRver;	par.proj = 'geog';		par.PTparams = PTparams;
		
		D = gmtmex(sprintf('mbgetdata -I%s', handles.fnameMB));
		xyz = [D(1).data(:) D(2).data(:) D(3).data(:)];
		bb = [min(xyz)' max(xyz)']';
		fname = [handles.path_tmp 'lixo.sd'];
		write_flederFiles('points', fname, xyz, 'first', bb(:), par);
		comm = [' -data ' fname ' &'];
		if (handles.whichFleder),	fcomm = ['iview4d' comm];		% Free viewer
		else						fcomm = ['fledermaus' comm];	% The real thing
		end
		try
			if (isunix)				% Stupid linux doesn't react to a non-existant iview4d
				resp = unix(fcomm);
				if (resp == 0)
					errordlg('I could not find Fledermaus. Hmmm, do you have it?','Error')
				end
			elseif (ispc)
				if (strcmp(handles.TDRver, '2.0')),	fcomm(6) = '3';		end		% Call iview3d
				delete(handles.figure1)		% It's useless now because Fleder will block it and
				if (handles.no_file && ishandle(handles.hMirFig)),	delete(handles.hMirFig),	end
				dos(fcomm);		% s is always 0 (success) as long as iview4d is accessible, even when the comm fails
			else
				errordlg('Unknown platform.','Error'),	return
			end
		catch
			errordlg('I could not find Fledermaus. Hmmm, do you have it?','Error')
		end

	else
		opt_Z = ' -Z5';
		if (get(handles.radio_imgSimple, 'Val')),		opt_Z = ' -Z2';
		elseif (get(handles.radio_imgShaded, 'Val')),	opt_Z = ' -Z3';
		elseif (get(handles.radio_imgAmp, 'Val')),		opt_Z = ' -Z4';
		end
		I = flipdim(gmtmex(['mbimport -I' handles.fnameMB opt_Z]),1);	% THE FLIP MUST BE DONE BY MBIMPORT
		move2side(handles.figure1,'right');		% Move it out of the way the best we can
		mirone(I)
	end
	
%-------------------------------------------------------------------------------------
function figure1_KeyPressFcn(hObject, eventdata)
	if isequal(get(hObject,'CurrentKey'),'escape')
		delete(hObject);
	end

% ----------------------------------------------------------------------
function show_MB_LayoutFcn(h1)

	set(h1, 'Position',[748 1173 179.5 266],...
		'Color',get(0,'factoryUicontrolBackgroundColor'),...
		'KeyPressFcn',@figure1_KeyPressFcn,...
		'MenuBar','none',...
		'Name','Show MB',...
		'NumberTitle','off',...
		'DoubleBuffer','on',...
		'Resize','off',...
		'HandleVisibility','Call',...
		'Tag','figure1');

	uicontrol('Parent',h1,'Position',[10 156 161 76],'Style','frame');

	uicontrol('Parent',h1, 'Position',[25 167 70 18],...
		'String',{'Point'; 'Cube'; 'Sphere'; 'Dyamond'; 'Cylinder'},...
		'Style','popupmenu',...
		'Value',1,...
		'BackgroundColor',[1 1 1],...
		'TooltipString','Symbol used in the point cloud',...
		'Tag','popup_symb');

	uicontrol('Parent',h1, 'Position',[101 167 49 19],...
		'String','0.01',...
		'Style','edit',...
		'BackgroundColor',[1 1 1],...
		'Call',@showMB_uiCB,...
		'TooltipString','Symbol size in unknown unites',...
		'Tag','edit_symbSize');

	uicontrol('Parent',h1, 'Position',[30 187 47 13],...
		'String','Symbol',...
		'Style','text',...
		'Tag','text_Symb');

	uicontrol('Parent',h1, 'Position',[109 187 30 13],...
		'String','Size',...
		'Style','text',...
		'Tag','text_size');

	uicontrol('Parent',h1,'Position',[10 31 161 122],'Style','frame');

	uicontrol('Parent',h1, 'Position',[18.5 203 75 17],...
		'String','Point cloud',...
		'Style','radiobutton',...
		'Value',1,...
		'Call',@showMB_uiCB,...
		'TooltipString','Show the point cloud.',...
		'Tag','radio_pointCloud');

	uicontrol('Parent',h1, 'Position',[18.5 122 75 17.5],...
		'String','Plot Image',...
		'Style','radiobutton',...
		'Call',@showMB_uiCB,...
		'TooltipString','Plot a quick view image',...
		'Tag','radio_image');

	uicontrol('Parent',h1, 'Position',[39 102 110 17],...
		'String','Simple bathymetry',...
		'Style','radiobutton',...
		'Call',@showMB_uiCB,...
		'Enable','off',...
		'TooltipString','Plain color image',...
		'Tag','radio_imgSimple');

	uicontrol('Parent',h1, 'Position',[39 80.5 110 17],...
		'String','Shaded bathymetry',...
		'Style','radiobutton',...
		'Call',@showMB_uiCB,...
		'Enable','off',...
		'TooltipString','Shaded illuminated color image',...
		'Tag','radio_imgShaded');

	uicontrol('Parent',h1, 'Position',[39 58.5 70 17],...
		'String','Amplitude',...
		'Style','radiobutton',...
		'Call',@showMB_uiCB,...
		'Enable','off',...
		'TooltipString','Gray scale amplitude plot',...
		'Tag','radio_imgAmp');

	uicontrol('Parent',h1, 'Position',[39 36 70 17],...
		'String','Side Scan',...
		'Style','radiobutton',...
		'Call',@showMB_uiCB,...
		'Enable','off',...
		'TooltipString','Gray scale Side Scan plot',...
		'Tag','radio_imgSS');

	uicontrol('Parent',h1, 'Position',[4.5 244 170 17],...
		'String','How to show this data?',...
		'Style','text',...
		'Tag','text2',...
		'FontSize',10,...
		'FontAngle','italic');

	uicontrol('Parent',h1, 'Position',[10 5 30 21],...
		'String','?',...
		'Call',@showMB_uiCB,...
		'TooltipString','Help/Info',...
		'Tag','push_help',...
		'FontSize',9,...
		'FontWeight','bold');

	uicontrol('Parent',h1, 'Position',[95.5 5 75.5 21],...
		'String','OK',...
		'Call',@showMB_uiCB,...
		'Tag','push_OK',...
		'FontSize',9,...
		'FontWeight','bold');

% ----------------------------------------------------------------------------------
function showMB_uiCB(hObject, evt)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));

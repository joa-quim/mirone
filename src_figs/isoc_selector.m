function varargout = isoc_selector(varargin)
% Helper win to select magnetic isochrons from private collections

%	Copyright (c) 2004-2018 by J. Luis
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

	if (isempty(varargin))
		errordlg('ISOC_SELECTOR: wrong number of input arguments.','Error'),	return
	end

	hObject = figure('Tag','figure1','Visible','off');
	isoc_selector_LayoutFcn(hObject);
	handles = guihandles(hObject);
	handlesMir = guidata(varargin{1});
	if (handlesMir.IamCompiled)
		move2side(hObject,'center');
		WindowAPI(hObject, 'TopMost')
	else
		move2side(handlesMir.figure1, hObject,'right');
	end

	% Load the directory list stored in mirone_pref
	s = load([handlesMir.path_data 'mirone_pref.mat'], 'isochrons_dir');
	try				% STUPID API returns a struct with no fields instead of an EMPTY s
		isocs_dir = s.isochrons_dir;
	catch
		isocs_dir = '';
	end
	set(handles.edit_isocsDir, 'Str', isocs_dir)

	handles.path_tmp  = handlesMir.path_tmp;
	handles.path_data = handlesMir.path_data;
	handles.hMirFig   = handlesMir.figure1;
	handles.last_dir  = handlesMir.last_dir;

	set(handles.listbox_isocs, 'Str', {'0' '2' '2a' '3' '3a' '4' '4a' '5' '5c' '6' '9' '13' '18' '20' '21' ...
		'22' '23' '24' '25' '26' '28' '29' '30' '32' '33' '33r' 'M0' 'M1' 'M5' 'M10' 'M16' 'M22' 'M25'})
	plates = {'' 'Africa' 'Acores' 'Eurasia' 'Greenland' 'Iberia' 'Jan Mayen' 'Kings Trough' 'Kings South' 'Lomonosov' ...
		'Mohns' 'North America' 'XX1' 'XX2' 'XXX'};
	set(handles.popup_plate1, 'Str', plates, 'Val', 12)
	set(handles.popup_plate2, 'Str', plates, 'Val', 1)
	guidata(hObject, handles);
	set(hObject,'Visible','on');
	if (nargout),	varargout{1} = hObject;		end

% --------------------------------------------------------------------------------------
function radio_oneOnly_CB(hObject, handles)
	if (~get(hObject,'Val')),	set(hObject,'Val',1),	return,		end
	set(handles.radio_both, 'Val',0)

% --------------------------------------------------------------------------------------
function radio_both_CB(hObject, handles)
	if (~get(hObject,'Val')),	set(hObject,'Val',1),	return,		end
	set(handles.radio_oneOnly, 'Val',0)

% --------------------------------------------------------------------------------------
function push_isocsDir_CB(hObject, handles)
% Selelect the directory of isochrons and save current choice ... unless it's faulty
	if (strcmp(computer, 'PCWIN'))
		isochrons_dir = uigetfolder_win32('Select a directory', cd);
		if (~isochrons_dir),	isochrons_dir = '';		end		% To make it compatible with the other branch
	else            % This guy doesn't let to be compiled
		isochrons_dir = uigetdir(cd, 'Select a directory');
	end
	if (isempty(isochrons_dir)),	return,		end
	isochrons_dir(end+1) = '/';

	fname = [handles.path_data 'mirone_pref.mat'];
	version7 = version;			% Find matlab version in use. Only interested to know if R13 or >= R14
	V7 = (sscanf(version7(1),'%f') > 6);
	if (~V7)					% R <= 13
		save(fname,'isochrons_dir', '-append')
	else
		save(fname,'isochrons_dir', '-append', '-v6')
	end

	edit_isocsDir_CB(handles.edit_isocsDir, handles, isochrons_dir)

% --------------------------------------------------------------------------------------
function edit_isocsDir_CB(hObject, handles, opt)
% Set path of dir where isochrons files lieve
	if (nargin == 3)					% Set by push_pato_CB()
		set(hObject, 'Str', opt)
	else
		str = get(hObject, 'Str');
		if (exist(str, 'dir') ~= 2)	% If entered manually, check calidity
			errordlg('This directory does not exist. Please pay attention and provide a good one.', 'Error')
			set(hObject, 'Str', '')
		end
	end

% --------------------------------------------------------------------------------------
function push_OK_CB(hObject, handles)
% Check if all is right and do the work.

	isocsDir = get(handles.edit_isocsDir, 'Str');
	if (isempty(isocsDir))
		errordlg('Hey, isochrons directory info is empty. How can I know where your isochron files are? Please tell me that via the edit box above.', 'Error')
		return
	end
	if (isocsDir(end) ~= '/' || isocsDir(end) ~= '\'),		isocsDir(end+1) = '/';		end

	strIsocs = get(handles.listbox_isocs, 'Str');
	selection = get(handles.listbox_isocs, 'Val');

	str = get(handles.popup_plate1, 'Str');
	plate_1 = str{get(handles.popup_plate1, 'Val')};
	if (strfind(plate_1, ' '))
		[t, r] = strtok(plate_1);
		r = ddewhite(r);
		plate_1 = ['_' upper(t(1)) upper(r(1))];
	else
		plate_1 = ['_' upper(plate_1(1:2))];
	end

	% Get second Plate. If empty pick all plates that share a border with Plate1
	str = get(handles.popup_plate2, 'Str');
	plate_2 = str{get(handles.popup_plate2, 'Val')};
	if (~isempty(plate_2))
		if (strfind(plate_2, ' '))
			[t, r] = strtok(plate_2);
			r = ddewhite(r);
			plate_2 = ['_' upper(t(1)) upper(r(1))];
		else
			plate_2 = ['_' upper(plate_2(1:2))];
		end
	end

	% Load the selected isocs
	msg = cell(numel(selection), 1);	c = false(numel(selection), 1);
	for (k = 1:numel(selection))
		msg{k} = load_this_isoc(handles, isocsDir, strIsocs{selection(k)}, plate_1, plate_2);
		if (isempty(msg{k})),	c(k) = true;	end
	end
	msg(c) = [];
	if (~isempty(msg)),		warndlg(msg, 'Warning'),	end

	% Reset the last_dir that was meanwhile set to tmp by the load_xyz() call
	handMir = guidata(handles.hMirFig);
	handMir.last_dir = handles.last_dir;
	guidata(handles.hMirFig, handMir)

% -----------------------------------------------------------------------------------
function msg = load_this_isoc(handles, isocsDir, isoc, plate_1, plate_2)
% Load Isochron ISOC or return an warning message if not found in ISOCDIR
	msg = '';
	all = dir([isocsDir 'c' isoc '_*.dat']);
	c = false(numel(all), 1);
	if (isempty(plate_2))				% Here we want all combinations that contain the Plate_1
		for (k = 1:numel(all))
			if (~isempty(strfind(all(k).name, plate_1)))
				c(k) = true;
			end
		end
	else								% Here we want the combination Plate1_Plate2 and/or Plate2_Plate1
		AB = [plate_1 plate_2];		BA = [plate_2 plate_1];		% Select only these two pairs
		if (get(handles.radio_both, 'Val'))		% Either AB or BA is good.
			for (k = 1:numel(all))
				if (~isempty(strfind(all(k).name, AB)) || ~isempty(strfind(all(k).name, BA)))
					c(k) = true;
				end
			end
		else							% Retain only the AB combination
			for (k = 1:numel(all))
				if (~isempty(strfind(all(k).name, AB)))
					c(k) = true;
				end
			end
		end
	end

	all = all(c);			% Retain only those that passed the above tests
	if (isempty(all))
		msg = ['Found no isochron ' isoc ' in the specified path.'];
		return
	end

	list = cell(numel(all), 1);
	for (k = 1:numel(all))
		list{k} = [isocsDir all(k).name];
	end
	fid = fopen([handles.path_tmp 'isoc_pairs.txt'], 'wt');
	if (fid < 0)
		errordlg(['Error writing file ' handles.path_tmp 'isoc_pairs.txt'], 'Error')
		return
	end
	fprintf(fid, '>HAVE_INCLUDES +proj=longlat\n');
	for (k = 1:numel(list))
		fprintf(fid, '> INCLUDE=%s\n', list{k});
	end
	fclose(fid);

	% Now load the selected Isochrons
	load_xyz(guidata(handles.hMirFig), [handles.path_tmp 'isoc_pairs.txt'], 'Isochron')

% -----------------------------------------------------------------------------------
function figure1_KeyPressFcn(hObject, eventdata)
	handles = guidata(hObject);
	if isequal(get(hObject,'CurrentKey'),'escape')
		delete(handles.figure1);
	end

% --------------------------------------------------------------------------------------
function h1 = isoc_selector_LayoutFcn(h1)

set(h1, 'Position',[520 608 240 220],...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'KeyPressFcn',@figure1_KeyPressFcn,...
'MenuBar','none',...
'Name','Load Isocs',...
'NumberTitle','off',...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1');

uicontrol('Parent',h1, 'Position',[10 196 68 14], 'FontSize',9,...
'String','Isochrons',...
'Style','text');

uicontrol('Parent',h1, 'Position',[10 10 69 181],...
'BackgroundColor',[1 1 1],...
'String',{'2'},...
'Style','listbox',...
'max',2,...
'Value',1,...
'Tag','listbox_isocs');

uicontrol('Parent',h1, 'Position',[115 196 68 14], 'FontSize',9,...
'String','Plate A',...
'Style','text');

uicontrol('Parent',h1, 'Position',[90 171 141 20],...
'BackgroundColor',[1 1 1],...
'String','',...
'Style','popupmenu',...
'TooltipString','Selct plate',...
'Value',1,...
'Tag','popup_plate1');

uicontrol('Parent',h1, 'Position',[115 145 68 14], 'FontSize',9,...
'String','Plate B',...
'Style','text');

uicontrol('Parent',h1, 'Position',[90 121 141 20],...
'BackgroundColor',[1 1 1],...
'Style','popupmenu',...
'TooltipString','Leave empty to select "Plate1 against all its neighbors"',...
'Value',1,...
'Tag','popup_plate2');

uicontrol('Parent',h1, 'Pos',[100 84 44 20],...
'Callback',@isoc_selector_uiCB,...
'String','Both',...
'Style','radiobutton',...
'Value', 1, ...
'Tooltip','Pick both combinations (AB and BA)',...
'Tag','radio_both')

uicontrol('Parent',h1, 'Pos',[151 84 65 20],...
'Callback',@isoc_selector_uiCB,...
'String','Only one',...
'Style','radiobutton',...
'Tooltip','Pick ony the selected combination (AB or BA)',...
'Tag','radio_oneOnly')

uicontrol('Parent',h1, 'Position',[90 40 121 21],...
'Style','edit',...
'Callback',@isoc_selector_uiCB,...
'TooltipString','Path of the isochrons files',...
'Tag','edit_isocsDir');

uicontrol('Parent',h1, 'Position',[211 40 21 21],...
'String','...',...
'Callback',@isoc_selector_uiCB,...
'FontSize',9,...
'FontWeight','bold', ...
'Tag','push_isocsDir');

uicontrol('Parent',h1, 'Position',[151 10 80 22], 'FontSize',9, 'FontWeight','bold',...
'Callback',@isoc_selector_uiCB,...
'String','OK',...
'Tag','push_OK');

function isoc_selector_uiCB(hObject, eventdata)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));
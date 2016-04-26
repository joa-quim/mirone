function varargout = update_gmt(varargin)
% Helper fig to update a GMT5 installation by generate batch to download only the needed files

%	Copyright (c) 2004-2016 by J. Luis
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

% $Id: update_gmt.m 7879 2016-04-26 16:12:02Z j $

	hObject = figure('Tag','figure1','Visible','off');
	update_gmt_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hObject,'right')

	if (nargin)			% Than varargin must hold the Mirone handles
		handles.path_tmp = varargin{1}.path_tmp;
		handles.home_dir = varargin{1}.home_dir;
		handles.IamCompiled = varargin{1}.IamCompiled;
	else
		mir_dirs = getappdata(0,'MIRONE_DIRS');
		if (~isempty(mir_dirs))
			handles.path_tmp = [mir_dirs.home_dir '/tmp/'];
			handles.home_dir = mir_dirs.home_dir;
		else
			handles.path_tmp = [pwd '/'];
			handles.home_dir = pwd;
		end
		handles.IamCompiled = false;		% Though we actually do not know
	end

	% Try to find which GMT5 is currently in use and find if it's a 64 or 32 bits version.
	try
		[s, w] = mat_lyies('gmt --show-bindir');
	catch			% Falls here if no GMT5 around
		s = 1;
	end
	if (s == 0)
		try
			patoGMT = w(1:end-4);		% Last char is the \n
			[s, w2] = mat_lyies(['pesnoop ' patoGMT '/bin/psxy.exe /PE_CD']);
			if (s == 0)
				ind = strfind(w2, 'bit Portable');
				if (~isempty(ind))
					if (w2(ind(1)-1) == '4')
						set(handles.edit_path64, 'Str', patoGMT)
					else
						set(handles.edit_path32, 'Str', patoGMT)
					end
				end
			end
		catch
			disp(lasterr)
		end
	end

	guidata(hObject, handles);
	set(hObject,'Visible','on');
	if (nargout),	varargout{1} = hObject;		end

% ---------------------------------------------------------------------
function push_path64_CB(hObject, handles)
	gmt_dir = uigetfolder_win32('Select the GMT5 root directory', pwd);
	if (isempty(gmt_dir)),	return,		end
	edit_path64_CB(handles.edit_path64, handles, gmt_dir)

% ---------------------------------------------------------------------
function edit_path64_CB(hObject, handles, opt)
	if (nargin == 3),	pato = opt;
	else				pato = get(hObject, 'Str');
	end
	if (~isempty(pato))
		pato = ddewhite(pato);	set(hObject, 'Str', pato)		% Make sure no blanks to f us up.
		sanitize_gmt_path(handles.edit_path32, pato)
	end

% ---------------------------------------------------------------------
function push_path32_CB(hObject, handles)
	gmt_dir = uigetfolder_win32('Select the GMT5 root directory', pwd);
	if (isempty(gmt_dir)),	return,		end
	edit_path64_CB(handles.edit_path32, handles, gmt_dir)

% ---------------------------------------------------------------------
function edit_path32_CB(hObject, handles, opt)
	if (nargin == 3),	pato = opt;
	else				pato = get(hObject, 'Str');
	end
	if (~isempty(pato))
		pato = ddewhite(pato);	set(hObject, 'Str', pato)		% Make sure no blanks to f us up.
		sanitize_gmt_path(handles.edit_path64, pato)
	end

% ---------------------------------------------------------------------
function sanitize_gmt_path(hObject, pato)
% Common code used by both edit_path??_CB functions to do some sanity checks on the path provided
	if (exist(pato, 'dir') ~= 7)
		errordlg('This directory does not exist.','Error')
		set(hObject, 'Str', ''),		return
	else
		if (~isempty(get(hObject, 'Str')))		% Attention, this 'hObject' is that of the complementary bitage edit box
			warndlg('Cannot update more than one version at same time. Make up your mind 64 or 32 bits?', 'Warn')
			set(hObject, 'Str', ''),	return
		end
	end
	if (exist([pato '/bin/psxy.exe'], 'file') ~= 2)
		errordlg('This is not a root directory of a GMT installation.','Error')
		set(hObject, 'Str', '')
	end

% ---------------------------------------------------------------------
function radio_DL_PT_CB(hObject, handles)
	set(hObject,'Val', 1)		% Only 'option' by now

% ---------------------------------------------------------------------
function radio_DL_US_CB(hObject, handles)
	warndlg('Sorry, not yet available', 'Warning')
	set(hObject,'Val', 0)		% Not yet possible

% ---------------------------------------------------------------------
function push_OK_CB(hObject, handles)
% ...
	do_64 = true;
	patoGMT = get(handles.edit_path64, 'Str');
	if (isempty(patoGMT))
		patoGMT = get(handles.edit_path32, 'Str');
		do_64 = false;
	end
	if (isempty(patoGMT))
		errordlg('Yes, update what?', 'Chico Clever'),	return
	end
	patoGMT = [ddewhite(patoGMT) '/'];

	dest_fiche = [handles.path_tmp 'apudeitaGMT.txt'];		url = 'w3.ualg.pt/~jluis/GMTshiny_w64/';
	if (~do_64),	url(end-2:end-1) = '32';	end			% Redirect to 32 bits dir
	dos(['wget "' url 'apudeitaGMT.txt' '" -q --tries=2 --connect-timeout=5 -O ' dest_fiche]);
	finfo = dir(dest_fiche);
	if (finfo.bytes == 0)
		builtin('delete',dest_fiche);
		msgbox('Nothing to update right now.','Nothing New1'),	return
	end
	fid = fopen(dest_fiche,'rt');
	todos = fread(fid,'*char');		fclose(fid);
	[MD5, nomes] = strread(todos,'%s %s');
	builtin('delete',dest_fiche);			% Remove this one right away
	namedl = cell(numel(nomes),1);			% To hold the names of the to-be-updated files

	c = true(numel(nomes),1);
	for (k = 1:numel(nomes))
		nome = [patoGMT nomes{k}];
		if (exist(nome, 'file') ~= 2)		% If it does not exist
			namedl{k} = [url nomes{k}];		% If file does not exist locally, it's for sure one to update
			c(k) = false;
		else
			localMD5 = CalcMD5(nome, 'file');
			if (~strcmpi(MD5{k}, localMD5))	% OK, we have a new version of this guy
				namedl{k} = [url nomes{k}];	% File name to update
				c(k) = false;
			end		
		end
	end
	namedl(c) = [];		nomes(c) = [];	MD5(c) = [];	% Remove names that don't need to be updated
	if (isempty(namedl))
		msgbox('Your GMT is updated to the latest.','Nothing New2'),	return
	end

	fid1 = fopen([patoGMT 'dowload.bat'],'wt');	% Create the downloading batch
	fprintf(fid1, '@echo off\nREM Batch file to download the latest GMT files needed to update your instalation\n\n');
	fprintf(fid1, 'REM Create a tmp directory to receive the downloaded files\n');
	t = strrep([patoGMT 'tmp'], '/', '\');		% Otherwise, stupid Windows says "The syntax of the command is incorrect"
	fprintf(fid1, 'IF NOT EXIST %s md %s\n\necho downloading files ...\n', t, t);
	fid2 = fopen([patoGMT 'update.bat'],'wt');	% Create the updating batch
	fprintf(fid2, '@echo off\nREM Batch file to move the dowloaded files to their finally destiny\n\n');
	fprintf(fid2, 'IF NOT EXIST share\\postscriptlight md share\\postscriptlight\n');
	fprintf(fid2, 'IF NOT EXIST share\\spotter md share\\spotter\n');
	fprintf(fid2, 'cd tmp\n');
	prefix = '';
	if (handles.IamCompiled),		prefix = [handles.home_dir '/'];	end
	for (k = 1:numel(namedl))
		[pato, nome, ext] = fileparts(namedl{k});
		if (~strcmp(nome, 'VERSION'))
			fprintf(fid1, '%swget %s -q --tries=2 --connect-timeout=5 -O %stmp/%s\n', prefix, namedl{k}, patoGMT, [nome ext]);
			fprintf(fid2, 'move /Y %s\t../%s\n', [nome ext], nomes{k});
		else
			fprintf(fid2, 'echo %s > ../share/VERSION\n', MD5{k});
		end
	end
	fprintf(fid1, 'echo Finished dowload. Now manually run the updtate batch to finish the update\npause\n');
	fprintf(fid2, 'echo Finished. Your GMT5 should be updated now\npause\n');
	fclose(fid1);	fclose(fid2);
	
	delete(handles.figure1)
	msgbox(sprintf('OK, now you must ran manually these two batch files:\n\n%s\n\n%s', ...
		[patoGMT 'dowload.bat'], [patoGMT 'update.bat']), 'Update Info')

% ------------------------------------
function h1 = update_gmt_LayoutFcn(h1)

set(h1, 'Position',[520 589 500 211],...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'MenuBar','none',...
'Name','Update GMT5',...
'NumberTitle','off',...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1');

uicontrol('Parent',h1, 'Position',[10 103 91 16],...
'FontSize',9,...
'HorizontalAlignment','left',...
'String','GMT5 64 path',...
'Style','text',...
'Tag','text1');

uicontrol('Parent',h1, 'Position',[100 100 371 22],...
'BackgroundColor',[1 1 1],...
'Call',@updateGMT_uiCB,...
'HorizontalAlignment','left',...
'String',blanks(0),...
'Style','edit',...
'TooltipString','Path to the GMT5 64 bits root directory.',...
'Tag','edit_path64');

uicontrol('Parent',h1, 'Position',[470 99 23 24],...
'Call',@updateGMT_uiCB,...
'FontWeight','bold',...
'String','...',...
'TooltipString','Browse for the GMT5 64 bits path',...
'Tag','push_path64');

uicontrol('Parent',h1, 'Position',[10 63 91 16],...
'FontSize',9,...
'HorizontalAlignment','left',...
'String','GMT5 32 path',...
'Style','text',...
'Tag','text2');

uicontrol('Parent',h1, 'Position',[100 60 371 22],...
'BackgroundColor',[1 1 1],...
'Call',@updateGMT_uiCB,...
'HorizontalAlignment','left',...
'Style','edit',...
'TooltipString','Path to the GMT5 32 bits root directory.',...
'Tag','edit_path32');

uicontrol('Parent',h1, 'Position',[470 59 23 24],...
'Call',@updateGMT_uiCB,...
'FontWeight','bold',...
'String','...',...
'TooltipString','Browse for the GMT5 32 bits path',...
'Tag','push_path32');

uicontrol('Parent',h1, 'Position',[50 140 380 67],...
'FontSize',10,...
'HorizontalAlignment','left',...
'String',{'Update your GMT5 installation to the latest developing version.'; 'No documentation or coastlines are updated here, only binaries.'; 'At the end, this program generates two batch files that you will'; 'need to run manually to finish the process.' },...
'Style','text',...
'Tag','text3');

uicontrol('Parent',h1, 'Position',[10 18 120 23],...
'Call',@updateGMT_uiCB,...
'String','Download from PT',...
'Style','radiobutton',...
'TooltipString','Download from the University of Algarve, Portuguese server',...
'Value',1,...
'Tag','radio_DL_PT');

uicontrol('Parent',h1, 'Position',[150 18 141 23],...
'Call',@updateGMT_uiCB,...
'String','Download from Hawaii',...
'Style','radiobutton',...
'TooltipString','Download from the Hawaii GMT server',...
'Tag','radio_DL_US');

uicontrol('Parent',h1, 'Position',[411 9 81 26],...
'Call',@updateGMT_uiCB,...
'FontSize',10,...
'FontWeight','bold',...
'String','OK',...
'Tag','push_OK');

function updateGMT_uiCB(hObject, evt)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));
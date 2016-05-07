function varargout = snapshot(varargin)
% GUI interface to the imcapture function
%
%   snapshot(H) operates on the image contents of the UNIQUE image in Figure whose handle is H
%   snapshot(H,'whatever') as above but captures image and frame.
%   snapshot(H,'frame') as above but captures image and frame.
%   snapshot(H,'noname') Form used for saving Georefed images where format is not selectable here
%   snapshot(H,'img') same as snapshot(H)

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

% $Id$

	if (isempty(varargin))
	    errordlg('SNAPSHOT: wrong number of input arguments.','Error'),		return
	end

	hObject = figure('Vis','off');
	snapshot_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hObject,'center')

	handles.hCallingFig = varargin{1};
	handlesMir = guidata(handles.hCallingFig);
	if (handlesMir.no_file)
		errordlg('You didn''t even load a file. What are you expecting then?','Error')
		delete(hObject),	return
	end

	handles.imgOnly = 1;            % Flag that indicates pure image capture
	handles.noname = false;
	if (numel(varargin) > 1 && strcmp(varargin{2},'frame'))
		handles.imgOnly = 0;        % OK, so we will have a image + frame capture
	elseif (numel(varargin) > 1 && strcmp(varargin{2},'noname'))
		handles.noname = true;
	elseif (numel(varargin) > 1 && strcmp(varargin{2},'img'))
		% Already the default
	end

	% Add this figure handle to the carraças list
	plugedWin = getappdata(handles.hCallingFig,'dependentFigs');
	plugedWin = [plugedWin hObject];
	setappdata(handles.hCallingFig,'dependentFigs',plugedWin);

	% -------------- Find if we have any graphical objects ploted
	ls = findobj(handles.hCallingFig,'Type','line');
	ps = findobj(handles.hCallingFig,'Type','patch');
	ts = findobj(handles.hCallingFig,'Type','text');
	fs = findobj(handles.hCallingFig, '-depth',1, 'Style','frame', 'Tag', 'FancyFrame');
	handles.imgIsClean = false;
	if (isempty(ls) && isempty(ps) && isempty(ts) && isempty(fs))
		handles.imgIsClean = true;		% We won't need to do screen capture if resolution is one-to-one
	end

	% -------------- Fill the format popup list
	str1 = {'JPEG image (*.jpg)'; ...
		'Portable Network Graphics (*.png)'; ...
		'Tagged Image File (*.tif)'; ...
		'Windows Bitmap (*.bmp)'; ...
		'EPS files (*.eps)'; ...
		'Adobe Illustrator (*.ai)'; ...
		'Enhanced metafiles (*.emf)'; ...
		'Hieralchical Data Format (*.hdf)'; ...
		'GIF image (*.gif)'; ...
		'Windows Paintbrush (*.pcx)'; ...
		'Portable Anymap (*.pnm)'; ...
		'SUN rasterfile (*.ras)'; ...
		'Raw RGB format (*.raw)'};

	handles.exts = {'.jpg' '.png' '.tif' '.bmp' '.eps' '.ai' '.emf' '.hdf' '.gif' '.pcx' '.pnm' '.ras' '.raw'};

	% ------------ Get the image dimensions. Note that the dimensions are those deffined
	%              by the axes size and not the true image size.
	axUnit = get(handlesMir.axes1,'Units');         set(handlesMir.axes1,'Units','pixels')
	axPos = get(handlesMir.axes1,'pos');            set(handlesMir.axes1,'Units',axUnit)
	handles.imSize = size(get(handlesMir.hImg,'CData'));
	handles.imAxSize = [axPos(4) axPos(3)];                     % Image size as reinterpolated to fit inside axes.

	if (~handles.noname)
		if (~handles.imgIsClean || numel(handles.imSize) == 3)  % Output will be RGB for sure. So remove formats
			ind = false(numel(str1),1);
			ind(9:11) = true;               % Remove gif, pcx & pnm
			str1(ind) = [];         handles.exts(ind) = [];
		end
		set(handles.popup_fileType,'String',str1)
	else
		handles.imAxSize = handles.imSize(1:2);                 % Geo-referenced files. Make it multiples of original size.
	end

	% ---------------- Get info to allow guessing imgAx image size
	PU = get(handles.hCallingFig,'paperunits');     set(handles.hCallingFig,'paperunits','inch')
	pp = get(handles.hCallingFig,'paperposition');  set(handles.hCallingFig,'paperunits',PU)
	dpi = round(handles.imSize(2) / pp(3));
	handles.pp = [pp(3) handles.imSize(1) / dpi];       % PaperPosition in inches as it will be set inside imcapture

	% ---------------- Fill the edit image size with the default values
	if (handles.imgOnly)
		nRows = round(handles.imAxSize(1));     nCols = round(handles.imAxSize(2));
		origMegs = prod(handles.imSize) / 1048576;   % Original image size in Mb (use prod so it works for indexed too)
		sizeOrigUnits = ' Megs';
		if (origMegs < 1)
			sizeOrigUnits = ' Kb';  origMegs = origMegs * 1024;
		end
		handles.txtOrigSize = sprintf('(%dx%d) %.1f%s',handles.imSize(1),handles.imSize(2),origMegs,sizeOrigUnits);
	else
		nRows = round(handles.pp(1) * 150);     nCols = round(handles.pp(2) * 150);
	end
	Megs = nRows * nCols * 3 / 1048576;
	sizeUnits = ' Megs';
	if (Megs < 1)               % Report file size in Kbytes
		sizeUnits = ' Kb';      Megs = Megs * 1024;
	end
	handles.txtThisSize = sprintf('(%dx%d) %.1f%s',nRows,nCols,Megs,sizeUnits);

	% ---------------- Fill the edit file name with a default value ----------------------------
	if (~handles.noname)
		fname = get(handlesMir.figure1,'Name');
		ind = strfind(fname, ' @ ');
		if (fname(end) == '%' && ~isempty(ind))		% Get rid of eventual ' @ xxx%' in the name
			fname(ind(1):end) = [];
			while (fname(end) == ' '),	fname(end) = [];	end		% ugly but may execute only once or twice.
		end
		fname = strrep(fname,' ','_');
		[pato, fname] = fileparts(fname);
		if (islogical(get(handlesMir.hImg,'CData')) && handles.imgOnly)   % Best proposition when we have a logical (mask) image
			fname = [handlesMir.last_dir fname '.png'];
			set(handles.popup_fileType,'Val',2)
			set(handles.checkbox_origSize,'Val',1)
			set(handles.slider_quality,'Visible','off')
			set(handles.text_Quality,'Visible','off')
			set(handles.text_qualityLev,'Visible','off')
			handles.txtThisSize = handles.txtOrigSize;
		else
			prefix = handlesMir.last_dir;
			if (prefix(end) ~= '\' && prefix(end) ~= '/')		% For a bloody reason compiled version looses '\' at the end
				prefix = [prefix '/'];		% '/' works for Windows as well
			end
			fname = [prefix fname '.jpg'];
		end
		set(handles.edit_fname,'String',fname);
	else
		set(handles.edit_fname,'String','Not used','Enable','off');
		set(handles.popup_fileType,'String','Nikles','Enable','off')
		set(handles.push_outFile,'Enable','off')
		set(handles.edit_imgSize,'String',handles.txtOrigSize)
		set(handles.slider_quality,'Visible','off')
		set(handles.text_Quality,'Visible','off')
		set(handles.text_qualityLev,'Visible','off')
		set(handles.checkbox_origSize,'Val',1,'Enable','inactive')
	end
	set(handles.edit_imgSize,'String',handles.txtThisSize);

	% ----------------- Set the slider with the apropriate range
	if (handles.imgOnly)
		sliderRange(handles,1,1)		% Sets in Magnification mode
	else
		set(handles.checkbox_origSize,'Enable','off')
		sliderRange(handles,0,150)		% Sets in DPI mode
	end
	if (handles.noname)
		set(handles.slider_mag,'Val',1,'Enable','inactive')	% The above call to sliderRange insists in make it active
	end

	handles.hImg = handlesMir.hImg;
	handles.hCallingAx = handlesMir.axes1;
	handles.currMag = 1;			% Current magnification
	handles.quality = 75;			% Current quality level. Only applyes to the jpeg format
	handles.currDPI = 150;			% Current DPI for rasters
	handles.currVecDPI = 300;		% Current DPI for vector graphics
	handles.vecGraph = 0;			% To signal when we are dealing with vector graphics
	handles.handlesMir = handlesMir;

	% ----------------- Choose default command line output for snapshot
	handles.output = [];
	guidata(hObject, handles);
	set(hObject,'Visible','on')

	if (nargout)
		uiwait(handles.figure1);
		handles = guidata(handles.figure1);
		varargout{1} = handles.output;
		delete(handles.figure1)
	end

% -----------------------------------------------------------------------------
function sliderRange(handles,magnification,val)
% Set the slider range and values for the magnification or DPI cases
	if (magnification)
		set(handles.slider_mag,'Min',1,'Max',20,'Val',val,'SliderStep',[1 2]/19)
		set(handles.edit_mag,'String',val,'Tooltipstring','')
		set(handles.text_mag,'String','Magnification')
	else
		set(handles.slider_mag,'Min',50,'Max',1000,'Val',val,'SliderStep',[50 100]/950)
		set(handles.edit_mag,'String',val,'Tooltipstring','Do capture at this DPI resolution')
		set(handles.text_mag,'String','Resolution')
	end

% -----------------------------------------------------------------------------
function edit_fname_CB(hObject, handles)
	fname = get(hObject,'String');
	[pato,fname,ext] = fileparts(fname);
	if (~any(strcmpi(ext, handles.exts)))
		errordlg('You cannot choose a different file format than the ones offered to you.','ERROR')
		set(hObject,'String','')
	else
		handles.handlesMir.last_dir = [pato];
		guidata(handles.handlesMir.figure1, handles.handlesMir)		% Update the Mirone handles
	end

% -----------------------------------------------------------------------------
function push_outFile_CB(hObject, handles)
	contents = get(handles.popup_fileType,'String');
	val = get(handles.popup_fileType,'Val');    
	str = {handles.exts{val} contents{val}};
	handMir = guidata(handles.hCallingFig);
	[FileName,PathName] = put_or_get_file(handMir, str,'Select File name','put', handles.exts{val}); % FORCE ext (compiler bug)
	if isequal(FileName,0),		return,		end

	[pato,fname,ext] = fileparts(FileName);
	if (~isempty(ext) && ~strcmpi(ext,handles.exts{val}))
		errordlg(['You cannot choose a different file format here. It has to be of type ' handles.exts{val}],'ERROR')
		return
	end
	set(handles.edit_fname,'String',[PathName FileName])

% -----------------------------------------------------------------------------
function popup_fileType_CB(hObject, handles)
	ext = handles.exts{get(hObject,'Value')};
	fname = get(handles.edit_fname,'String');
	fname = stripExt(fname);
	set(handles.edit_fname,'String',[fname ext]);
	if (strcmp(ext,'.jpg'))
		set(handles.slider_quality,'Visible','on');         set(handles.text_Quality,'Visible','on')
		set(handles.text_qualityLev,'Visible','on');        set(handles.edit_imgSize,'String',handles.txtThisSize);
		val = handles.currDPI;      magn = 0;               % Depending on raster capture mode the current dpi is different
		if (handles.imgOnly)
			val = handles.currMag;      magn = 1;
			set(handles.checkbox_origSize,'Enable','on')
		end
		sliderRange(handles,magn,val);						handles.vecGraph = 0;
	elseif ( any(strcmp(ext,{'.ps' '.eps' '.ai' '.emf'})) )
		set(handles.slider_quality,'Visible','off');        set(handles.text_Quality,'Visible','off')
		set(handles.text_qualityLev,'Visible','off');       set(handles.edit_imgSize,'String','Don''t know')
		set(handles.checkbox_origSize,'Enable','off')
		sliderRange(handles,0,handles.currVecDPI);          handles.vecGraph = 1;
	elseif ( ~get(handles.checkbox_origSize,'Val') && strcmp(ext,'.gif'))
		warndlg('Gif files do not support RGB. Check the "Preserve Image ..." checkbox before using this option.','Warning')
		set(handles.edit_fname,'String',[fname '.png']);	% Replace the censured gif by a png
		set(hObject,'Val',2)
		set(handles.slider_quality,'Visible','off');        set(handles.text_Quality,'Visible','off')
		set(handles.text_qualityLev,'Visible','off');       set(handles.edit_imgSize,'String',handles.txtThisSize);
		val = handles.currDPI;      magn = 0;				% Depending on raster capture mode the current dpi is different
		if (handles.imgOnly)
			val = handles.currMag;      magn = 1;
			set(handles.checkbox_origSize,'Enable','on')
		end
		sliderRange(handles,magn,val);						handles.vecGraph = 0;
	else
		set(handles.slider_quality,'Visible','off');		set(handles.text_Quality,'Visible','off')
		set(handles.text_qualityLev,'Visible','off');		set(handles.edit_imgSize,'String',handles.txtThisSize);
		val = handles.currDPI;      magn = 0;				% Depending on raster capture mode the current dpi is different
		if (handles.imgOnly)
			val = handles.currMag;      magn = 1;
			set(handles.checkbox_origSize,'Enable','on')
		end
		sliderRange(handles,magn,val);						handles.vecGraph = 0;
	end
	guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------
function slider_mag_CB(hObject, handles)
% Callback to control the image magnetization factor slider
	mag = get(hObject,'Value');
	set(handles.edit_mag,'String',mag);
	if (~handles.vecGraph)
		if (handles.imgOnly)		% Image size as reinterpolated to fit inside axes
			nRows = handles.imAxSize(1) * mag;      nCols = handles.imAxSize(2) * mag;
		else						% handles.pp = PaperPosition
			nRows = round(handles.pp(1) * mag);     nCols = round(handles.pp(2) * mag);
		end
		Megs = nRows * nCols * 3 / 1048576;
		sizeUnits = ' Megs';
		if (Megs < 1)               % Report file size in Kbytes
			sizeUnits = ' Kb';      Megs = Megs * 1024;
		end
		handles.txtThisSize = [sprintf('(%d',nRows) 'x' sprintf('%d) ',nCols) sprintf('%.1f',Megs) sizeUnits];
		set(handles.edit_imgSize,'String',handles.txtThisSize);
	end
	if (handles.vecGraph),          handles.currVecDPI = mag;
	elseif (handles.imgOnly),       handles.currMag = mag;
	else                            handles.currDPI = mag;
	end
	guidata(handles.figure1,handles)
    
% -----------------------------------------------------------------------------
function checkbox_origSize_CB(hObject, handles)
% When checked, image size reflects the original size, which cannot be changed.
	if (get(hObject,'Val'))
		set(handles.edit_imgSize,'String',handles.txtOrigSize)
		set(handles.edit_mag,'String',1)
		set(handles.slider_mag,'Val',1,'Enable','inactive')
	else
		set(handles.edit_imgSize,'String',handles.txtThisSize)
		set(handles.edit_mag,'String',handles.currMag)
		set(handles.slider_mag,'Val',handles.currMag,'Enable','on')
	end

% -----------------------------------------------------------------------------
function slider_quality_CB(hObject, handles)
	handles.quality = round(get(hObject,'Value') * 100);
	set(handles.text_qualityLev,'String',[sprintf('%d',handles.quality) ' %']);
	guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------
function push_save_CB(hObject, handles)
	try
		im = get(handles.hImg,'CData');
	catch
		errordlg('Something stupid occured. Probably on your side. Did you kill the Figure?','Error')
		return
	end

	fname = get(handles.edit_fname,'String');
	[pato, name, EXT] = fileparts(fname);
	if (isempty(name))
		errordlg('Saving in a no name! You think you are funny?','Error');			return
	end

	dims = handles.imAxSize * handles.currMag;
	if (get(handles.checkbox_origSize,'Val')),      dims = handles.imSize(1:2);		end
	handsStBar = getappdata(handles.hCallingFig,'CoordsStBar');
	set(handsStBar,'Visible','off')
	set(handles.figure1,'Visible','off');       pause(0.01)     % To remove an cripled visual effect
	% The next flipdim() cases are a result of MANY hours trying to figure out where the guys in TMW
	% f... had the head when they invented that MESS of image fliped up/down depending on axes 'YDir'
	% is normal or STUPIDLY Y->down and, if that was not enough, make the result vary from image capture
	% or get(hImage,'CData'). Someone should really be punished for this.
	try
		if (handles.vecGraph)           % PS, etc ...
			ctype = 'img';
			if (~handles.imgOnly),   ctype = 'imgAx';      end
			imcapture(handles.hCallingFig,ctype,fname,get(handles.edit_mag,'String'));
			% We are done here, bye bye
			set(handsStBar(2:end),'Visible','on')
			return
		else                            % RASTER FORMATS
			captura = true;
			if (handles.imgOnly)            % Image only
				if ( handles.imgIsClean && get(handles.checkbox_origSize,'Val') && (handles.currMag == 1) )
					% Here we don't need to do any screen capture which forces img to be RGB
					img = im;
					captura = false;
				else
					img = imcapture(handles.hCallingFig,'img',dims);
				end
			else                            % Image and frames
				img = imcapture(handles.hCallingFig,'imgAx',get(handles.edit_mag,'String'));
			end
			if (handles.noname)             % We are in the output img (modal) mode
				if (captura && strncmp(get(handles.hCallingAx,'YDir'),'nor',3))
					img = flipdim(img,1);
				end
				handles.output = img;
				guidata(handles.figure1,handles)
				uiresume(handles.figure1);
				set(handsStBar(2:end),'Visible','on')
				return
			end
			if (~captura && strncmp(get(handles.hCallingAx,'YDir'),'nor',3))    % For imwrite
				img = flipdim(img,1);
			end
		end
	catch
		errordlg(lasterr,'Error');    return
	end
	set(handsStBar(2:end),'Visible','on')

	if ( strcmpi(EXT,'.jpg') || strcmpi(EXT,'.jpeg') )
		if (ndims(img) == 2),   img = ind2rgb8(img,get(handles.hCallingFig,'Colormap'));    end
		imwrite(img,fname,'Quality',handles.quality);
	elseif (strcmpi(EXT,'.raw'))
		if (ndims(img) == 2),   img = ind2rgb8(img,get(handles.hCallingFig,'Colormap'));    end
		fid = fopen(fname,'wb');
		[nl,nc,np] = size(img);                 l = 1;      pix = repmat(uint8(0),nl*nc*3,1);
		m = nl:-1:1;
		if (~strcmp(get(handles.hCallingAx,'Ydir'),'normal')),    m = 1:nl;     end
		for (i = m)
			for (j = 1:nc)
				for (k=1:3),	pix(l) = img(i,j,k);    l = l + 1;     end
			end
		end
		fwrite(fid,pix,'uint8');        fclose(fid);
	elseif strcmpi(EXT,'.gif')              % Non existent in < R14
		if (islogical(img)),        writegif(img,fname);
		else                        writegif(img,get(handles.hCallingFig,'Colormap'),fname);
		end
	elseif (strcmpi(EXT,'.png') || strcmpi(EXT,'.tif'))
		if (ndims(img) == 2)
			if (islogical(img)),    imwrite(img,fname);
			else                    imwrite(img,get(handles.hCallingFig,'Colormap'),fname);
			end
		else
			imwrite(img,fname);
		end
	else
		imwrite(img,fname);
	end
	delete(handles.figure1)

% -----------------------------------------------------------------------------
function push_cancel_CB(hObject, handles)
	figure1_CloseRequestFcn(hObject, [])

% -----------------------------------------------------------------------------
function [fname,ext] = stripExt(fname)
% Remove extension from filename FNAME and optionaly return the extension as well
	ind = max(find(fname == '.'));
	if (~isempty(ind))
		fname(ind:end) = [];
		if (nargout == 2)
			ext = fname(ind:end);
		end
	end

%-------------------------------------------------------------------------------------
% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata)
	handles = guidata(hObject);
	if (exist('OCTAVE_VERSION','builtin'))		% To know if we are running under Octave
		do_uiresume = ( isprop(hObject, '__uiwait_state__') && strcmp(get(hObject, '__uiwait_state__'), 'active') );
	else
		do_uiresume = strcmp(get(handles.figure1, 'waitstatus'), 'waiting');
	end
	if (do_uiresume)		% The GUI is still in UIWAIT, us UIRESUME
		guidata(handles.figure1, handles);	uiresume(handles.figure1);
	else					% The GUI is no longer waiting, just close it
		delete(handles.figure1);
	end

% --- Creates and returns a handle to the GUI figure. 
function snapshot_LayoutFcn(h1)

set(h1,'PaperUnits',get(0,'defaultfigurePaperUnits'),...
'Position',[520 625 427 175],...
'CloseRequestFcn',@figure1_CloseRequestFcn,...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'MenuBar','none',...
'Name','Snapshot',...
'NumberTitle','off',...
'RendererMode','manual',...
'Resize','off',...
'HandleVisibility','Call',...
'Tag','figure1');

uicontrol('Parent',h1, 'Position',[70 150 331 20],...
'BackgroundColor',[1 1 1],...
'Call',@main_uiCB,...
'HorizontalAlignment','left',...
'String','',...
'Style','edit',...
'Tag','edit_fname');

uicontrol('Parent',h1, 'Position',[400 149 24 21],...
'Call',@main_uiCB,...
'FontSize',12,...
'FontWeight','bold',...
'String','...',...
'Tag','push_outFile');

uicontrol('Parent',h1, 'Position',[10 153 60 15],...
'HorizontalAlignment','left',...
'String','Filename',...
'Style','text',...
'Tag','text1');

uicontrol('Parent',h1, 'Position',[10 128 70 15],...
'HorizontalAlignment','left',...
'String','Image Type',...
'Style','text',...
'Tag','text2');

uicontrol('Parent',h1, 'Position',[10 103 70 15],...
'HorizontalAlignment','left',...
'String','Image Size',...
'Style','text',...
'Tag','text_imgSize');

uicontrol('Parent',h1, 'Position',[10 77 65 15],...
'HorizontalAlignment','left',...
'String','Magnification',...
'Style','text',...
'Tag','text_mag');

uicontrol('Parent',h1, 'Position',[70 100 201 23],...
'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'),...
'Enable','inactive',...
'HorizontalAlignment','left',...
'Style','edit',...
'String','',...
'TooltipString','Estimated, uncompressed file size ... if RGB',...
'Tag','edit_imgSize');

uicontrol('Parent',h1, 'Position',[76 73 32 21],...
'BackgroundColor',[1 1 1],...
'Enable','inactive',...
'String','1',...
'Style','edit',...
'Tag','edit_mag');

uicontrol('Parent',h1, 'Position',[70 125 201 20],...
'BackgroundColor',[1 1 1],...
'Call',@main_uiCB,...
'Style','popupmenu',...
'String',{''},...
'Value',1,...
'Tag','popup_fileType');

uicontrol('Parent',h1, 'Position',[108 77 163 14],...
'BackgroundColor',[0.9 0.9 0.9],...
'Call',@main_uiCB,...
'Max',20,'Min',1,...
'Style','slider',...
'SliderStep',[0.05 0.1],...
'Val',1,...
'Tag','slider_mag');

uicontrol('Parent',h1, 'Position',[10 54 190 15],...
'Call',@main_uiCB,...
'String','Preserve Image original size',...
'Style','checkbox',...
'TooltipString','Use this option when you want that the output image has exactly the same same size as input',...
'Tag','checkbox_origSize');

uicontrol('Parent',h1, 'Position',[51 32 261 14],...
'BackgroundColor',[0.9 0.9 0.9],...
'Call',@main_uiCB,...
'Style','slider',...
'TooltipString','Higher numbers mean higher quality, and larger file size',...
'Value',0.75,...
'Tag','slider_quality');

uicontrol('Parent',h1, 'Position',[321 31 40 15],...
'String','75 %',...
'Style','text',...
'Tag','text_qualityLev');

uicontrol('Parent',h1, 'Position',[230 5 91 21],...
'Call',@main_uiCB,...
'String','Save',...
'Tag','push_save');

uicontrol('Parent',h1, 'Position',[332 5 91 21],...
'Call',@main_uiCB,...
'String','Cancel',...
'Tag','push_cancel');

uicontrol('Parent',h1, 'Position',[10 32 41 15],...
'HorizontalAlignment','left',...
'String','Quality',...
'Style','text',...
'Tag','text_Quality');

function main_uiCB(hObject, eventdata)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));

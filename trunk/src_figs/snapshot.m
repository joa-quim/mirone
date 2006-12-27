function varargout = snapshot(varargin)
% GUI interface to the imcapture function
%
%   snapshot(H) operates on the image contents of the UNIQUE image in Figure whose handle is H
%   snapshot(H,'whatever') as above but captures image and frame.
%
% M-File changed by desGUIDE 
%
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
snapshot_LayoutFcn(hObject,handles);
handles = guihandles(hObject);
movegui(hObject,'north')

if (~isempty(varargin))
    handles.hCallingFig = varargin{1};
    handlesMir = guidata(handles.hCallingFig);
else
    errordlg('SNAPSHOT: wrong number of arguments.','Error')
    delete(hObject);    return
end

if (handlesMir.no_file)
    errordlg('You didn''t even load a file. What are you expecting then?','Error')
    delete(hObject);    return
end

handles.imgOnly = 1;            % Flag that indicates pure image capture
if (numel(varargin) == 2)
    handles.imgOnly = 0;        % OK, so we will have a image + frame capture
end

% Add this figure handle to the carraças list
plugedWin = getappdata(handles.hCallingFig,'dependentFigs');
plugedWin = [plugedWin hObject];
setappdata(handles.hCallingFig,'dependentFigs',plugedWin);

% -------------- Find if we have any graphical objects ploted
ls = findobj(handles.hCallingFig,'Type','line');
ps = findobj(handles.hCallingFig,'Type','patch');
ts = findobj(handles.hCallingFig,'Type','text');
handles.imgIsClean = 0;
if (isempty(ls) && isempty(ps) && isempty(ts))
    handles.imgIsClean = 1;     % We won't need to screen capture if resolution is one-to-one
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
    'SUN rasterfile (*.ras)'; ...
    'Raw RGB format (*.raw)'};

handles.exts = {'.jpg' '.png' '.tif' '.bmp' '.eps' '.ai' '.emf' '.hdf' '.gif' '.pcx' '.ras' '.raw'};

% ------------ Get the image dimensions. Note that the dimensions are those deffined
%              by the axes size and not the true image size.
axUnit = get(handlesMir.axes1,'Units');         set(handlesMir.axes1,'Units','pixels')
axPos = get(handlesMir.axes1,'pos');            set(handlesMir.axes1,'Units',axUnit)
handles.imSize = size(get(handlesMir.hImg,'CData'));
handles.imAxSize = [axPos(4) axPos(3)];                     % Image size as reinterpolated to fit inside axes.

if (~handles.imgIsClean || numel(handles.imSize) == 3)      % Output will be RGB for sure. So remove formats
    ind = false(numel(str1),1);
    ind(9:10) = true;
    str1(ind) = [];       handles.exts(ind) = [];
end

% ---
set(handles.popup_fileType,'String',str1)

% ---------------- Get info to allow guessing imgAx image size
PU = get(handles.hCallingFig,'paperunits');     set(handles.hCallingFig,'paperunits','inch')
pp = get(handles.hCallingFig,'paperposition');  set(handles.hCallingFig,'paperunits',PU)
dpi = round(handles.imSize(2) / pp(3));
handles.pp = [pp(3) handles.imSize(1) / dpi];       % PaperPosition in inches as it will be set inside imcapture

% ---------------- Fill the edit image size with the default values
if (handles.imgOnly)
    nRows = handles.imAxSize(1);        nCols = handles.imAxSize(2);
    origMegs = handles.imSize(1)*handles.imSize(2) / 1048576;   % Original image size in Mb
    sizeOrigUnits = ' Megs';
    if (origMegs < 1)
        sizeOrigUnits = ' Kb';  origMegs = origMegs * 1024;
    end
    handles.txtOrigSize = ['(' sprintf('%d',handles.imSize(1)) 'x' sprintf('%d',handles.imSize(2)) ') ' ...
            sprintf('%.1f',origMegs) sizeOrigUnits];
else
    nRows = round(handles.pp(1) * 150);     nCols = round(handles.pp(2) * 150);
end
Megs = nRows * nCols * 3 / 1048576;
sizeUnits = ' Megs';
if (Megs < 1)               % Report file size in Kbytes
    sizeUnits = ' Kb';      Megs = Megs * 1024;
end
handles.txtThisSize = ['(' sprintf('%d',nRows) 'x' sprintf('%d',nCols) ') ' sprintf('%.1f',Megs) sizeUnits];
set(handles.edit_imgSize,'String',handles.txtThisSize,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));

% ---------------- Fill the edit file name with a default value
fname = get(handlesMir.figure1,'Name');
[pato,fname,ext] = fileparts(fname);
if (islogical(get(handlesMir.hImg,'CData')) && handles.imgOnly)   % Best proposition when we have a logical (mask) image
    fname = [handlesMir.work_dir filesep fname '.png'];
    set(handles.popup_fileType,'Val',2)
    set(handles.checkbox_origSize,'Val',1)
else
    fname = [handlesMir.work_dir filesep fname '.jpg'];
end
set(handles.edit_fname,'String',fname);

% ----------------- Set the slider with the apropriate range
if (handles.imgOnly)
    sliderRange(handles,1,1)
else
    set(handles.checkbox_origSize,'Enable','off')
    sliderRange(handles,0,150)
end

handles.hImg = handlesMir.hImg;
handles.hCallingAx = handlesMir.axes1;
handles.currMag  = 1;       % Current magnification
handles.quality = 75;       % Current quality level. Only applyes to the jpeg format
handles.currDPI  = 150;     % Current DPI for rasters
handles.currVecDPI = 300;   % Current DPI for vector graphics
handles.vecGraph = 0;       % To signal when we are dealing with vector graphics
handles.work_dir = handlesMir.work_dir;
handles.home_dir = handlesMir.home_dir;

% ----------------- Choose default command line output for snapshot_export
handles.output = hObject;
guidata(hObject, handles);

% UIWAIT makes snapshot_export wait for user response (see UIRESUME)
% uiwait(handles.figure1);


set(hObject,'Visible','on');
% NOTE: If you make uiwait active you have also to uncomment the next three lines
% handles = guidata(hObject);
% out = snapshot_OutputFcn(hObject, [], handles);
% varargout{1} = out;

% --- Outputs from this function are returned to the command line.
function varargout = snapshot_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

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
function edit_fname_Callback(hObject, eventdata, handles)
    fname = get(hObject,'String');
    [pato,fname,ext] = fileparts(fname);
    if (~strmatch(lower(ext),handles.exts))
        errordlg('You cannot choose a different file format than the ones offered to you.','ERROR')
        set(hObject,'String','')
    end

% -----------------------------------------------------------------------------
function push_outFile_Callback(hObject, eventdata, handles)
    contents = get(handles.popup_fileType,'String');
    val = get(handles.popup_fileType,'Val');    
    str = {handles.exts{val} contents{val}};
    cd(handles.work_dir)
    [FileName,PathName] = uiputfile(str,'Select file name');
    cd(handles.home_dir);       % allways go home
    if isequal(FileName,0);     return;     end
    [pato,fname,ext] = fileparts(FileName);
    if (~strcmp(lower(ext),handles.exts{val}))
        errordlg(['You cannot choose a different file format here. It has to be of type ' handles.exts{val}],'ERROR')
        return
    end
    set(handles.edit_fname,'String',[PathName FileName])

% -----------------------------------------------------------------------------
function popup_fileType_Callback(hObject, eventdata, handles)
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
        sliderRange(handles,magn,val);                      handles.vecGraph = 0;
    elseif (strmatch(ext,{'.ps' '.eps' '.ai' '.emf'}))
        set(handles.slider_quality,'Visible','off');        set(handles.text_Quality,'Visible','off')
        set(handles.text_qualityLev,'Visible','off');       set(handles.edit_imgSize,'String','Don''t know')
        set(handles.checkbox_origSize,'Enable','off')
        sliderRange(handles,0,handles.currVecDPI);          handles.vecGraph = 1;
    elseif ( ~get(handles.checkbox_origSize,'Val') && strcmp(ext,'.gif'))
        warndlg('Gif files do not support RGB. Check the "Preserve Image ..." checkbox before using this option.','Warning')
        set(handles.edit_fname,'String',[fname '.png']);    % Replace the censured gif by a png
        set(hObject,'Val',2)
        set(handles.slider_quality,'Visible','off');        set(handles.text_Quality,'Visible','off')
        set(handles.text_qualityLev,'Visible','off');       set(handles.edit_imgSize,'String',handles.txtThisSize);
        val = handles.currDPI;      magn = 0;               % Depending on raster capture mode the current dpi is different
        if (handles.imgOnly)
            val = handles.currMag;      magn = 1;
            set(handles.checkbox_origSize,'Enable','on')
        end
        sliderRange(handles,magn,val);                      handles.vecGraph = 0;
    else
        set(handles.slider_quality,'Visible','off');        set(handles.text_Quality,'Visible','off')
        set(handles.text_qualityLev,'Visible','off');       set(handles.edit_imgSize,'String',handles.txtThisSize);
        val = handles.currDPI;      magn = 0;               % Depending on raster capture mode the current dpi is different
        if (handles.imgOnly)
            val = handles.currMag;      magn = 1;
            set(handles.checkbox_origSize,'Enable','on')
        end
        sliderRange(handles,magn,val);                      handles.vecGraph = 0;
    end
    guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------
function slider_mag_Callback(hObject, eventdata, handles)
    % Callback to control the image magnetization factor slider
    mag = get(hObject,'Value');
    set(handles.edit_mag,'String',mag);
    if (~handles.vecGraph)
        if (handles.imgOnly)
            nRows = handles.imAxSize(1) * mag;      nCols = handles.imAxSize(2) * mag;
        else
            nRows = round(handles.pp(1) * mag);     nCols = round(handles.pp(2) * mag);
        end
        Megs = nRows * nCols * 3 / 1048576;
        sizeUnits = ' Megs';
        if (Megs < 1)               % Report file size in Kbytes
            sizeUnits = ' Kb';      Megs = Megs * 1024;
        end
        handles.txtThisSize = ['(' sprintf('%d',nRows) 'x' sprintf('%d',nCols) ') ' sprintf('%.1f',Megs) sizeUnits];
        set(handles.edit_imgSize,'String',handles.txtThisSize);
    end
    if (handles.vecGraph),          handles.currVecDPI = mag;
    elseif (handles.imgOnly),       handles.currMag = mag;
    else                            handles.currDPI = mag;
    end
    guidata(handles.figure1,handles)
    
% -----------------------------------------------------------------------------
function checkbox_origSize_Callback(hObject, eventdata, handles)
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
function slider_quality_Callback(hObject, eventdata, handles)
    handles.quality = round(get(hObject,'Value') * 100);
    set(handles.text_qualityLev,'String',[sprintf('%d',handles.quality) ' %']);
    guidata(handles.figure1,handles)

% -----------------------------------------------------------------------------
function push_save_Callback(hObject, eventdata, handles)
    try
        im = get(handles.hImg,'CData');
    catch
        errordlg('Something stupid occured. Probably on your side. Did you kill the Figure?','Error')
        return
    end
    
    fname = get(handles.edit_fname,'String');
    [pato, name, EXT] = fileparts(fname);
    if (isempty(name))
        errordlg('Saving in a no name! You think you are funny?','Error');        return
    end

    dims = handles.imAxSize * handles.currMag;
    if (get(handles.checkbox_origSize,'Val'))
        dims = handles.imSize(1:2);
    end
    try
        if (handles.vecGraph)           % PS, etc ...
            ctype = 'img';
            if (~handles.imgOnly),   ctype = 'imgAx';      end
            imcapture(handles.hCallingFig,ctype,fname,get(handles.edit_mag,'String'));
            % We are done here, bye bye
            return
        else                            % RASTER FORMATS
            if (handles.imgOnly)            % Image only
                if (get(handles.checkbox_origSize,'Val') && handles.imgIsClean)
                    % Here we don't need to do any screen capture which forces img to be RGB
                    img = im;
                else
                    img = imcapture(handles.hCallingFig,'img',dims);
                end
            else                            % Image and frames
                img = imcapture(handles.hCallingFig,'imgAx',get(handles.edit_mag,'String'));
            end
        end
    catch
        errordlg(lasterr,'Error');    return
    end

    if ( strcmpi(EXT,'.jpg') || strcmpi(EXT,'.jpeg') )
        if (ndims(img) == 2),   img = ind2rgb8(img,get(handles.hCallingFig,'Colormap'));    end
        imwrite(img,fname,'Quality',handles.quality);
    elseif (strcmpi(EXT,'.raw'))
        if (ndims(img) == 2),   img = ind2rgb8(img,get(handles.hCallingFig,'Colormap'));    end
        fid = fopen(fname,'wb');
        [nl,nc,np] = size(img);                 l = 1;      pix = repmat(uint8(0),nl*nc*3,1);
        m = nl:-1:1;
        if (~strcmp(get(handles.hCallingAx,'Ydir'),'normal')),    m = 1:nl;     end
        for i=m
            for j=1:nc
                for k=1:3;  pix(l) = img(i,j,k);    l = l + 1;     end
            end
        end
        fwrite(fid,pix,'uint8');        fclose(fid);
    elseif strcmpi(EXT,'.gif')              % Non existent in < R14
        if (strcmp(get(handles.hCallingAx,'Ydir'),'normal')),   img = flipud(img);      end     % That f... tendency to invertion
        if (islogical(img)),        writegif(img,fname);
        else                        writegif(img,get(handles.hCallingFig,'Colormap'),fname);
        end
    elseif (strcmpi(EXT,'.png') || strcmpi(EXT,'.tif'))
        if (ndims(img) == 2)
            if (strcmp(get(handles.hCallingAx,'Ydir'),'normal')),   img = flipud(img);   end    % That f... tendency to invertion
            if (islogical(img)),    imwrite(img,fname);
            else                    imwrite(img,get(handles.hCallingFig,'Colormap'),fname);
            end
        else
            imwrite(img,fname);
        end
    else
        imwrite(img,fname);
    end

% -----------------------------------------------------------------------------
function push_cancel_Callback(hObject, eventdata, handles)
    delete(handles.figure1)

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

% --- Creates and returns a handle to the GUI figure. 
function snapshot_LayoutFcn(h1,handles);

set(h1,'PaperUnits',get(0,'defaultfigurePaperUnits'),...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'MenuBar','none',...
'Name','Snapshot',...
'NumberTitle','off',...
'PaperPosition',get(0,'defaultfigurePaperPosition'),...
'PaperSize',[20.98404194812 29.67743169791],...
'Position',[520 625 427 175],...
'RendererMode','manual',...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1');

h2 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@snapshot_uicallback,h1,'edit_fname_Callback'},...
'HorizontalAlignment','left',...
'Position',[70 150 331 20],...
'Style','edit',...
'Tag','edit_fname');

h3 = uicontrol('Parent',h1,...
'Callback',{@snapshot_uicallback,h1,'push_outFile_Callback'},...
'FontSize',12,...
'FontWeight','bold',...
'Position',[400 149 24 21],...
'String','...',...
'Tag','push_outFile');

h4 = uicontrol('Parent',h1,...
'HorizontalAlignment','left',...
'Position',[10 153 51 15],...
'String','Filename',...
'Style','text',...
'Tag','text1');

h5 = uicontrol('Parent',h1,...
'HorizontalAlignment','left',...
'Position',[10 128 58 15],...
'String','Image Type',...
'Style','text',...
'Tag','text2');

h6 = uicontrol('Parent',h1,...
'HorizontalAlignment','left',...
'Position',[10 103 55 15],...
'String','Image Size',...
'Style','text',...
'Tag','text_imgSize');

h7 = uicontrol('Parent',h1,...
'HorizontalAlignment','left',...
'Position',[10 77 65 15],...
'String','Magnification',...
'Style','text',...
'Tag','text_mag');

h8 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Enable','inactive',...
'HorizontalAlignment','left',...
'Position',[70 100 201 20],...
'Style','edit',...
'TooltipString','Estimated, uncompressed file size ... if RGB',...
'Tag','edit_imgSize');

h9 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Enable','inactive',...
'Position',[76 73 32 21],...
'String','1',...
'Style','edit',...
'Tag','edit_mag');

h10 = uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@snapshot_uicallback,h1,'popup_fileType_Callback'},...
'Position',[70 125 201 20],...
'Style','popupmenu',...
'Value',1,...
'Tag','popup_fileType');

h11 = uicontrol('Parent',h1,...
'BackgroundColor',[0.9 0.9 0.9],...
'Callback',{@snapshot_uicallback,h1,'slider_mag_Callback'},...
'Max',20,'Min',1,...
'Position',[108 77 163 14],...
'Style','slider',...
'SliderStep',[0.05 0.1],...
'Tag','slider_mag');

h12 = uicontrol('Parent',h1,...
'Callback',{@snapshot_uicallback,h1,'checkbox_origSize_Callback'},...
'Position',[10 54 160 15],...
'String','Preserve Image original size',...
'Style','checkbox',...
'TooltipString','Use this option when you want that the output image has exactly the same same size as input',...
'Tag','checkbox_origSize');

h13 = uicontrol(...
'Parent',h1,...
'BackgroundColor',[0.9 0.9 0.9],...
'Callback',{@snapshot_uicallback,h1,'slider_quality_Callback'},...
'Position',[51 32 261 14],...
'Style','slider',...
'TooltipString','Higher numbers mean higher quality, and larger file size',...
'Value',0.75,...
'Tag','slider_quality');

h14 = uicontrol('Parent',h1,...
'Position',[321 31 40 15],...
'String','75 %',...
'Style','text',...
'Tag','text_qualityLev');

h15 = uicontrol('Parent',h1,...
'Callback',{@snapshot_uicallback,h1,'push_save_Callback'},...
'Position',[230 5 91 23],...
'String','Save',...
'Tag','push_save');

h16 = uicontrol('Parent',h1,...
'Callback',{@snapshot_uicallback,h1,'push_cancel_Callback'},...
'Position',[332 5 91 23],...
'String','Cancel',...
'Tag','push_cancel');

h17 = uicontrol('Parent',h1,...
'HorizontalAlignment','left',...
'Position',[10 32 41 15],...
'String','Quality',...
'Style','text',...
'Tag','text_Quality');

function snapshot_uicallback(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
feval(callback_name,hObject,[],guidata(h1));

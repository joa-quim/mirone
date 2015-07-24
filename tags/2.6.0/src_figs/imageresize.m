function varargout = imageresize(varargin)
% Helper Window to do image resizing

%	Copyright (c) 2004-2013 by J. Luis
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

	if (isempty(varargin)),		return,		end

	hObject = figure('Vis','off');
	imageresize_LayoutFcn(hObject);
	handles = guihandles(hObject);

	handles.hMirFig = varargin{1};

	% Get the Mirone handles. We need it here
	handMir = guidata(handles.hMirFig);
	if (handMir.no_file)
		errordlg('You didn''t even load a file. What are you expecting then?','Error')
		delete(hObject);    return
	end

	move2side(handles.hMirFig, hObject,'right')

	handles.hImage = handMir.hImg;
	handles.imgSize = size(get(handles.hImage,'CData'));

	handles.head = handMir.head;
	handles.DefLineThick = handMir.DefLineThick;
	handles.DefLineColor = handMir.DefLineColor;
	handles.geog = handMir.geog;
	handles.image_type = handMir.image_type;
	if (numel(handles.imgSize) == 2)
		handles.parentCmap = get(handMir.figure1,'Colormap');
	end

	handles.pixWidth  = handles.imgSize(2);
	handles.pixHeight = handles.imgSize(1);
	handles.resolution = 72;
	handles.unitFact  = 2.54;
	handles.resolutionFact = 1 / handles.resolution * handles.unitFact;
	handles.constrainProp = 1;
	handles.intepMethod = 'nearest';
	handles.isPercent = 0;
	handles.IamCompiled = handMir.IamCompiled;

	% Fill in the necessary uicontrols
	set(handles.edit_pixHeight,'String',handles.imgSize(1))
	set(handles.edit_pixWidth,'String',handles.imgSize(2))
	set(handles.edit_docResolution,'String',handles.resolution)
	handles = pix2size(handles,'w');           % Fill the "Document" edits
	handles = pix2size(handles,'h');

	%------------ Give a Pro look (3D) to the frame boxes  --------
	new_frame3D(hObject, [handles.text_docSize handles.text_pixDim]);
	%------------- END Pro look (3D) ------------------------------

	% Add this figure handle to the carraças list
	plugedWin = getappdata(handMir.figure1,'dependentFigs');
	plugedWin = [plugedWin hObject];
	setappdata(handMir.figure1,'dependentFigs',plugedWin);

	set(hObject,'Vis','on');
	if (nargout),		varargout{1} = hObject;		end
	guidata(hObject, handles);

% ----------------------------------------------------------------------------------
function hand = pix2size(handles,opt)
% Convert from pixels to "Document Size" unities
	if ( opt == 'w' )       % Width
		handles.docWidth = handles.pixWidth * handles.resolutionFact;
		set(handles.edit_docWidth,'String',handles.docWidth)
	else                    % Height
		handles.docHeight = handles.pixHeight * handles.resolutionFact;
		set(handles.edit_docHeight,'String',handles.docHeight)
	end
	guidata(handles.figure1, handles);
	if (nargout),   hand = handles;     end

% ----------------------------------------------------------------------------------
function handles = size2pix(handles,opt)
% Convert from "Document Size" to pixels unities
	if ( opt == 'w' )       % Width
		handles.pixWidth = handles.docWidth / handles.resolutionFact;
		set(handles.edit_pixWidth,'String',round(handles.pixWidth))
	else                    % Height
		handles.pixHeight = handles.docHeight / handles.resolutionFact;
		set(handles.edit_pixHeight,'String',round(handles.pixHeight))
	end
	guidata(handles.figure1, handles);

% ----------------------------------------------------------------------------------
function edit_pixWidth_CB(hObject, handles)
    xx = round( str2double(get(hObject,'String')) );        p = xx / 100;
    if (isnan(xx)),     set(hObject,'String',handles.pixWidth);    return;      end
    if (handles.isPercent),     xx = handles.imgSize(2) * p;       end
    pixWidthOld = handles.pixWidth;
    handles.pixWidth  = xx;
    pix2size(handles,'w')           % Update the "Document" edit
    if (handles.constrainProp)      % Recompute height
        if (~handles.isPercent)
            handles.pixHeight = handles.pixHeight * handles.pixWidth / pixWidthOld;
            set(handles.edit_pixHeight,'String',round(handles.pixHeight))
        else
            handles.pixHeight = handles.imgSize(1) * p;
            set(handles.edit_pixHeight,'String',round(p*100))
        end
        pix2size(handles,'h')       % Update the "Document" edit
    end
    guidata(handles.figure1, handles);

% ----------------------------------------------------------------------------------
function edit_pixHeight_CB(hObject, handles)
    xx = round( str2double(get(hObject,'String')) );        p = xx / 100;
    if (handles.isPercent),     xx = round(handles.imgSize(1) * xx / 100);      end
    if (isnan(xx)),     set(hObject,'String',handles.pixHeight);    return;     end
    pixHeightOld = handles.pixHeight;
    handles.pixHeight  = xx;
    pix2size(handles,'h')           % Update the "Document" edit
    if (handles.constrainProp)      % Recompute height
        if (~handles.isPercent)
            handles.pixWidth = handles.pixWidth * handles.pixHeight / pixHeightOld;
            set(handles.edit_pixWidth,'String',round(handles.pixWidth))
        else
            handles.pixWidth = handles.imgSize(2) * p;
            set(handles.edit_pixWidth,'String',round(p*100))
        end
        pix2size(handles,'w')       % Update the "Document" edit
    end
    guidata(handles.figure1, handles);

% ----------------------------------------------------------------------------------
function edit_docWidth_CB(hObject, handles)
    xx = round( str2double(get(hObject,'String')) );
    if (isnan(xx)),     set(hObject,'String',handles.docWidth);    return;     end
    docWidthOld = handles.docWidth;
    handles.docWidth  = xx;
    handles = size2pix(handles,'w');% Update the "Document" edit
    if (handles.constrainProp)      % Recompute height
        handles.docHeight = handles.docHeight * handles.docWidth / docWidthOld;
        set(handles.edit_docHeight,'String',handles.docHeight)
        size2pix(handles,'h');      % Update the "Document" edit
    end

% ----------------------------------------------------------------------------------
function edit_docHeight_CB(hObject, handles)
    xx = round( str2double(get(hObject,'String')) );
    if (isnan(xx)),     set(hObject,'String',handles.docHeight);    return;     end
    docHeightOld = handles.docHeight;
    handles.docHeight  = xx;
    handles = size2pix(handles,'h');% Update the "Pixel" edit
    if (handles.constrainProp)      % Recompute height
        handles.docWidth = handles.docWidth * handles.docHeight / docHeightOld;
        set(handles.edit_docWidth,'String',handles.docWidth)
        size2pix(handles,'w');      % Update the "Pixel" edit
    end

% ----------------------------------------------------------------------------------
function edit_docResolution_CB(hObject, handles)
% Note that even if the resolution is set to pixels/cm we store in pixels/inch
    xx = str2double(get(hObject,'String'));
    if (isnan(xx)),     set(hObject,'String',handles.resolution);    return;     end
    if (get(handles.popup_docResolution,'Value') == 1)     % pixel/inch
        handles.resolution  = round(xx);
    else
        handles.resolution  = round(xx*2.54);
    end
    handles.resolutionFact = 1 / handles.resolution * handles.unitFact;
    handles = size2pix(handles,'w');        % Update the "Pixel" edits
    size2pix(handles,'h');                  % handles is updated here

% ----------------------------------------------------------------------------------
function popup_pixWidth_CB(hObject, handles)
% Change between pixels and percent, but only on the display
	val = get(hObject,'Value');
	switch val
		case 1
			w = handles.pixWidth;            h = handles.pixHeight;
			handles.isPercent = 0;
		case 2
			w = handles.pixWidth / handles.imgSize(2) * 100;
			h = handles.pixHeight / handles.imgSize(1) * 100;
			handles.isPercent = 1;
	end
	set(handles.edit_pixWidth,'String',round(w))
	set(handles.edit_pixHeight,'String',round(h))
	set(handles.popup_pixHeight,'Value',val)
	guidata(handles.figure1,handles)

% ----------------------------------------------------------------------------------
function popup_pixHeight_CB(hObject, handles)
% Let the popup_pixWidth do the work
    set(handles.popup_pixWidth,'Value',get(hObject,'Value'))
    popup_pixWidth_CB(handles.popup_pixWidth, handles)

% ----------------------------------------------------------------------------------
function popup_docWidth_CB(hObject, handles)
	switch get(hObject,'Value')
		case 1,     handles.unitFact  = 2.54;
		case 2,     handles.unitFact  = 254;
		case 3,     handles.unitFact  = 1;
		case 4,     handles.unitFact  = 1 / 72;
	end
	set(handles.popup_docHeight,'Value',get(hObject,'Value'))
	handles.resolutionFact = 1 / handles.resolution * handles.unitFact;
	guidata(handles.figure1,handles)
	pix2size(handles,'w')
	pix2size(handles,'h')
    
% ----------------------------------------------------------------------------------
function popup_docHeight_CB(hObject, handles)
% Let the popup_docWidth do the work
    set(handles.popup_docWidth,'Value',get(hObject,'Value'))
    popup_pixWidth_CB(handles.popup_docWidth, handles)

% ----------------------------------------------------------------------------------
function popup_docResolution_CB(hObject, handles)
% In fact we never change the resolution unites out of DPI
	switch get(hObject,'Value')
		case 1,     resolution = handles.resolution;
		case 2,     resolution = handles.resolution / 2.54;
	end
	set(handles.edit_docResolution,'String',resolution)

% ----------------------------------------------------------------------------------
function check_constProportions_CB(hObject, handles)
	if (get(hObject,'Value'))
		handles.constrainProp = 1;
	else
		handles.constrainProp = 0;
	end
	guidata(handles.figure1,handles)

% ----------------------------------------------------------------------------------
function popup_resampMethod_CB(hObject, handles)
% Change the interpolation method
	switch get(hObject,'Value')
		case 1,    handles.intepMethod = 'nearest';
		case 2,    handles.intepMethod = 'bilinear';
		case 3,    handles.intepMethod = 'bicubic';
		case 4,    handles.intepMethod = 'area';
	end
	guidata(handles.figure1,handles)

% ----------------------------------------------------------------------------------
function push_OK_CB(hObject, handles)
% Do the job

	img = get(handles.hImage,'CData');
	set(handles.figure1,'pointer','watch');
	if (handles.IamCompiled && (strcmp(handles.intepMethod,'bicubic') || strcmp(handles.intepMethod,'bilinear')) )
		img = img_fun('imresize',img,round([handles.pixHeight handles.pixWidth]),handles.intepMethod);
	else
		img = cvlib_mex('resize',img,round([handles.pixHeight handles.pixWidth]),handles.intepMethod);
	end
	set(handles.figure1,'pointer','arrow');

	if (handles.image_type == 2)
		if (ndims(img) == 2)
			setappdata(0,'CropedColormap',handles.parentCmap);
		end
		h = mirone(img);
		set(h, 'Name', 'Resized Image')
	else
		tmp.head = [handles.head(1:6) 0];
		tmp.head(8) = diff(handles.head(1:2)) / (round(handles.pixWidth) - ~handles.head(7));
		tmp.head(9) = diff(handles.head(3:4)) / (round(handles.pixHeight) - ~handles.head(7));
		tmp.X = handles.head(1:2);      tmp.Y = handles.head(3:4);
		tmp.name = 'Resized Image';
		if (ndims(img) == 2)
			tmp.cmap = handles.parentCmap;
		end
		mirone(img,tmp);
	end
    
% ----------------------------------------------------------------------------------
function push_cancel_CB(hObject, handles)
	delete(handles.figure1)

% --- Creates and returns a handle to the GUI figure. 
function imageresize_LayoutFcn(h1)

set(h1,...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'MenuBar','none',...
'Name','Resize Image',...
'NumberTitle','off',...
'Position',[520 539 328 261],...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1');

uicontrol('Parent',h1,'Position',[7 177 241 75],'Style','frame');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@imageresize_uiCB,...
'Position',[79 214 71 21],...
'Style','edit',...
'Tooltip','Describes the width in pixels','Tag','edit_pixWidth');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@imageresize_uiCB,...
'Position',[79 187 71 21],...
'Style','edit',...
'Tooltip','Describes the height in pixels','Tag','edit_pixHeight');

uicontrol('Parent',h1,'FontSize',9,...
'Position',[29 219 41 15],...
'String','Width',...
'Style','text',...
'Tooltip','Describes the width in pixels');

uicontrol('Parent',h1,'FontSize',9,...
'Position',[29 190 41 18],...
'String','Height',...
'Style','text',...
'Tooltip','Describes the height in pixels');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@imageresize_uiCB,...
'Position',[158 214 81 22],...
'String',{'pixels'; 'percent'},...
'Style','popupmenu',...
'Value',1,...
'Tag','popup_pixWidth');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@imageresize_uiCB,...
'Position',[158 187 81 22],...
'String',{'pixels'; 'percent'},...
'Style','popupmenu',...
'Value',1,...
'Tag','popup_pixHeight');

uicontrol('Parent',h1,'FontSize',9,...
'Position',[28 242 120 17],...
'String',' Pixel Dimensions:',...
'Style','text','Tag','text_pixDim');

uicontrol('Parent',h1,'Position',[7 60 241 100],'Style','frame');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@imageresize_uiCB,...
'Position',[79 122 71 21],...
'Style','edit',...
'Tooltip','Set the document width','Tag','edit_docWidth');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@imageresize_uiCB,...
'Position',[79 95 71 21],...
'Style','edit',...
'Tooltip','Set the document height','Tag','edit_docHeight');

uicontrol('Parent',h1,'FontSize',9,...
'Position',[30 127 41 15],...
'String','Width',...
'Style','text',...
'Tooltip','Set the document width','Tag','text4');

uicontrol('Parent',h1,'FontSize',9,...
'Position',[30 98 41 18],...
'String','Height',...
'Style','text',...
'Tooltip','Set the document height','Tag','text5');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@imageresize_uiCB,...
'Position',[159 122 81 22],...
'String',{'cm'; 'mm'; 'inch'; 'points' },...
'Style','popupmenu',...
'Value',1,'Tag','popup_docWidth');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@imageresize_uiCB,...
'Position',[159 95 81 22],...
'String',{'cm'; 'mm'; 'inch'; 'points' },...
'Style','popupmenu',...
'Value',1,'Tag','popup_docHeight');

uicontrol('Parent',h1,'FontSize',9,...
'Position',[29 150 110 17],...
'String',' Document Size:',...
'Style','text','Tag','text_docSize');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@imageresize_uiCB,...
'Position',[79 69 71 21],...
'Style','edit',...
'Tooltip','Set the document resolution','Tag','edit_docResolution');

uicontrol('Parent',h1,'FontSize',9,...
'Position',[13 71 65 18],...
'String','Resolution',...
'Style','text',...
'Tooltip','Set the document resolution');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@imageresize_uiCB,...
'Position',[158 68 83 22],...
'String',{  'pixels/inch'; 'pixels/cm' },...
'Style','popupmenu',...
'Tooltip','Set the document resolution',...
'Value',1,'Tag','popup_docResolution');

uicontrol('Parent',h1,...
'Call',@imageresize_uiCB,...
'FontSize',9,'Position',[7 36 151 16],...
'String','Constrain Proportions',...
'Style','checkbox',...
'Tooltip','Constrain aspect ratio',...
'Value',1,'Tag','check_constProportions');

uicontrol('Parent',h1,'BackgroundColor',[1 1 1],...
'Call',@imageresize_uiCB,...
'Position',[118 10 135 22],...
'String',{'Nearest Neighbor'; 'Bilinear'; 'Bicubic'; 'Area Relation' },...
'Style','popupmenu',...
'Value',1,'Tag','popup_resampMethod');

uicontrol('Parent',h1,...
'Call',@imageresize_uiCB,...
'FontSize',10,...
'Position',[257 229 66 21],...
'String','OK','Tag','push_OK');

uicontrol('Parent',h1,...
'Call',@imageresize_uiCB,...
'FontSize',10,...
'Position',[257 189 66 21],...
'String','Cancel','Tag','push_cancel');

uicontrol('Parent',h1,'FontSize',9,...
'Position',[7 11 106 18],...
'String','Resample Method:',...
'Style','text');

function imageresize_uiCB(hObject, eventdata)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));

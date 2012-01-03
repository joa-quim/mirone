function varargout = rotatetool(varargin)
% ROTATETOOL is an interface to the IMROTATE function.  
%
%   ROTATETOOL(HANDLE) rotates a single image contained in HANDLE.
%   HANDLE can be a handle to a figure, axes, or image object.
%
%   To retrieve rotated image use the GETIMAGE function on the appropriate
%   figure window.
%
% USAGE:
%   [OUT] = ROTATETOOL(GRD), where GRD is a MxN of any class "above" UINT8 rotates the GRD
%   array and returns it in OUT. The GRD image's representation is computed with a linear
%   scaling on the [0 255] interval.
%       ( img = uint8(round( ((GRD - GRD_z_min) / (GRD_z_max - GRD_z_min))*255 )) )
%
%   [OUT,HDR_OUT] = ROTATETOOL(GRD,HDR_IN), where HDR_IN is a header row vector with the GRD
%   coordinate referencing, outputs also an updated version of the header. This header has the
%   the form used by the GMT MEX programs, which is:
%       [x_min x_max y_min y_max z_min z_max registratrion x_inc y_inc]
%       registratrion is either 0 (for grid registration) or 1 to pixel registration.
%       Interested users that want to know more must consult the GMT manual.
%   NOTE: This is not the ML way of working with grids, but it's ALOT more efficient than the
%   standard use of MESHGRID which implies triplicating the memory use.
%
%   [OUT] = ROTATETOOL(IMG), IMG is a uint8 indexed or RGB array rotates the IMG
%   array and returns it in OUT.
%
%   [OUT] = ROTATETOOL(IMG,CMAP), uses the colormap CMAP in the preview image. If CMAP is
%   not given in input and the image is indexed, a cmap = jet(256) is used.
%
%   [OUT,HDR_OUT] = ROTATETOOL(IMG,HDR_IN,CMAP), IMG is a uint8 indexed or RGB array. Same as
%   above, but with a header info vector like in the GRD case. Use this for example if the
%   IMG comes from a georeferenced image.
%
%   [...] = ROTATETOOL(GRD|IMG,HANDLES), GRD & IMG as above and HANDLES is a Mirone handles
%
%   The "Apply n return" button performs the rotation on the input data and returns
%   the output if it was requested.
%
%   The "Cancel" button just closes this figure returns nothing
%
% EXAMPLES:
%   I = checkerboard(50);
%   h = imshow(I);
%   rotatetool(h)   % When satisfied push 'Apply n return' then call GETIMAGE
%   B = getimage(gcf);
%
%   Operate on a double array, and get it rotated
%   Z = peaks(128);
%   grd = rotatetool(Z,jet(256));
%
%   Let us pretend that the Z obove is a referenced grid. So contruct first the HDR_IN vector
%   hdr_in = [1 128 1 128 min(Z(:)) max(Z(:)) 0 1 1];
%   [grd,hdr] = rotatetool(Z,hdr_in);
%
%   Demonstration mode
%   rotatetool      % displays 'peppers.png' and opens GUI.
%
%   THIS GUI NEEDS THE IMAGE PROCESSING TOOLBOX TO RUN
%
%   CREDITS
%       This GUI is a pure imitation of ROTATEGUI from Birju Patel and uses some code from
%       that function. However, it has extended functionalities and was re-writen using R13
%       syntax. Which means that it can be compiled with the V3.x compiler.
%

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
	rotatetool_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hObject,'right')
	handles.IAmOctave = (exist('OCTAVE_VERSION','builtin') ~= 0);	% To know if we are running under Octave

	% Some initializations
	handles.frame_axesPos = get(handles.frame_axes,'Pos');
	handles.cancel = 0;				% It will be set to one if the cancel button is hit or the fig is killed
	handles.angle = 0;
	handles.do_flipLR = false;
	handles.do_flipUD = false;
	handles.hOrigFigure = [];
	handles.hOrigAxes = [];
	handles.hOrigImage = [];
	handles.OrigGrd = [];
	handles.in_hdr = [];			% To hold an eventual input GMT style row vector header
	handles.out_header = [];		% To hold an eventual output GMT style row vector header
	cmap = [];

	n_argout = nargout;     n_argin = numel(varargin);
	if (n_argin >= 1)
		if (ishandle(varargin{1}))			% First argument is an image handle, or ... error
			if ( strcmp(get(varargin{1},'type'),'axes') || strcmp(get(varargin{1},'type'),'fig'))
				handles.hOrigImage = findobj(varargin{1},'type','image');
			elseif ( strcmp(get(varargin{1},'type'),'image') )
				handles.hOrigImage = varargin{1};
			end
			msg = [];
			if (isempty(handles.hOrigImage) || numel(handles.hOrigImage) > 1)
				msg = 'ERROR: None or more than one image objects found in the figure.';
			end
			if (n_argout)
				msg = 'ERROR: Output arguments are not allowed when the input is a handle.';
			end
			if (~isempty(msg))
				errordlg(msg,'ERROR');      delete(hObject);    return
			end
			handles.hOrigAxes = get(handles.hOrigImage,'Parent');
			handles.hOrigFigure = get(handles.hOrigAxes,'Parent');
		else
			if (isa(varargin{1},'uint8') || isa(varargin{1},'logical'))			% First argument is an image array
				handles.OrigImage = varargin{1};
			else                                    % First arg is a "grid"
				handles.OrigGrd = varargin{1};
				handles.OrigImage = scaleto8(handles.OrigGrd);  % Compute grid's image
				cmap = jet(256);
			end
		end

		if (n_argin == 2 && isstruct(varargin{2}))
			% This is a new form (as on 6-3-07) ROTATETOOL(GRD|IMG,HANDLES)
			handles.hOrigAxes = varargin{2}.axes1;
			handles.in_hdr = varargin{2}.head;
			cmap = get(varargin{2}.figure1,'Colormap');

			% Add this Fig to the carra?as list
			plugedWin = getappdata(varargin{2}.figure1,'dependentFigs');
			plugedWin = [plugedWin hObject];
			setappdata(varargin{2}.figure1,'dependentFigs',plugedWin);
		else
			if (n_argin == 2 && size(varargin{2},1) == 1)       % Second argument is a GMT style vector header
				handles.in_hdr = varargin{2};
			elseif (n_argin == 2 && size(varargin{2},2) == 3)   % Second argument is a colormap
				cmap = varargin{2};
			end
			if (n_argin == 3)                           % There is no choice in args order here
				handles.in_hdr = varargin{2};
				cmap = varargin{3};
			end
		end
		% One last test on CMAP & HEADER likelihood
		msg = [];
		if (~isempty(handles.in_hdr) && size(handles.in_hdr,2) < 4)
			msg = 'ERROR: The GMT header vector is invalid';
		end
		if (~isempty(cmap) && size(cmap,2) ~= 3)
			msg = 'ERROR: Invalid colormap in imput';
		end
		if (~isempty(msg))
			errordlg(msg,'ERROR');      delete(hObject);    return
		end
	else        % Demo mode
		h = figure;
		handles.hOrigImage = imshow('peppers.png');
		handles.hOrigAxes = get(handles.hOrigImage,'Parent');
		handles.hOrigFigure = h;
	end

	% Determine resize needed so that preivew image fits in defined area
	% ReductionFactor determines how to scale original image so
	% that it fits in the defined preview area.
	pImageConstr = min(handles.frame_axesPos(3),handles.frame_axesPos(4));    %Constraining dimension in image panel
	if (~isempty(handles.hOrigImage))       % A image handle was transmited in argument
		[oriImageHeight,oriImageWidth,k] = size(get(handles.hOrigImage,'Cdata'));
	else
		[oriImageHeight,oriImageWidth,k] = size(handles.OrigImage);
	end
	%Length of Image diagonal in pixels. +5 added to allow for some panel border effects
	dImage = ceil(sqrt(oriImageHeight^2 + oriImageWidth^2)) + 5;
	reductionFactor = (pImageConstr/dImage);

	if (~isempty(handles.hOrigImage))       % A image handle was transmited in argument
		handles.previewImage = img_fun('imresize',get(handles.hOrigImage,'Cdata'),reductionFactor);
	else
		handles.previewImage = img_fun('imresize',handles.OrigImage,reductionFactor);
	end
	[imHeight,imWidth,k] = size(handles.previewImage);

	% Center axes containing preview image in the "image panel" so rotation will look like a pin-wheel
	[leftPosition,bottomPosition] = calcPosition(handles,imWidth,imHeight);
		% Make the position relative to frame_axes LL corner
	leftPosition = leftPosition + handles.frame_axesPos(1);
	bottomPosition = bottomPosition + handles.frame_axesPos(2);
	set(handles.axes1,'Pos',[leftPosition bottomPosition imWidth imHeight],...
		'DataAspectRatio', [1 1 1], 'PlotBoxAspectRatioMode', 'auto')

	%display image in axes
	handles.himage = displayPreviewImage(handles,handles.previewImage,cmap);

	% Store the f. ydir mess for be taken account later
	if (strcmp(get(handles.axes1,'YDir'),'normal'))
		handles.y_dir = 1;
	else
		handles.y_dir = -1;
	end

	%------------ Give a Pro look (3D) to the frame boxes  -------------------------------
		new_frame3D(hObject, NaN)
	%------------- END Pro look (3D) -----------------------------------------------------

	% Choose default command line output for rotatetool
	handles.output = hObject;
	if (n_argout == 1 && isempty(handles.OrigGrd))		% This the only case where the Fig handle can go out
		varargout{1} = hObject;
	end
	if (n_argout == 2)
		handles.out_header = 1;
		set(handles.popup_bboxmenu,'Enable','off')		% Don't want the head-ache of computing new limits
	end
	handles.n_argout = n_argout;
	guidata(hObject, handles);
	set(hObject,'Vis','on');

	% Check if we have to do a uiwait and what goes out (if anything)
	if (n_argout > 0)						% If called with output, we must wait
		uiwait(hObject);
		handles = guidata(hObject);
		if (handles.cancel)					% Don't try to output eventual non-existing variables
			if (n_argout == 1),     varargout{1} = [];      end
			if (n_argout == 2),     varargout{1} = [];      varargout{2} = [];  end
			delete(handles.figure1);		% The figure can be deleted now
			return
		end
		varargout{1} = handles.output_grd;
		if (n_argout == 2)
			varargout{2} = handles.out_header;
		end
		delete(handles.figure1);			% The figure can be deleted now
	end

% ----------------------------------------------------------------------------
function push_cw90_CB(hObject, handles)
    % When the rotation button is pushed repeatedly, the rotation
    % applied to the image is cumulative.  Pushing the cw button n
    % times will rotate preview image 90*n degrees.
    updateAngleInfo(handles,-90);
    updateImage(handles);

% ----------------------------------------------------------------------------
function push_ccw90_CB(hObject, handles)
    updateAngleInfo(handles,90);
    updateImage(handles);

% ----------------------------------------------------------------------------
function push_flipLR_CB(hObject, handles)
    newCdata = flipdim(get(handles.himage,'CData'),2);
    set(handles.himage,'Cdata',newCdata);
    handles.do_flipLR = ~handles.do_flipLR;     % For the Apply button
    guidata(handles.figure1,handles)

% ----------------------------------------------------------------------------
function push_flipUD_CB(hObject, handles)
    newCdata = flipdim(get(handles.himage,'CData'),1);
    set(handles.himage,'Cdata',newCdata);
    handles.do_flipUD = ~handles.do_flipUD;
    guidata(handles.figure1,handles)

% ----------------------------------------------------------------------------
function slider_CB(hObject, handles)
    % Affordance for manual rotation. User can drag slider to rotate preview
    % image in small increments. The slider boundaries are [-180, 180].
    % Text box next to slider will show values within this range.
    handles.angle = get(hObject,'Value');
    guidata(handles.figure1,handles)
    updateImage(handles);
    set(handles.edit_angRot,'String',num2str(handles.angle));

% ----------------------------------------------------------------------------
function edit_angRot_CB(hObject, handles)
    % User can enter in a single scalar number (rotation angle in degrees) and image
    % will be rotated by that amount.  If the input angle is out of the [-180,180]
    % range, the value will be wrapped into this range using the checkAngle function.
    textvalue = str2double(get(hObject,'String'));
    % Expression not allowe in text field.
    if (~isempty(textvalue) && ~isnan(textvalue) && isreal(textvalue))
        handles.angle = checkAngle(textvalue);
        guidata(handles.figure1,handles)
        updateImage(handles);
        set(handles.slider,'Value',handles.angle);
        set(hObject,'String',num2str(handles.angle));
    else
        errordlg('Text field must contain a single scalar value.')
    end

% ----------------------------------------------------------------------------
function popup_bboxmenu_CB(hObject, handles)
    % changes bounding box option for IMROTATE
    updateImage(handles);

% ----------------------------------------------------------------------------
function popup_interpmenu_CB(hObject, handles)
    % changes to interp method won't be noticeable in preview image so the preview
    % image will always use nearest interp. But we need to record the selected
    % option in case the user exits without further rotations to image
    val = get(hObject,'Value');
    switch val
        case 1,     handles.interpmethod = 'nearest';
        case 2,     handles.interpmethod = 'bilinear';
        case 3,     handles.interpmethod = 'bicubic';
    end
    guidata(handles.figure1,handles)

% ----------------------------------------------------------------------------
function push_cancel_CB(hObject, handles)
	if (~handles.IAmOctave)
		do_uiresume = strcmp(get(handles.figure1, 'waitstatus'), 'waiting');
	else
		do_uiresume = ( isprop(hObject, '__uiwait_state__') && strcmp(get(hObject, '__uiwait_state__'), 'active') );
	end
	if (do_uiresume)		% The GUI is still in UIWAIT, us UIRESUME
		handles.cancel = 1;		% User gave up, return nothing
		guidata(handles.figure1, handles);    uiresume(handles.figure1);
	else					% The GUI is no longer waiting, just close it
		delete(handles.figure1);
	end

% ----------------------------------------------------------------------------
function push_apply_CB(hObject, handles)
    switch get(handles.popup_interpmenu,'Value')
        case 1,     interpmethod = 'nearest';
        case 2,     interpmethod = 'bilinear';
        case 3,     interpmethod = 'bicubic';
    end
    if (~isempty(handles.hOrigFigure))  % Original data is in a figure
        A = get(handles.hOrigImage,'CData');
        if (handles.do_flipLR),     A = flipdim(A,2);   end
        if (handles.do_flipUD),     A = flipdim(A,1);   end
        rotCdata = transform_fun('imrotate',A,-handles.angle*handles.y_dir,interpmethod,handles.bbox);
        set(handles.hOrigImage,'CData',rotCdata);
    else
        if (~isempty(handles.OrigGrd))
            A = handles.OrigGrd;
            if (handles.do_flipLR),     A = flipdim(A,2);   end
            if (handles.do_flipUD),     A = flipdim(A,1);   end
            rotCdata = transform_fun('imrotate',A,-handles.angle*handles.y_dir,interpmethod,'loose');
        else
            A = handles.OrigImage;
            if (handles.do_flipLR),     A = flipdim(A,2);   end
            if (handles.do_flipUD),     A = flipdim(A,1);   end
            rotCdata = transform_fun('imrotate',A,-handles.angle*handles.y_dir,interpmethod,handles.bbox);
        end
        handles.output_grd = rotCdata;
        if (handles.out_header)
            nx = size(rotCdata,2);      ny = size(rotCdata,1);
            x_inc = (handles.in_hdr(2) - handles.in_hdr(1)) / (nx - ~handles.in_hdr(7));
            y_inc = (handles.in_hdr(4) - handles.in_hdr(3)) / (ny - ~handles.in_hdr(7));
            handles.out_header = [handles.in_hdr(1:4) double(min(rotCdata(:))) ...
                double(max(rotCdata(:))) handles.in_hdr(7) x_inc y_inc];
        end
        guidata(handles.figure1,handles)
    end

	if (~handles.IAmOctave)
		do_uiresume = strcmp(get(handles.figure1, 'waitstatus'), 'waiting');
	else
		do_uiresume = ( isprop(hObject, '__uiwait_state__') && strcmp(get(hObject, '__uiwait_state__'), 'active') );
	end
	if (do_uiresume)		% The GUI is still in UIWAIT, us UIRESUME
		uiresume(handles.figure1);
	else					% The GUI is no longer waiting, just close it
		delete(handles.figure1);
	end

%=============================
function updateAngleInfo(handles,angle)
    % used for the +/- 90 rotation callbacks, which rotate current
    % preview image in increments of +/- 90 degrees
    updatedAngle = checkAngle(angle + handles.angle);
    handles.angle = updatedAngle;
    set(handles.edit_angRot,'String',num2str(updatedAngle));
    set(handles.slider,'Value',updatedAngle);
    guidata(handles.figure1,handles)

%=======================================
function updatedAngle = checkAngle(angle)
    % used to wrap rotation angle between +/- 180
    updatedAngle = angle;
    if ~(abs(angle) <= 180) % wrap angle if necessary
        r = rem(angle,360);
        if (abs(r) <= 180),     updatedAngle = r;
        else                    updatedAngle = r-sign(r)*360;
        end
    end

%==============================================
function [leftpos,bottompos] = calcPosition(handles,w,h)
    % calculate position values to center axes(image) in preview area.
    % Rotation will look like a pin wheel.
    PanelPos = handles.frame_axesPos;
    leftpos   = (PanelPos(3) - w)/2;
    bottompos = (PanelPos(4) - h)/2;

%==============================================================
function [hout] = displayPreviewImage(handles,newCdata,cmap)
    % Create a preview image from original image.

    [h,w,k] = size(newCdata);    xdata = [1 w];    ydata = [1 h];

    if (~isempty(handles.hOrigImage))
        clim = get(get(handles.hOrigImage,'Parent'),'Clim');
        cdatamapping = get(handles.hOrigImage,'CdataMapping');
        cmap = get(handles.hOrigFigure,'Colormap');             % Those will not overwrite the input because
        y_dir = get(get(handles.hOrigImage,'Parent'),'YDir');   % cmap is not know yet (that is: cmap = [])
        hout = image(xdata,ydata,newCdata, ...
            'BusyAction', 'cancel', 'CDataMapping', cdatamapping, ...
            'Parent', handles.axes1, 'Interruptible', 'off');
    else
        hout = image(xdata,ydata,newCdata, 'BusyAction', 'cancel', ...
            'Parent', handles.axes1, 'Interruptible', 'off');
        clim = [];
        if (~isempty(handles.hOrigAxes))
            y_dir = get(handles.hOrigAxes,'YDir');
        else
            y_dir = 'reverse';      % Comes here when called with only grd|img and no info on gca
        end
    end

    % Set axes and figure properties if necessary to display the image object correctly.
    axesPosition = get(handles.axes1,'Position');
    set(handles.axes1, 'Xlim',[1 w], 'Ylim',[1 h], 'Ydir', y_dir, ...
        'Position',[axesPosition(1:2) w h], 'Visible','off');

    if (~isempty(cmap)),    set(handles.figure1, 'Colormap', cmap);   end
    if (~isempty(clim)),    set(handles.axes1, 'CLim', clim);   end

    isIndexedUint16Image = strcmpi(get(hout,'CDataMapping'),'direct') && size(cmap,1) > 256;

    if (isIndexedUint16Image && ispc)
        set(handles.figure1,'Renderer','Zbuffer');
    end
        
%=========================
function updateImage(handles)
    % updateImage performs the rotates preview image
    handles = guidata(handles.figure1);
    interpmethod = 'nearest';                       % changes to the interp method will
    switch get(handles.popup_bboxmenu,'Value')      % not be noticeable in the preview image.
        case 1,     outsize = 'loose';
        case 2,     outsize = 'crop';
    end
    %assuming cw rotation is positive (need to use -angle because
    %IMROTATE assumes cw rotation is negative)
    rotCdata = transform_fun('imrotate',handles.previewImage,-handles.angle*handles.y_dir,interpmethod,outsize);
    [h,w,k] = size(rotCdata);
    [left,bottom] = calcPosition(handles,w,h);
    left = left + handles.frame_axesPos(1);
    bottom = bottom + handles.frame_axesPos(2);
    %set(handles.axes1, 'Xlim',[1 w], 'Ylim',[1 h], 'Position',[left bottom w h])
    set(handles.axes1, 'Position',[left bottom w h])
    set(handles.himage,'Cdata',rotCdata)
    %Store rotation specific information needed to rotate the original
    %image after the rotate GUI is closed or the user hits apply.

    handles.bbox = outsize;
    guidata(handles.figure1,handles)

% ----------------------------------------------------------------------------
function figure1_CloseRequestFcn(hObject, eventdata)
handles = guidata(hObject);
push_cancel_CB(hObject, handles)

% ----------------------------------------------------------------------------
function figure1_KeyPressFcn(hObject, eventdata)
if isequal(get(hObject,'CurrentKey'),'escape')
	handles = guidata(hObject);
	push_cancel_CB(hObject, handles)
end

% ----------------------------------------------------------------------------
% --- Creates and returns a handle to the GUI figure. 
function rotatetool_LayoutFcn(h1)

set(h1,...
'PaperUnits',get(0,'defaultfigurePaperUnits'),...
'CloseRequestFcn',@figure1_CloseRequestFcn,...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'DoubleBuffer','on',...
'KeyPressFcn',@figure1_KeyPressFcn,...
'MenuBar','none',...
'Name','rotatetool',...
'NumberTitle','off',...
'Position',[520 418 451 375],...
'RendererMode','manual',...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1');

uicontrol('Parent',h1,'Position',[134 40 314 108],'Style','frame');
uicontrol('Parent',h1,'Position',[2 4 446 35],'Style','frame');
uicontrol('Parent',h1,'Position',[2 40 131 332],'Style','frame');

uicontrol('Parent',h1,...
'Call',@rotatetool_uiCB,...
'Position',[7 339 120 20],...
'String',['Clockwise 90' char(186)],...
'Tag','push_cw90');

axes('Parent',h1,...
'Units','pixels',...
'CameraPosition',[0.5 0.5 9.16025403784439],...
'CameraPositionMode',get(0,'defaultaxesCameraPositionMode'),...
'Position',[135 150 311 221],...
'Tag','axes1');

uicontrol('Parent',h1,...
'Call',@rotatetool_uiCB,...
'Position',[7 309 120 20],...
'String',['Counterclockwise 90' char(186)],...
'Tag','push_ccw90');

uicontrol('Parent',h1,...
'Call',@rotatetool_uiCB,...
'Position',[7 279 120 20],...
'String','Flip Left/Right',...
'Tag','push_flipLR');

uicontrol('Parent',h1,...
'Call',@rotatetool_uiCB,...
'Position',[7 249 120 20],...
'String','Flip Up/Down',...
'Tag','push_flipUD');

uicontrol('Parent',h1,...
'Call',@rotatetool_uiCB,...
'Max',180,...
'Min',-180,...
'SliderStep',[1/360 5/360],...
'Position',[150 122 181 18],...
'Style','slider',...
'Tag','slider');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@rotatetool_uiCB,...
'Position',[350 120 71 21],...
'String','0',...
'Style','edit',...
'Tag','edit_angRot');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@rotatetool_uiCB,...
'Position',[210 79 231 22],...
'String',{'Expanded to Fit Rotated Input Image'; 'Same as Input Image'},...
'Style','popupmenu', 'Value',1,...
'Tag','popup_bboxmenu');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Call',@rotatetool_uiCB,...
'Position',[210 50 140 22],...
'String',{'Nearest-neighbor'; 'Bilinear'; 'Bicubic'},...
'Style','popupmenu', 'Value',1,...
'Tag','popup_interpmenu');

uicontrol('Parent',h1,'Position',[141 83 67 15],...
'String','Output Size:','Style','text');

uicontrol('Parent',h1,'Position',[141 54 68 15],...
'String','Interpolation:','Style','text');

uicontrol('Parent',h1,...
'Call',@rotatetool_uiCB,...
'Position',[262 9 60 21],...
'String','Cancel',...
'Tag','push_cancel');

uicontrol('Parent',h1,...
'Call',@rotatetool_uiCB,...
'Position',[341 9 100 21],...
'String','Apply n Return',...
'Tooltip','Apply the rotation to original data and quit',...
'Tag','push_apply');

uicontrol('Parent',h1,...
'HorizontalAlignment','left',...
'Position',[147 106 30 15],...
'String',['-180' char(186)],...
'Style','text');

uicontrol('Parent',h1,...
'HorizontalAlignment','left',...
'Position',[309 106 30 15],...
'String',['180' char(186)],...
'Style','text');

uicontrol('Parent',h1,...
'HorizontalAlignment','left',...
'Position',[241 106 10 15],...
'String',['0' char(186)],...
'Style','text');

uicontrol('Parent',h1,'Position',[134 149 314 223], 'Style','frame', 'Tag','frame_axes');

function rotatetool_uiCB(hObject, eventdata)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));

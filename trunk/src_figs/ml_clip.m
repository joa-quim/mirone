function varargout = ml_clip(varargin)
% Do grid cippling based on user selected values.
%
% ml_clip(handles) 
%		where handles is the Mirone handles do the clipping and display it on a separate figure
%
% [Zout, head] = ml_clip(handles, Z)
%		Apply the cipping on the array Z and assume that handles.head holds the array's header
%		On output the head(5:6) is updated to new ZOUT min/max
%		Figure is set to 'modal' when this option is used.

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

	if isempty(varargin)
		errordlg('GRDCLIP: wrong number of arguments.','Error'),	return
	end

	handMir = varargin{1};

	if (handMir.no_file)
		errordlg('GRDCLIP: You didn''t even load a file. What are you expecting then?','ERROR')
		return
	end
	if (~handMir.validGrid)
        errordlg('GRDCLIP: This operation is deffined only for images derived from DEM grids.','ERROR')
		return
	end

	hObject = figure('Vis','off');
	ml_clip_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(handMir.figure1, hObject)

	if (nargin == 1)
		handles.Z = getappdata(handMir.figure1,'dem_z');
	else
		handles.Z = varargin{2};		% Operate on a given array
	end
	handles.have_nans = handMir.have_nans;

	if (isempty(handles.Z))
        errordlg('GRDCLIP: Grid was not saved in memory. Increase "Grid max size" and start over.','ERROR')
        delete(hObject);    return
	end

    handles.hMirFig = handMir.figure1;
	handles.version6 = ~handMir.version7;
	handles.have_nans = handMir.have_nans;
	handles.head = handMir.head;
	handles.above_val = [];
	handles.below_val = [];
    handles.z_min = handles.head(5);
    handles.z_max = handles.head(6);
    handles.above = handles.z_max;
    handles.below = handles.z_min;
    
    set(handles.edit_above,'String',sprintf('%.4g',handles.z_max))
    set(handles.edit_below,'String',sprintf('%.4g',handles.z_min))

	% Add this figure handle to the carraças list
	plugedWin = getappdata(handles.hMirFig,'dependentFigs');
	plugedWin = [plugedWin hObject];
	setappdata(handles.hMirFig,'dependentFigs',plugedWin);

	%------------ Give a Pro look (3D) to the frame boxes  -------------------------------
	new_frame3D(hObject, handles.text_statHammer, handles.frame1)
	%------------- END Pro look (3D) -----------------------------------------------------

	guidata(hObject, handles);
	set(hObject,'Visible','on');
	
	if (nargin == 2)
		set(hObject,'WindowStyle','modal')
		uiwait(hObject);			% Waita aí
		handles = guidata(hObject);
		if (~isempty(handles.head))
			varargout{1} = handles.Z;		varargout{2} = handles.head;
		else
			[varargout{1:nargout}] = deal([]);
		end
		delete(hObject)
	elseif (nargout)
		varargout{1} = hObject;
	end

% -------------------------------------------------------------------------------------
function edit_above_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	if (~isnan(xx) && xx < handles.z_max),	handles.above = xx;
	else									set(hObject,'String',handles.z_max);
	end
	guidata(hObject,handles)

% -------------------------------------------------------------------------------------
function edit_Ab_val_CB(hObject, handles)
	handles.above_val = str2double(get(hObject,'String'));
	guidata(handles.figure1,handles)

% -------------------------------------------------------------------------------------
function edit_below_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	if (~isnan(xx) && xx > handles.z_min),	handles.below = xx;
	else									set(hObject,'String',handles.z_min);
	end
	guidata(hObject,handles)

% -------------------------------------------------------------------------------------
function edit_Bl_val_CB(hObject, handles)
	handles.below_val = str2double(get(hObject,'String'));
	guidata(handles.figure1,handles)

% -------------------------------------------------------------------------------------
function push_OK_CB(hObject, handles)
% Do the actual work

	% The following is due to unbelievable BUG in R14 that does not know how to deal correctly with single(NaNs)
	if (handles.version6 && handles.have_nans)
		indNaN = isnan(handles.Z);
	end

	if ( ~isempty(handles.above_val) && ~isempty(handles.below_val) && ...
			(handles.above_val < handles.below) | (handles.below_val > handles.above) ) %#ok (NEED that |)
		% Need special care to not clip the already clipped values
		ind1 = handles.Z > handles.above;
		ind2 = handles.Z < handles.below;
		handles.Z(ind1) = handles.above_val;	clear ind1
		handles.Z(ind2) = handles.below_val;	clear ind2
	else
		if ~isempty(handles.above_val)				% Clip above
			handles.Z(handles.Z > handles.above) = handles.above_val;
		end
		if ~isempty(handles.below_val)				% Clip below
			handles.Z(handles.Z < handles.below) = handles.below_val;
		end
	end

	if (handles.version6 && handles.have_nans)			% Shame on you TMW
		handles.Z(indNaN) = single(NaN);		clear indNaN
	end
	zz = grdutils(handles.Z,'-L');       handles.head(5:6) = zz(1:2);

	if (strcmp(get(handles.figure1,'WindowStyle'),'modal'))		% In this mode we just return the clipped array (in handles)
		guidata(handles.figure1, handles)
		uiresume(handles.figure1)
		return
	end

    tmp.X = linspace(handles.head(1),handles.head(2),size(handles.Z,2));
    tmp.Y = linspace(handles.head(3),handles.head(4),size(handles.Z,1));
    tmp.head = handles.head;
    tmp.name = 'Clipped grid';
    mirone(handles.Z,tmp);
	if (get(handles.check_delFig, 'Val'))
		delete(handles.figure1)
	end

% -------------------------------------------------------------------------------------
function edit_percent_CB(hObject, handles)
	xx = str2double(get(hObject, 'String'));
	if (isnan(xx)),		set(hObject, 'String', ''),		return,		end
	set([handles.edit_nSigma handles.edit_mad], 'String', '')

% -------------------------------------------------------------------------------------
function edit_nSigma_CB(hObject, handles)
	xx = str2double(get(hObject, 'String'));
	if (isnan(xx)),		set(hObject, 'String', ''),		return,		end
	set([handles.edit_percent handles.edit_mad], 'String', '')

% -------------------------------------------------------------------------------------
function edit_mad_CB(hObject, handles)
	xx = str2double(get(hObject, 'String'));
	if (isnan(xx)),		set(hObject, 'String', ''),		return,		end
	set([handles.edit_percent handles.edit_nSigma], 'String', '')

% -------------------------------------------------------------------------------------
function push_okUP_CB(hObject, handles)
% ...
	up = [];
	xx = abs(str2double(get(handles.edit_percent, 'String'))) * 0.01;
	if (~isnan(xx))
		if (handles.have_nans)
			ind = ~isnan(handles.Z);
			s = sort(handles.Z(ind(:)));
		else
			s = sort(handles.Z(:));
		end
		n_out = round(numel(s) * xx/2);
		low = s(n_out);		up = s(numel(s) - n_out);
	end

	xx = abs(str2double(get(handles.edit_nSigma, 'String')));
	if (~isnan(xx))
		med_std = grdutils(handles.Z, '-S');
		media = double(med_std(1));
		stdv = double(med_std(2));
		low = media - xx * stdv;	up = media + xx * stdv;
	end
	xx = abs(str2double(get(handles.edit_mad, 'String')));
	if (~isnan(xx))
		z = handles.Z(:);
		if (handles.have_nans),		z(isnan(z)) = [];	end
		med = double(median(z));
		z = cvlib_mex('absDiffS', z, med);
		z = sort(z);
		n = numel(z);
		if rem(n,2) == 1
			mad = 1.4826 * double(z((n+1)/2));
		else
			mad = 1.4826 * 0.5 * ( double(z(n/2)) + double(z(n/2+1)) );
		end
		low = med - xx * mad;		up = med + xx * mad;
	end
	
	if (isempty(up)),	return,		end			% No entries in any of the option boxes
	set(handles.edit_above, 'String', up)
	set(handles.edit_below, 'String', low)
	handles.above = up;
	handles.below = low;
	guidata(handles.figure1,handles)

% -------------------------------------------------------------------------------------
function pushbutton_help_CB(hObject, handles)

% -------------------------------------------------------------------------------------
function m = mad(handles, x)
% MAD Median Absolute Deviation
%	MAD (A) returns the robust estimate of the deviation
%	about the median.

% 	x = sort(abs(x - median(x)));
	med = median(x);
	x = cvlib_mex('SubS', x, double(med));
	cvlib_mex('abs', x);
	x = sort(x);
	n = length(x);
	if rem(n,2) == 1
		m = 1.4826 * double(x((n+1)/2));
	else
		n2 = n / 2;
		m = 1.4826 * 0.5 * ( double(x(n2)) + double(x(n2+1)) );
	end

% -------------------------------------------------------------------------------------
function figure1_KeyPressFcn(hObject, eventdata)
	if isequal(get(hObject,'CurrentKey'),'escape')
		delete(hObject);
	end

% -------------------------------------------------------------------------------------
function figure1_CloseRequestFcn(hObject, eventdata)
	if (strcmp(get(hObject,'WindowStyle'),'modal'))
		handles = guidata(hObject);
		handles.head = [];			% Secret signal that user aborted
		guidata(handles.figure1, handles)
		uiresume(handles.figure1)
	end

% -------------------------------------------------------------------------------------
% --- Creates and returns a handle to the GUI figure. 
function ml_clip_LayoutFcn(h1)
set(h1, 'Position',[520 612 285 185],...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'KeyPressFcn',@figure1_KeyPressFcn,...
'CloseRequestFcn',@figure1_CloseRequestFcn,...
'MenuBar','none',...
'Name','Clipp Grid',...
'NumberTitle','off',...
'Resize','off',...
'Tag','figure1');

uicontrol('Parent',h1, 'Position',[4 145 75 21],...
'BackgroundColor',[1 1 1],...
'Call',@ml_clip_uiCB,...
'HorizontalAlignment','left',...
'Style','edit',...
'Tooltip','Grid nodes higher than this will be replaced "Value"',...
'Tag','edit_above');

uicontrol('Parent',h1, 'Position',[82 145 55 21],...
'BackgroundColor',[1 1 1],...
'Call',@ml_clip_uiCB,...
'HorizontalAlignment','left',...
'Style','edit',...
'Tooltip','Grid nodes > "Above" will be replaced by this value (''NaN'' is a valid string)',...
'Tag','edit_Ab_val');

uicontrol('Parent',h1, 'Position',[147 145 75 21],...
'BackgroundColor',[1 1 1],...
'Call',@ml_clip_uiCB,...
'HorizontalAlignment','left',...
'Style','edit',...
'Tooltip','Grid nodes lower than this will be replaced "Value"',...
'Tag','edit_below');

uicontrol('Parent',h1, 'Position',[224 145 55 21],...
'BackgroundColor',[1 1 1],...
'Call',@ml_clip_uiCB,...
'HorizontalAlignment','left',...
'Style','edit',...
'Tooltip','Grid nodes < "Below" will be replaced by this value (''NaN'' is a valid string)',...
'Tag','edit_Bl_val');

uicontrol('Parent',h1, 'Position',[21 166 41 15], 'FontName','Helvetica',...
'String','Above', 'Style','text');

uicontrol('Parent',h1, 'Position',[88 166 41 15], 'FontName','Helvetica',...
'String','Value', 'Style','text');

uicontrol('Parent',h1, 'Position',[166 166 41 15], 'FontName','Helvetica',...
'String','Below', 'Style','text');

uicontrol('Parent',h1, 'Position',[230 166 41 15], 'FontName','Helvetica',...
'String','Value', 'Style','text');

uicontrol('Parent',h1, 'Position',[4 116 175 21],...
'Style','checkbox',...
'Val', 1,...
'String','Delete Window at end',...
'ToolTip', 'When checked this window is deleted after hitting "Apply"',...
'Tag','check_delFig');

uicontrol('Parent',h1, 'Position',[215 117 66 21],...
'Call',@ml_clip_uiCB,...
'FontName','Helvetica',...
'FontSize',10,...
'String','Apply',...
'Tag','push_OK');

uicontrol('Parent',h1,'Position',[0 104 285 3],'Style','frame','Tag','frame1');

uicontrol('Parent',h1, 'Position',[4 64 40 21],...
'BackgroundColor',[1 1 1],...
'Call',@ml_clip_uiCB,...
'Style','edit',...
'Tooltip',sprintf('These percentage of points on the lower and upper\nZ limits value are selected for clipping'),...
'Tag','edit_percent');

uicontrol('Parent',h1, 'Position',[45 67 85 16],...
'FontName','Helvetica',...
'String','% End members',...
'HorizontalAlignment','left',...
'Tooltip',sprintf('These percentage of points on the lower and upper\nZ limits value are selected for clipping'),...
'Style','text');

uicontrol('Parent',h1, 'Position',[4 34 40 21],...
'BackgroundColor',[1 1 1],...
'Call',@ml_clip_uiCB,...
'Style','edit',...
'Tooltip', 'Pick up limits based on mean ± n*STD',...
'Tag','edit_nSigma');

uicontrol('Parent',h1, 'Position',[46 38 51 15],...
'FontName','Helvetica',...
'String','n STD',...
'HorizontalAlignment','left',...
'Tooltip', 'Pick up limits based on mean ± n*STD',...
'Style','text');

uicontrol('Parent',h1, 'Position',[4 4 40 21],...
'BackgroundColor',[1 1 1],...
'Call',@ml_clip_uiCB,...
'Style','edit',...
'Tooltip', 'Pick up limits based on median ± n*MAD',...
'Tag','edit_mad');

uicontrol('Parent',h1, 'Position',[46 8 51 15],...
'FontName','Helvetica',...
'HorizontalAlignment','left',...
'String','n MAD',...
'Tooltip', 'Pick up limits based on median ± n*MAD',...
'Style','text');

uicontrol('Parent',h1, 'Position',[221 3 50 75],...
'Call',@ml_clip_uiCB,...
'FontName','Helvetica',...
'FontSize',10,...
'String','UP',...
'Tooltip','Translate into clipping values  ',...
'Tag','push_okUP');

uicontrol('Parent',h1, 'Position',[60 97 150 18],...
'FontAngle','italic',...
'FontName','Helvetica',...
'FontSize',10,...
'String','Statistical Hammering',...
'Style','text',...
'Tag','text_statHammer');

function ml_clip_uiCB(hObject, eventdata)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));

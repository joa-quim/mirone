function varargout = ml_clip(varargin)
% M-File changed by desGUIDE 
% varargin   command line arguments to ml_clip (see VARARGIN)

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
	ml_clip_LayoutFcn(hObject);
	handles = guihandles(hObject);
	movegui(hObject,'north')
 
	if ~isempty(varargin)
		handMir  = varargin{1};
		handles.Z = getappdata(handMir.figure1,'dem_z');
		handles.have_nans = handMir.have_nans;
	else
        errordlg('GRDCLIP: wrong number of arguments.','Error')
        delete(hObject);    return
	end
    
	if (handMir.no_file)
		errordlg('GRDCLIP: You didn''t even load a file. What are you expecting then?','ERROR')
        delete(hObject);    return
	end
	if (~handMir.validGrid)
        errordlg('GRDCLIP: This operation is deffined only for images derived from DEM grids.','ERROR')
        delete(hObject);    return
	end
	if (isempty(handles.Z))
        errordlg('GRDCLIP: Grid was not saved in memory. Increase "Grid max size" and start over.','ERROR')
        delete(hObject);    return
	end

    handles.hMirFig = handMir.figure1;
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
	bgcolor = get(0,'DefaultUicontrolBackgroundColor');
	framecolor = max(min(0.65*bgcolor,[1 1 1]),[0 0 0]);
	h_f = handles.frame1;
	for i = 1:numel(h_f)
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
	delete(handles.text_statHammer)		% Recreate this text. Better than call the stupid uistack function
	handles.text_statHammer = uicontrol('Parent',hObject,'FontName','Helvetica','FontSize',10,'Position',[60 97 150 18],...
		'FontAngle','italic', 'String','Statistical Hammering','Style','text');
	%------------- END Pro look (3D) -------------------------------------------------------

	guidata(hObject, handles);
	set(hObject,'Visible','on');
	if (nargout),   varargout{1} = hObject;     end

% -------------------------------------------------------------------------------------
function edit_above_Callback(hObject, eventdata, handles)
	xx = str2double(get(hObject,'String'));
	if (~isnan(xx) && xx < handles.z_max),	handles.above = xx;
	else									set(hObject,'String',num2str(handles.z_max));
	end
	guidata(hObject,handles)

% -------------------------------------------------------------------------------------
function edit_Ab_val_Callback(hObject, eventdata, handles)
	handles.above_val = str2double(get(hObject,'String'));
	guidata(handles.figure1,handles)

% -------------------------------------------------------------------------------------
function edit_below_Callback(hObject, eventdata, handles)
	xx = str2double(get(hObject,'String'));
	if (~isnan(xx) & xx > handles.z_min),	handles.below = xx;
	else									set(hObject,'String',num2str(handles.z_min));
	end
	guidata(hObject,handles)

% -------------------------------------------------------------------------------------
function edit_Bl_val_Callback(hObject, eventdata, handles)
	handles.below_val = str2double(get(hObject,'String'));
	guidata(handles.figure1,handles)

% -------------------------------------------------------------------------------------
function pushbutton_OK_Callback(hObject, eventdata, handles)
	Out{1} = handles.above;             Out{3} = handles.below;
	Out{2} = handles.above_val;         Out{4} = handles.below_val;
    
	if ~isempty(handles.above_val)     % Clip above
        handles.Z(handles.Z > handles.above) = handles.above_val;
    end
	if ~isempty(handles.below_val)     % Clip below
        handles.Z(handles.Z < handles.below) = handles.below_val;
    end

	zz = grdutils(handles.Z,'-L');       handles.head(5:6) = zz(1:2);
    tmp.X = linspace(handles.head(1),handles.head(2),size(handles.Z,2));
    tmp.Y = linspace(handles.head(3),handles.head(4),size(handles.Z,1));
    tmp.head = handles.head;
    tmp.name = 'Clipped grid';
    mirone(handles.Z,tmp);
    delete(handles.figure1)

% -------------------------------------------------------------------------------------
function edit_percent_Callback(hObject, eventdata, handles)
	xx = str2double(get(hObject, 'String'));
	if (isnan(xx)),		set(hObject, 'String', ''),		return,		end
	set([handles.edit_nSigma handles.edit_mad], 'String', '')

% -------------------------------------------------------------------------------------
function edit_nSigma_Callback(hObject, eventdata, handles)
	xx = str2double(get(hObject, 'String'));
	if (isnan(xx)),		set(hObject, 'String', ''),		return,		end
	set([handles.edit_percent handles.edit_mad], 'String', '')

% -------------------------------------------------------------------------------------
function edit_mad_Callback(hObject, eventdata, handles)
	xx = str2double(get(hObject, 'String'));
	if (isnan(xx)),		set(hObject, 'String', ''),		return,		end
	set([handles.edit_percent handles.edit_nSigma], 'String', '')

% -------------------------------------------------------------------------------------
function push_okUP_Callback(hObject, eventdata, handles)
	xx = abs(str2double(get(handles.edit_percent, 'String'))) * 0.01;
	if (~isnan(xx))
		s = sort(handles.Z(:));
		n_out = round(numel(s) * xx/2);
		low = s(n_out);		up = s(numel(s) - n_out);
	end
	xx = abs(str2double(get(handles.edit_nSigma, 'String')));
	if (~isnan(xx))
		med_std = grdutils(handles.Z, '-S');		media = double(med_std(1));		stdv = double(med_std(2));
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
	set(handles.edit_above, 'String', up)
	set(handles.edit_below, 'String', low)
	handles.above = up;
	handles.below = low;
	guidata(handles.figure1,handles)

% -------------------------------------------------------------------------------------
function pushbutton_help_Callback(hObject, eventdata, handles)

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
% --- Creates and returns a handle to the GUI figure. 
function ml_clip_LayoutFcn(h1);
set(h1,...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'KeyPressFcn',@figure1_KeyPressFcn,...
'MenuBar','none',...
'Name','Clipp Grid',...
'NumberTitle','off',...
'Position',[520 612 285 185],...
'Resize','off',...
'Tag','figure1');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@ml_clip_uicallback,h1,'edit_above_Callback'},...
'HorizontalAlignment','left',...
'Position',[4 145 75 21],...
'Style','edit',...
'TooltipString','Grid nodes higher than this will be replaced "Value"',...
'Tag','edit_above');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@ml_clip_uicallback,h1,'edit_Ab_val_Callback'},...
'HorizontalAlignment','left',...
'Position',[82 145 55 21],...
'Style','edit',...
'TooltipString','Grid nodes > "Above" will be replaced by this value (''NaN'' is a valid string)',...
'Tag','edit_Ab_val');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@ml_clip_uicallback,h1,'edit_below_Callback'},...
'HorizontalAlignment','left',...
'Position',[147 145 75 21],...
'Style','edit',...
'TooltipString','Grid nodes lower than this will be replaced "Value"',...
'Tag','edit_below');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@ml_clip_uicallback,h1,'edit_Bl_val_Callback'},...
'HorizontalAlignment','left',...
'Position',[224 145 55 21],...
'Style','edit',...
'TooltipString','Grid nodes < "Below" will be replaced by this value (''NaN'' is a valid string)',...
'Tag','edit_Bl_val');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'Position',[21 168 41 15],...
'String','Above',...
'Style','text');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'Position',[88 168 41 15],...
'String','Value',...
'Style','text');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'Position',[166 168 41 15],...
'String','Below',...
'Style','text');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'Position',[230 168 41 15],...
'String','Value',...
'Style','text');

uicontrol('Parent',h1,...
'Callback',{@ml_clip_uicallback,h1,'pushbutton_OK_Callback'},...
'FontName','Helvetica',...
'FontSize',10,...
'Position',[215 117 66 21],...
'String','Apply',...
'Tag','pushbutton_OK');

uicontrol('Parent',h1,'Position',[0 104 285 3],'Style','frame','Tag','frame1');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@ml_clip_uicallback,h1,'edit_percent_Callback'},...
'Position',[4 64 40 21],...
'Style','edit',...
'Tooltip',sprintf('These percentage of points on the lower and upper\nZ limits value are selected for clipping'),...
'Tag','edit_percent');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'Position',[45 67 85 16],...
'String','% End members',...
'HorizontalAlignment','left',...
'Tooltip',sprintf('These percentage of points on the lower and upper\nZ limits value are selected for clipping'),...
'Style','text');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@ml_clip_uicallback,h1,'edit_nSigma_Callback'},...
'Position',[4 34 40 21],...
'Style','edit',...
'Tooltip', 'Pick up limits based on mean ± n*STD',...
'Tag','edit_nSigma');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'Position',[46 38 51 15],...
'String','n STD',...
'HorizontalAlignment','left',...
'Tooltip', 'Pick up limits based on mean ± n*STD',...
'Style','text');

uicontrol('Parent',h1,...
'BackgroundColor',[1 1 1],...
'Callback',{@ml_clip_uicallback,h1,'edit_mad_Callback'},...
'Position',[4 4 40 21],...
'Style','edit',...
'Tooltip', 'Pick up limits based on median ± n*MAD',...
'Tag','edit_mad');

uicontrol('Parent',h1,...
'FontName','Helvetica',...
'HorizontalAlignment','left',...
'Position',[46 8 51 15],...
'String','n MAD',...
'Tooltip', 'Pick up limits based on median ± n*MAD',...
'Style','text');

uicontrol('Parent',h1,...
'Callback',{@ml_clip_uicallback,h1,'push_okUP_Callback'},...
'FontName','Helvetica',...
'FontSize',10,...
'Position',[221 3 50 75],...
'String','UP',...
'TooltipString','Translate into clipping values  ',...
'Tag','push_okUP');

uicontrol('Parent',h1,...
'FontAngle','italic',...
'FontName','Helvetica',...
'FontSize',10,...
'Position',[60 97 150 18],...
'String','Statistical Hammering',...
'Style','text',...
'Tag','text_statHammer');

function ml_clip_uicallback(hObject, eventdata, h1, callback_name)
% This function is executed by the callback and than the handles is allways updated.
feval(callback_name,hObject,[],guidata(h1));

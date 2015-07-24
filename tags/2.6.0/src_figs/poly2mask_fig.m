function varargout = poly2mask_fig(varargin)
% Helper figure to create a mask image from lines or polygons in the calling fig 

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

	if (numel(varargin) < 2)
		errordlg('POLY2MASK: wrong number of input arguments.','Error')
		return
	end
 
	hObject = figure('Vis','off');
	poly2mask_fig_LayoutFcn(hObject);
	handles = guihandles(hObject);
 
	handles.head = varargin{1}.head;
	handles.hMirAxes = varargin{1}.axes1;
	handles.hMirImg = varargin{1}.hImg;
	handles.image_type = varargin{1}.image_type;
	handles.hMirFig = varargin{1}.figure1;
	handles.have_nans = varargin{1}.have_nans;
	handles.hPoly_current = varargin{2};

	move2side(varargin{1}.figure1, hObject)

	handles.inputIsLine = false;	% This was the original & unique behavior but now we accept to mask plines too
	handles.grid_in  = 1;
	handles.grid_out = 0;

	if (strcmp(get(handles.hPoly_current, 'Type'), 'line'))		% If called "from a line" (we'll do a screencature)
		x = get(handles.hPoly_current,'XData');		y = get(handles.hPoly_current,'YData');
		if (numel(x) < 3 || x(1) ~= x(end) || y(1) ~= y(end) )
			handles.inputIsLine = true;
		end
	end
	if (~handles.inputIsLine)		% Get handles of closed objs (polygons & patches)
		handles.hPoly = find_closed(handles);
		if (numel(handles.hPoly) > 1)
			set(handles.check_allPolygs,'Vis','on')
		end
	end

	if (~handles.have_nans)			% If no NaNs than better to not even show these two
		set([handles.check_trinary handles.edit_NaN],'Vis','off')
	end

	% Fill the image size boxes
	[m, n, k] = size(get(handles.hMirImg,'CData'));
	set(handles.edit_nRows, 'String', m)
	set(handles.edit_nCols, 'String', n)
	handles.nRows = m;
	handles.nCols = n;

	set(hObject,'Vis','on');

	% Add this figure handle to the carraças list
	plugedWin = getappdata(varargin{1}.figure1,'dependentFigs');
	plugedWin = [plugedWin hObject];
	setappdata(varargin{1}.figure1,'dependentFigs',plugedWin);

	guidata(hObject, handles);
	if (nargout),	varargout{1} = 	hObject;	end

% -------------------------------------------------------------------------
function radio_binary_CB(hObject, handles)
	if (~get(hObject,'Value')),		set(hObject,'Value',1),		return,		end
	set(handles.radio_float,'Val',0)
	set([handles.edit_in handles.edit_out],'Enable','off')
	set([handles.radio_in handles.radio_out],'Enable','on')
	if (strcmp(get(handles.check_trinary,'Vis'), 'on'))	% They are visible when is grid and has NaNs
		set([handles.check_trinary handles.edit_NaN],'Enable','off')
	end

% -------------------------------------------------------------------------
function radio_float_CB(hObject, handles)
	if (~get(hObject,'Value')),		set(hObject,'Value',1),		return,		end
	set(handles.radio_binary,'Val',0)
	set([handles.edit_in handles.edit_out],'Enable','on')
	set([handles.radio_in handles.radio_out],'Enable','off')
	if (strcmp(get(handles.check_trinary,'Vis'), 'on'))	% They are visible when is grid and has NaNs
		set([handles.check_trinary handles.edit_NaN],'Enable','on')
	end


% -------------------------------------------------------------------------
function radio_in_CB(hObject, handles)
	if (~get(hObject,'Value')),		set(hObject,'Value',1),		return,		end
	set(handles.radio_out,'Val',0)

% -------------------------------------------------------------------------
function radio_out_CB(hObject, handles)
	if (~get(hObject,'Value')),		set(hObject,'Value',1),		return,		end
	set(handles.radio_in,'Val',0)

% -------------------------------------------------------------------------
function edit_in_CB(hObject, handles)
	handles.grid_in = str2double(get(hObject,'String'));
	guidata(handles.figure1, handles)

% -------------------------------------------------------------------------
function edit_out_CB(hObject, handles)
	handles.grid_out = str2double(get(hObject,'String'));
	guidata(handles.figure1, handles)

% -------------------------------------------------------------------------
function edit_NaN_CB(hObject, handles)
	handles.grid_NaN = str2double(get(hObject,'String'));
	guidata(handles.figure1, handles)

% -------------------------------------------------------------------------
function edit_nRows_CB(hObject, handles)
	xx = str2double(get(hObject, 'String'));
	if (isnan(xx) || xx < 10)
		set(hObject, 'String', handles.nRows)
	end

% -------------------------------------------------------------------------
function edit_nCols_CB(hObject, handles)
	xx = str2double(get(hObject, 'String'));
	if (isnan(xx) || xx < 10)
		set(hObject, 'String', handles.nCols)
	end

% -------------------------------------------------------------------------
function push_ok_CB(hObject, handles)

	% See if a size request different from original image
	nRows = fix(str2double(get(handles.edit_nRows, 'String')));
	nCols = fix(str2double(get(handles.edit_nCols, 'String')));
	if ( (nRows == handles.nRows) && (nCols == handles.nCols) )
		nRows = 0;		nCols = 0;
	end

	if (~handles.inputIsLine)

		if (get(handles.check_allPolygs,'Val'))		% We have more than one closed poly
			hLine = handles.hPoly;
			hLine = hLine(ishandle(hLine));			% Some might have been killed meanwhile
			handles.hPoly = hLine;
			guidata(handles.figure1, handles);
		else
			hLine = handles.hPoly_current;
		end

		if (nRows == 0)					% No different size request
			I = get(handles.hMirImg,'CData');
		else
			I = false([nRows nCols]);	% ...
		end
		limits = getappdata(handles.hMirAxes,'ThisImageLims');
		x = get(hLine(1),'XData');   y = get(hLine(1),'YData');
		mask = img_fun('roipoly_j',limits(1:2),limits(3:4),I,x,y);
		mask2 =  [];
		for (k = 2:numel(hLine))
			x = get(hLine(k),'XData');   y = get(hLine(k),'YData');
			mask2 = img_fun('roipoly_j',limits(1:2),limits(3:4),I,x,y);
			mask = mask | mask2;		% Combine all masks
		end
		if (~isempty(mask2)),	clear mask2;	end
	else

		set(handles.hMirImg, 'Vis', 'off')
		if (nRows == 0)					% No different size request
			mask = imcapture(handles.hMirAxes, 'img', 0);
		else
			mask = imcapture(handles.hMirAxes, 'img', [nRows nCols]);
		end
		set(handles.hMirImg, 'Vis', 'on')
		mask = mask(:,:,1);
		mask = (mask ~= 255);
		if (strcmp(get(handles.hMirAxes, 'YDir'), 'normal'))
			mask = flipud(mask);
		end
	end

	if (get(handles.radio_float,'Val'))			% Float mask
		Z = single(false(size(mask)));
		Z(mask)  = single(handles.grid_in);
		Z(~mask) = single(handles.grid_out);
		if (get(handles.check_trinary,'Val'))
			[X,Y,Zg] = load_grd(guidata(handles.hMirFig));
			mask = isnan(Zg);
			Z(mask) = single( str2double(get(handles.edit_NaN,'Str')) );
		end
	else
		if (get(handles.radio_out,'Val')),		mask = ~mask;	end
		Z = mask;
	end

	if (handles.image_type ~= 2)
		tmp.head = handles.head;	tmp.name = 'Mask image';
		tmp.X = tmp.head(1:2);		tmp.Y = tmp.head(3:4);
		if (get(handles.radio_float,'Val'))
			tmp.name = 'Mask grid';
			tmp.X = linspace(tmp.X(1), tmp.X(2), size(Z,2));
			tmp.Y = linspace(tmp.Y(1), tmp.Y(2), size(Z,1));
		end
     	mirone(Z, tmp)
	else
		h = mirone(Z);
		set(h,'Name','Mask image')
	end

% -------------------------------------------------------------------------
function hPoly = find_closed(handles)
% Get handles of closed lines (lines & patches)

	hPoly = findobj(handles.hMirAxes,'Type','line');
	% Find if any of the eventual above line is closed
	if (~isempty(hPoly))
		vec = false(numel(hPoly),1);
		for (k = 1:numel(hPoly))
			x = get(hPoly(k),'XData');   y = get(hPoly(k),'YData');
			if (numel(x) >= 3 && x(1) == x(end) && y(1) == y(end) )
				vec(k) = true;
			end
		end
		hPoly = hPoly(vec);
	end
	
	hPoly = [hPoly; findobj(handles.hMirAxes,'Type','patch')];
	hPoly = unique([handles.hPoly_current; hPoly]);		% one of them is repeated

% --- the GUI figure. 
function poly2mask_fig_LayoutFcn(h1)

set(h1, 'Position',[520 660 301 130],...
'Color', get(0,'factoryUicontrolBackgroundColor'),...
'MenuBar','none',...
'Name','poly2mask',...
'NumberTitle','off',...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1');

uicontrol('Parent',h1, 'Position',[72 82 51 16],...
'Call',@poly2mask_fig_uiCB,...
'FontName','Helvetica',...
'String','Inside',...
'Style','radiobutton',...
'TooltipString','Polygon interior set to 1',...
'Value',1,...
'Tag','radio_in');

uicontrol('Parent',h1, 'Position',[5 82 60 16],...
'Call',@poly2mask_fig_uiCB,...
'FontName','Helvetica',...
'String','Outside',...
'Style','radiobutton',...
'TooltipString','Outside polygon set to 0',...
'Tag','radio_out');

uicontrol('Parent',h1, 'Position',[5 110 85 15],...
'Call',@poly2mask_fig_uiCB,...
'FontName','Helvetica',...
'String','Binary Mask',...
'Style','radiobutton',...
'TooltipString','Create a black and white mask image',...
'Value',1,...
'Tag','radio_binary');

uicontrol('Parent',h1, 'Position',[222 110 75 15],...
'Call',@poly2mask_fig_uiCB,...
'FontName','Helvetica',...
'String','Float Mask',...
'Style','radiobutton',...
'TooltipString','Create a float mask image',...
'Tag','radio_float');

uicontrol('Parent',h1, 'Position',[245 79 51 21],...
'BackgroundColor',[1 1 1],...
'Call',@poly2mask_fig_uiCB,...
'Enable','off',...
'String','1',...
'Style','edit',...
'TooltipString','Inside polygon value. (NaNs are alowed)',...
'Tag','edit_in');

uicontrol('Parent',h1, 'Position',[195 79 51 21],...
'BackgroundColor',[1 1 1],...
'Call',@poly2mask_fig_uiCB,...
'Enable','off',...
'String','0',...
'Style','edit',...
'TooltipString','Outside polygon value. (NaNs are alowed)',...
'Tag','edit_out');

uicontrol('Parent',h1, 'Position',[14 17 131 23],...
'Call',@poly2mask_fig_uiCB,...
'FontName','Helvetica',...
'String','Apply to all polygons',...
'Style','checkbox',...
'TooltipString','Apply settings to all polygons. If unchecked apply only to the selected one.',...
'Tag','check_allPolygs',...
'Visible','off');

uicontrol('Parent',h1, 'Position',[235 9 61 21],...
'Call',@poly2mask_fig_uiCB,...
'FontName','Helvetica',...
'FontSize',10,...
'FontWeight','bold',...
'String','OK',...
'Tag','push_ok');

uicontrol('Parent',h1, 'Position',[55 48 61 22],...
'BackgroundColor',[1 1 1],...
'Call',@poly2mask_fig_uiCB,...
'Style','edit',...
'TooltipString','Number of rows',...
'Tag','edit_nRows');

uicontrol('Parent',h1, 'Position',[235 48 61 22],...
'BackgroundColor',[1 1 1],...
'Call',@poly2mask_fig_uiCB,...
'Style','edit',...
'TooltipString','Number of columns',...
'Tag','edit_nCols');

uicontrol('Parent',h1, 'Position',[6 52 47 14],...
'HorizontalAlignment','right',...
'String','N rows',...
'Style','text',...
'Tag','text1');

uicontrol('Parent',h1, 'Position',[189 52 45 14],...
'HorizontalAlignment','right',...
'String','N cols',...
'Style','text',...
'Tag','text2');

uicontrol('Parent',h1, 'Position',[144 105 65 23],...
'Enable','off',...
'String','Trinary?',...
'Style','checkbox',...
'TooltipString','Set NaNs into a separate class',...
'Tag','check_trinary');

uicontrol('Parent',h1, 'Position',[145 79 51 21],...
'BackgroundColor',[1 1 1],...
'Call',@poly2mask_fig_uiCB,...
'Enable','off',...
'String','NaN',...
'Style','edit',...
'TooltipString','Replace Nans by this . (NaNs or ther values are alowed)',...
'Tag','edit_NaN');

function poly2mask_fig_uiCB(hObject, evt)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));

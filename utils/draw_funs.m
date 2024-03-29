function varargout = draw_funs(hand, varargin)
%function OUT = draw_funs(hand,opt,data)
%   This contains several functions necessary to the "Draw" menu of mirone
%   There are no error checking.
%   HAND    contains the handle to the graphical object (or eventually empty for direct access)
%   OPT     is a string for choosing what action to perform
%   DATA    contains data currently used in the volcanoes, fogspots and some other options
%
%	To use the direct access mode to any of the local functions that don't need to fish
%	the data from an object handle, call with HAND = []. E.g (in load_xyz)
%	draw_funs([], 'doSave_formated', x, y, z)

%	Copyright (c) 2004-2020 by J. Luis
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

% $Id: draw_funs.m 11441 2019-07-24 14:35:47Z j $

% A bit of strange tests but they are necessary for the cases when we use the new feval(fun,varargin{:}) 
opt = varargin{1};		% function name to evaluate (new) or keyword to select one (old form)
if (numel(varargin) > 1)
	data = varargin{2};
end

switch opt
	case 'line_uicontext',			set_line_uicontext(hand,'line')
	case 'setSHPuictx',				setSHPuictx(hand)
	case 'ContourLines',			set_ContourLines_uicontext(hand,data)
	case 'MBtrackUictx',			set_line_uicontext(hand,'MBtrack')
	case 'MBbarUictx',				set_bar_uicontext(hand)
	case 'CoastLineUictx',			setCoastLineUictx(hand)
	case 'deleteObj',				deleteObj(hand);
	case 'DrawCircleEulerPole',		draw_circleEulerPole(data(1),data(2));  
	case 'SessionRestoreCircle'			% Called by "FileOpenSession" or "DrawGeogCircle_CB"
		set_circleGeo_uicontext(hand)
	case 'SessionRestoreCircleCart'		% Called by "FileOpenSession" or "DrawGeogCircle_CB"
		set_circleCart_uicontext(hand)
	case 'DrawText'
		cmenuHand = uicontextmenu('parent',get(get(hand,'parent'),'parent'));		% We know 'hand' is a text handle
		set(hand, 'UIContextMenu', cmenuHand);
		cb_color = uictx_color(hand);	% there are 9 cb_color outputs
		uimenu(cmenuHand, 'Label', 'Change Font', 'Callback', @text_FontSize);
		item_fc = uimenu(cmenuHand, 'Label', 'Font Color');
		setLineColor(item_fc,cb_color)
		uimenu(cmenuHand, 'Label', 'Edit   text', 'Callback', 'set(gco, ''Editing'', ''on''); refresh', 'Sep','on');
		uimenu(cmenuHand, 'Label', 'Copy   text', 'Callback', @copy_text_object);
		uimenu(cmenuHand, 'Label', 'Delete text', 'Callback', 'delete(gco); refresh');
		uimenu(cmenuHand, 'Label', 'Move   text', 'Callback', @move_text);
		uimenu(cmenuHand, 'Label', 'Rotate text', 'Callback', @rotate_text);
		uimenu(cmenuHand, 'Label', 'Export text', 'Callback', @export_text);
		set(hand, 'ButtonDownFcn', 'move_obj(1)')
	case 'DrawSymbol',			set_symbol_uicontext(hand)
	case {'hotspot','volcano','ODP','City_major','City_other','Earthquakes','TideStation', 'Meteor', 'Hydro'}
		set_symbol_uicontext(hand,data)
	case 'PlateBoundPB',		set_PB_uicontext(hand,data)
	case 'ChngAxLabels',		changeAxesLabels(data)
	case 'SRTMrect',			set_SRTM_rect_uicontext(hand)
	case 'isochron',			set_isochrons_uicontext(hand,data)
	case 'gmtfile',				set_gmtfile_uicontext(hand,data)
	case 'country_patch',		set_country_uicontext(hand)
	case 'telhas_patch',		set_telhas_uicontext(hand)
	case 'save_xyz',			save_formated([],[],[], data)
	case 'tellAzim',			show_lineAzims([],[], hand);
	case 'tellLLength',			show_LineLength([],[], hand);
	case 'tellArea',			show_Area([],[], hand);
	otherwise
		if (nargout)
			[varargout{1:nargout}] = ...
				feval(opt, varargin{2:end});	% NEW. Eventualy, all calls should evolve to use this form
		else
			feval(opt, varargin{2:end});
		end
end
% Now short-cuted:
% 'DrawVector', 'magbarcode' 'Ctrl_v' 'DrawGreatCircle' 'DrawCartesianCircle' 'loc_quiver'

% -----------------------------------------------------------------------------------------
function Ctrl_v(h)
% Paste a line whose handle is h(1) in figure gcf that is different from parent(h(1))
% If H has two elements, the second should contain the CurrentAxes
	hLine = h(1);			% No testing. Do not fail
	if (numel(h) == 2),		hAx = h(2);
	else,					hAx = get(get(0,'CurrentFigure'), 'CurrentAxes');
	end
	x = get(hLine, 'xdata');	y = get(hLine, 'ydata');
	if (strcmp(get(hLine,'type'), 'line'))
		h = line('xdata',x, 'ydata',y, 'Parent', hAx, 'LineWidth', get(hLine,'LineWidth'), ...
			'LineStyle',get(hLine,'LineStyle'), 'Color',get(hLine,'Color'), 'Tag',get(hLine,'Tag') );
		marker = get(hLine, 'Marker');				markSize = get(hLine,'MarkerSize');
		markFC = get(hLine,'MarkerFaceColor');		markEC = get(hLine,'MarkerEdgeColor');
		% Set the markers only if they are different from the red square markers used for edition
		% However, this is a risky test because if one of them is ever changed the test will fail.
		if (~(strcmp(marker, 'square') && strcmp(markFC, 'none') && markSize == 6 && isequal(markEC,[1 0 0]) ))
			set(h, 'Marker', marker, 'MarkerSize',markSize, 'MarkerEdgeColor',markEC, 'MarkerFaceColor',markFC)
		end
	else
		h = patch('xdata',x, 'ydata',y, 'Parent', hAx, 'LineWidth', get(hLine,'LineWidth'), ...
			'LineStyle',get(hLine,'LineStyle'), 'EdgeColor',get(hLine,'EdgeColor'), 'FaceColor',get(hLine,'FaceColor'), ...
			 'FaceAlpha',get(hLine,'FaceAlpha'), 'Tag',get(hLine,'Tag') );
	end
	z = getappdata(hLine,'ZData');
	if (~isempty(z)),		setappdata(h, 'ZData', z);		end
	set_line_uicontext(h,'line')		% Set lines's uicontextmenu

% % -----------------------------------------------------------------------------------------
% function setUIcbs(item, labels, cbs)
% % Set uimenu uicontexts of graphic elements
% 	for (k = 1:numel(cbs))
% 		uimenu(item, 'Label', labels{k}, 'Callback', cbs{k});
% 	end

% -----------------------------------------------------------------------------------------
function setLineStyle(item,cbs)
% Set the line Style uicontexts of graphic elements
	uimenu(item, 'Label', 'solid', 'Callback', cbs{1});
	uimenu(item, 'Label', 'dashed', 'Callback', cbs{2});
	uimenu(item, 'Label', 'dotted', 'Callback', cbs{3});
	uimenu(item, 'Label', 'dash-dotted', 'Callback', cbs{4});
	if (numel(cbs) == 6)
		uimenu(item, 'Label', 'Show markers', 'Callback', cbs{5}, 'Sep', 'on');
		uimenu(item, 'Label', 'Hide markers', 'Callback', cbs{6});
	end

% -----------------------------------------------------------------------------------------
function setLineColor(item,cbs)
% Set the line color uicontexts of graphic elements
	uimenu(item, 'Label', 'Black', 'Callback', cbs{1});
	uimenu(item, 'Label', 'White', 'Callback', cbs{2});
	uimenu(item, 'Label', 'Red', 'Callback', cbs{3});
	uimenu(item, 'Label', 'Green', 'Callback', cbs{4});
	uimenu(item, 'Label', 'Blue', 'Callback', cbs{5});
	uimenu(item, 'Label', 'Yellow', 'Callback', cbs{6});
	uimenu(item, 'Label', 'Cyan', 'Callback', cbs{7});
	uimenu(item, 'Label', 'Magenta', 'Callback', cbs{8});
	uimenu(item, 'Label', 'Other...', 'Callback', cbs{9});

% -----------------------------------------------------------------------------------------
function setLineWidth(item,cbs)
% Set the line color uicontexts of graphic elements
	uimenu(item, 'Label', '1       pt', 'Callback', cbs{1});
	uimenu(item, 'Label', '2       pt', 'Callback', cbs{2});
	uimenu(item, 'Label', '3       pt', 'Callback', cbs{3});
	uimenu(item, 'Label', '4       pt', 'Callback', cbs{4});
	uimenu(item, 'Label', 'Other...', 'Callback', cbs{5});

% -----------------------------------------------------------------------------------------
function setSHPuictx(h,opt)
% h is a handle to a shape line object

	handles = guidata(h(1));
	for (i = 1:numel(h))
		cmenuHand = uicontextmenu('Parent',handles.figure1);      set(h(i), 'UIContextMenu', cmenuHand);
		uimenu(cmenuHand, 'Label', 'Save line', 'Callback', {@save_formated,h});
		uimenu(cmenuHand, 'Label', 'Delete this line', 'Callback', {@del_line,h(i)});
		uimenu(cmenuHand, 'Label', 'Delete class', 'Callback', 'delete(findobj(''Tag'',''SHPpolyline''))');

		if (~isempty(get(h(i), 'ZData')))				% If we have z info
			item = uimenu(cmenuHand, 'Label', 'Quick grid', 'Sep', 'on');
			uimenu(item, 'Label', 'auto', 'Callback', {@pt_quick_grd,h(i),'a'});
			uimenu(item, 'Label', 'set increments', 'Callback', {@pt_quick_grd, h(i), 'q'});
		end

		cb_solid   = 'set(gco, ''LineStyle'', ''-''); refresh';
		cb_dashed  = 'set(gco, ''LineStyle'', ''--''); refresh';
		cb_dotted  = 'set(gco, ''LineStyle'', '':''); refresh';
		cb_dashdot = 'set(gco, ''LineStyle'', ''-.''); refresh';

		item = uimenu(cmenuHand, 'Label', 'Line Width', 'Sep','on');
		uimenu(item, 'Label', 'Other...', 'Callback', {@other_LineWidth,h(i)});

		item = uimenu(cmenuHand, 'Label', 'Line Style');
		setLineStyle(item,{cb_solid cb_dashed cb_dotted cb_dashdot})

		item = uimenu(cmenuHand, 'Label', 'Line Color');
		uimenu(item, 'Label', 'Other...', 'Callback', {@other_color,h(i)});
		ui_edit_polygon(h(i))
		isPt = getappdata(h(1), 'isPoint');
		if (~isempty(isPt) && ~isPt)	% For points it makes no sense a 'Join lines'
			uimenu(cmenuHand, 'Label', 'Join lines', 'Callback', {@join_lines,handles.figure1});
		end
	end

% -----------------------------------------------------------------------------------------
function set_rect_uictx_PC(h, opt)
% Set the line uicontext for rectangles plotted in Plot Composer figure.
	if (isempty(h)),	return,		end
	hFig = getParentFigure(h(1));
	handles = guidata(hFig);		% Get PC handles

	for (k = 1:numel(h))
		cmenuHand = uicontextmenu('Parent',handles.figure1);
		set(h(k), 'UIContextMenu', cmenuHand);
		set(cmenuHand, 'UserData', h(k))			% And with this the cmenuHand knows to whom it belongs

		cb_LineWidth = uictx_LineWidth(h(k));		% there are 5 cb_LineWidth outputs
		cb_solid   = 'set(gco, ''LineStyle'', ''-''); refresh';
		cb_dashed  = 'set(gco, ''LineStyle'', ''--''); refresh';
		cb_dotted  = 'set(gco, ''LineStyle'', '':''); refresh';
		cb_dashdot = 'set(gco, ''LineStyle'', ''-.''); refresh';
		cb_color   = uictx_color(h(k));				% there are 9 cb_color outputs

		uimenu(cmenuHand,  'Label', 'Delete', 'Callback', {@del_line,h(k)});
		item_tools = uimenu(cmenuHand, 'Label', 'Trim with rect');
		uimenu(item_tools, 'Label', 'inside', 'Callback', {@trim_withPolygon,h,1});
		uimenu(item_tools, 'Label', 'outside', 'Callback', {@trim_withPolygon,h,0});
		uimenu(cmenuHand,  'Label', 'Copy', 'Callback', {@copy_line_object,handles.figure1,handles.axes1});
		uimenu(cmenuHand,  'Label', 'Rectangle limits (edit)', 'Sep','on', 'Callback', @rectangle_limits);
		CB = 'deal_opts([], ''assoc_gmt_symbol'', [], [])';
		uimenu(cmenuHand,  'Label', 'Insert a PostScript file here','Callback', CB);
		setLineWidth(uimenu(cmenuHand, 'Label', 'Line Width', 'Sep','on'), cb_LineWidth)
		setLineStyle(uimenu(cmenuHand, 'Label', 'Line Style'), {cb_solid cb_dashed cb_dotted cb_dashdot})
		setLineColor(uimenu(cmenuHand, 'Label', 'Line Color'), cb_color)
		set_stack_order(cmenuHand)      % Set options to change order in the stackpot
		ui_edit_polygon(h(k))			% Set edition functions
	end

% -----------------------------------------------------------------------------------------
function set_line_uicontext_XY(h, opt)
% Set the line uicontext for lines plotted in a profiler (Ecran) figure. They are much
% simpler that the Mirone ones, so put it in a separate function (this function).
	if (isempty(h)),	return,		end
	if (nargin == 1),	opt = '';	end
	hFig = getParentFigure(h(1));
	handles = guidata(hFig);		% Get Ecran handles

	for (k = 1:numel(h))
		cmenuHand = uicontextmenu('Parent',handles.figure1);
		set(h(k), 'UIContextMenu', cmenuHand);
		set(cmenuHand, 'UserData', h(k))			% And with this the cmenuHand knows to whom it belongs

		cb_LineWidth = uictx_LineWidth(h(k));		% there are 5 cb_LineWidth outputs
		cb_solid   = 'set(gco, ''LineStyle'', ''-''); refresh';
		cb_dashed  = 'set(gco, ''LineStyle'', ''--''); refresh';
		cb_dotted  = 'set(gco, ''LineStyle'', '':''); refresh';
		cb_dashdot = 'set(gco, ''LineStyle'', ''-.''); refresh';
		cb_markers_on = 'set(gco, ''Marker'', ''o'', ''MarkerFaceColor'', get(gco, ''Color''), ''MarkerSize'',4); refresh';
		cb_markers_off = 'set(gco, ''Marker'', ''none''); refresh';
		cb_color   = uictx_color(h(k));				% there are 9 cb_color outputs

		uimenu(cmenuHand, 'Label', 'Delete', 'Callback', {@del_line,h(k)});
		uimenu(cmenuHand, 'Label', 'Save line', 'Callback', {@save_formated,h(k)});
		uimenu(cmenuHand, 'Label', 'Copy line', 'Callback', 'setappdata(0,''CtrlCHandEcran'',gco)');
		if (strcmp(opt, 'main'))
			% Attention, if I ever change these labels I MUST do it also in ecran/finish_line_uictx()
			hh = uimenu(cmenuHand, 'Label', 'Shift origin here');
			uimenu('Parent',hh, 'Label', 'X origin only');
			uimenu('Parent',hh, 'Label', 'Y origin only');
			uimenu('Parent',hh, 'Label', 'XY origin');
			uimenu(cmenuHand,   'Label', 'Filter Outliers', 'Sep','on');
			uimenu(cmenuHand,   'Label', 'Filter line');
			uimenu(cmenuHand,   'Label', 'Show histogram',  'Sep','on');
			uimenu(cmenuHand,   'Label', 'Show Bar graph');
		else
			ui_edit_polygon(h(k))			% Set edition functions
			uimenu(cmenuHand, 'Label', 'Line length', 'Callback', @show_LineLength_XY)
		end
		setLineWidth(uimenu(cmenuHand, 'Label', 'Line Width', 'Sep','on'), cb_LineWidth)
		setLineStyle(uimenu(cmenuHand, 'Label', 'Line Style'), {cb_solid cb_dashed cb_dotted cb_dashdot cb_markers_on cb_markers_off})
		setLineColor(uimenu(cmenuHand, 'Label', 'Line Color'), cb_color)
	end

% -----------------------------------------------------------------------------------------
function set_line_uicontext(h, opt)
% h is a handle to a line object (that can be closed)
	if (isempty(h)),	return,		end
	
	if (numel(h) > 1)			% Many, make recursive and hope this doesn't choke the fragile beast
		for (k = 1:numel(h))
			set_line_uicontext(h(k), opt)
		end
		return
	end

	IS_SEISPOLYG = false;		% Seismicity Polygons have special options
	IS_SEISMICLINE = false;		% Seismicity Lines have special options
	LINE_ISCLOSED = false;		IS_RECTANGLE = false;	IS_PATCH = false;
	IS_ARROW = false;
	% Check to see if we are dealing with a closed polyline
	x = get(h,'XData');			y = get(h,'YData');
	if (isempty(x) || isempty(y)),		return,		end		% Line is totally out of the figure
	if ((x(1) == x(end)) && (y(1) == y(end)))
		LINE_ISCLOSED = true;
		if (length(x) == 5 && (x(1) == x(2)) && (x(3) == x(4)) && (y(1) == y(4)) && (y(2) == y(3)))
			IS_RECTANGLE = true;	
		end  
		if (strcmp(get(h,'Tag'),'SeismicPolyg')),	IS_SEISPOLYG = true;	end
	end
	if (strcmp(get(h,'Type'),'patch'))
		IS_PATCH = true;
		if (IS_PATCH && ~LINE_ISCLOSED),	LINE_ISCLOSED = true;	end
	elseif (strcmp(get(h,'Tag'),'SeismicLine'))
		IS_SEISMICLINE = true;
	end

	handles = guidata(get(h,'Parent'));		% Get Mirone handles

	IamSpectrum = false;
	if (handles.validGrid && ~isempty(strfind(get(handles.figure1,'Name'), 'spectrum')))
		IamSpectrum = true;
	end

	% Check to see if we are dealing with a multibeam track
	cmenuHand = uicontextmenu('Parent',handles.figure1);
	set(h, 'UIContextMenu', cmenuHand);
	set(cmenuHand, 'UserData', h)			% And with this the cmenuHand knows to whom it belongs
	switch opt
		case 'line'
			label_save = 'Save line';   label_length = 'Line length(s)';   label_azim = 'Line azimuth(s)';
			IS_LINE = true;		IS_MBTRACK = false;
			if (strcmp(get(h,'Tag'),'Seta')),		IS_ARROW = true;	end
		case 'MBtrack'
			label_save = 'Save track';   label_length = 'Track length';   label_azim = 'Track azimuth(s)';
			IS_LINE = false;	IS_MBTRACK = true;
	end
	cb_LineWidth = uictx_LineWidth(h);		% there are 5 cb_LineWidth outputs
	cb_solid = 'set(gco, ''LineStyle'', ''-''); refresh';
	cb_dashed = 'set(gco, ''LineStyle'', ''--''); refresh';
	cb_dotted = 'set(gco, ''LineStyle'', '':''); refresh';
	cb_dashdot = 'set(gco, ''LineStyle'', ''-.''); refresh';
	cb_color = uictx_color(h);				% there are 9 cb_color outputs

	if (IS_RECTANGLE)
		uimenu(cmenuHand, 'Label', 'Delete me', 'Callback', {@del_line,h});
		uimenu(cmenuHand, 'Label', 'Delete inside rect', 'Callback', {@del_insideRect,h});
		item_tools = uimenu(cmenuHand, 'Label', 'Trim with rect');
		uimenu(item_tools, 'Label', 'inside', 'Callback', {@trim_withPolygon,h,1});
		uimenu(item_tools, 'Label', 'outside', 'Callback', {@trim_withPolygon,h,0});
		ui_edit_polygon(h)
	elseif (IS_LINE)
		uimenu(cmenuHand, 'Label', 'Delete', 'Callback', {@del_line,h});
		ui_edit_polygon(h)			% Set edition functions
	elseif (IS_MBTRACK)				% Multibeam tracks, when deleted, have to delete also the bars
		uimenu(cmenuHand, 'Label', 'Delete track (left-click on it)', 'Callback', 'save_track_mb(1);');
		% Old style edit function. New edit is provided by ui_edit_polygon which doesn't work with mbtracks 
		uimenu(cmenuHand, 'Label', 'Edit track (left-click on it)', 'Callback', 'edit_track_mb');
	end
	uimenu(cmenuHand, 'Label', label_save, 'Callback', {@save_formated,h});
	if (~IS_SEISPOLYG && ~IS_MBTRACK && ~strcmp(get(h,'Tag'),'FaultTrace'))	% Those are not to allowed to copy
		if (~LINE_ISCLOSED && ~IS_ARROW)
			uimenu(cmenuHand, 'Label', 'Join lines', 'Callback', {@join_lines,handles.figure1});
		end
		uimenu(cmenuHand, 'Label', 'Copy', 'Callback', {@copy_line_object,handles.figure1,handles.axes1});
	end
	if (~IS_SEISPOLYG && ~IS_ARROW && ~IS_RECTANGLE)
		if (numel(x) > 2)
			uimenu(cmenuHand, 'Label', 'Spline Smooth', 'Callback', {@smooth_line,h})
		end
		uimenu(cmenuHand, 'Label', label_length, 'Callback', @show_LineLength)
		uimenu(cmenuHand, 'Label', label_azim,   'Callback', @show_lineAzims)
	end
	if (IS_MBTRACK),	uimenu(cmenuHand, 'Label', 'All MB-tracks length', 'Callback', @show_AllTrackLength);	end

	if (LINE_ISCLOSED)
		uimenu(cmenuHand, 'Label', 'Area under polygon', 'Callback', @show_Area);
		if (IS_PATCH && ~IS_SEISPOLYG)
			item8 = uimenu(cmenuHand, 'Label','Fill Color');
			setLineColor(item8, uictx_color(h, 'facecolor'))		% there are 9 cb_color outputs
			uimenu(item8, 'Label', 'None', 'Sep','on', 'Callback', 'set(gco, ''FaceColor'', ''none'');refresh');
			uimenu(cmenuHand, 'Label', 'Transparency', 'Callback', @set_transparency);
		end
		if (~IamSpectrum)
			uimenu(cmenuHand, 'Label', 'Create Mask', 'Callback', 'poly2mask_fig(guidata(gcbo),gco)');
		end
	end
	if (IS_RECTANGLE && handles.validGrid && ~IamSpectrum)
		uimenu(cmenuHand, 'Label', 'Make Chess board', 'Callback', {@chessify, h});
	end

	if (handles.image_type ~= 20 && ~LINE_ISCLOSED && strcmp(opt,'line'))
		%uimenu(cmenuHand, 'Label', 'Stacked profile', 'Callback', @stack_profiles);
		if ((ndims(get(handles.hImg,'CData')) == 2) || (handles.validGrid))	% Because Track of RGB doesn't know how to save
			if (handles.nLayers > 1)
				cbTrack = 'setappdata(gcf,''TrackThisLine'',gco); mirone(''ExtractProfile_CB'',guidata(gcbo),''3D'')';
				uimenu(cmenuHand, 'Label', '3D interpolation', 'Callback', cbTrack);
			end
			cbTrack = 'setappdata(gcf,''TrackThisLine'',gco); mirone(''ExtractProfile_CB'',guidata(gcbo),''point'')';
			uimenu(cmenuHand, 'Label', 'Point interpolation', 'Callback', cbTrack);
		end
		cbTrack = 'setappdata(gcf,''TrackThisLine'',gco); mirone(''ExtractProfile_CB'',guidata(gcbo))';
		uimenu(cmenuHand, 'Label', 'Extract profile', 'Callback', cbTrack);
	end

	if strcmp(opt,'MBtrack'),	uimenu(cmenuHand, 'Label', 'Show track''s Swath Ratio', 'Callback', {@show_swhatRatio,h});	end

	if (IS_RECTANGLE)
		uimenu(cmenuHand, 'Label', 'Rectangle limits (edit)', 'Sep','on', 'Callback', @rectangle_limits);
		if (handles.image_type ~= 20 && ~IamSpectrum), uimenu(cmenuHand, 'Label', 'Register Image', 'Callback', @rectangle_register_img);	end
		if (~handles.validGrid)
			uimenu(cmenuHand, 'Label', 'Transplant Image here','Callback', 'transplants(gco,''image'')');
		end
		if (handles.geog)
			uimenu(cmenuHand, 'Label', 'Grid/Image mosaicer', 'Callback', 'mosaicer(gco)');
			%uimenu(cmenuHand, 'Label', 'Get image from Web Map Server', 'Callback', 'wms_tool(gco)');
			uimenu(cmenuHand, 'Label', 'CMT Catalog (Web download)', 'Callback', 'globalcmt(gcf,gco)');
		end

		item_tools = uimenu(cmenuHand, 'Label', 'ROI Crop Tools','Sep','on');
		if (handles.validGrid)    % Option only available to recognized grids
			uimenu(item_tools, 'Label', 'Crop Grid', 'Callback', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''CropaGrid_pure'')');
			uimenu(item_tools, 'Label', 'Crop Image', 'Callback', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco)');
			uimenu(item_tools, 'Label', 'Crop Image (with coords)', 'Callback', ...
				'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''CropaWithCoords'')');
% 			BL = getappdata(handles.figure1,'BandList');
% 			if (~isempty(BL))
% 				uimenu(item_tools, 'Label', 'Crop stack', 'Callback', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''CropaStack'')');
% 			end
			uimenu(item_tools, 'Label', 'Set to constant', 'Callback', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''SetConst'')','Sep','on');
			item_tools3 = uimenu(item_tools, 'Label', 'Stats');
			uimenu(item_tools3,'Label', 'Mean',          'Callback', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''GetMean'')');
			uimenu(item_tools3,'Label', 'Median',        'Callback', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''GetMedian'')');
			uimenu(item_tools3,'Label', 'STD',           'Callback', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''GetSTD'')');
			uimenu(item_tools, 'Label', 'Clip grid',     'Callback', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''Clip'')');
			uimenu(item_tools, 'Label', 'Median filter', 'Callback', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''MedianFilter'')');
			uimenu(item_tools, 'Label', 'Spline smooth', 'Callback', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''SplineSmooth'')');
			uimenu(item_tools, 'Label', 'Histogram (grid)', 'Callback', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''CropaGrid_histo'')');
			uimenu(item_tools, 'Label', 'Detect Fronts', 'Callback', 'cayula_cornillon(guidata(gcbo),gco)');
			uimenu(item_tools, 'Label', 'Power',           'Callback', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''CropaGrid_power'')');
			uimenu(item_tools, 'Label', 'Autocorrelation', 'Callback', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''CropaGrid_autocorr'')');
			uimenu(item_tools, 'Label', 'FFT tool', 'Callback', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''CropaGrid_fftTools'')');
			item_fill = uimenu(item_tools, 'Label', 'Fill gaps');
			uimenu(item_fill,  'Label', 'Fill gaps (surface)','Callback', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''FillGaps'',''surface'')');
			uimenu(item_fill,  'Label', 'Fill gaps (cubic)',  'Callback', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''FillGaps'',''cubic'');');
			uimenu(item_fill,  'Label', 'Fill gaps (linear)', 'Callback', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''FillGaps'',''linear'');');
			uimenu(item_tools, 'Label', 'Fill sinks', 'Callback', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''FillSinks_pitt'')');
			uimenu(item_tools, 'Label', 'Slice peaks','Callback', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''FillSinks_peak'')');
		elseif (handles.image_type ~= 20)			% We have an Image
			uimenu(item_tools, 'Label', 'Crop Image', 'Callback', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco)');
			if (handles.image_type == 3)
					uimenu(item_tools, 'Label', 'Crop Image (with coords)', 'Callback', ...
						'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''CropaWithCoords'')');
			end
			uimenu(item_tools, 'Label', 'Histogram', 'Callback', 'image_histo(guidata(gcbo),gco)');
		end
		if (handles.validGrid || handles.image_type ~= 20)
			uimenu(item_tools, 'Label', 'Image process', 'Sep','on', 'Callback', 'mirone(''DrawClosedPolygon_CB'',guidata(gcbo),gco)');
		else
			delete(item_tools)		% Easier to create it at top and delete in the case where it had no use
		end
		if (~IamSpectrum)
			deal_opts({'MGG' 'MICROLEV' 'GMT_DB_IDS' 'GMT_SYMBOL' 'GRAVITY'}, cmenuHand);
		end
	end

	if (~IS_SEISPOLYG && LINE_ISCLOSED && ~IS_RECTANGLE && (handles.validGrid || handles.image_type ~= 20))
		item_tools2 = set_ROI_uicontext(handles, cmenuHand, IamSpectrum);
		uimenu(item_tools2, 'Label', 'Image process', 'Sep','on', 'Callback', 'mirone(''DrawClosedPolygon_CB'',guidata(gcbo),gco)');

		if (strcmp(get(h,'Tag'),'EulerTrapezium'))
			uimenu(cmenuHand, 'Label', 'Compute Euler Pole', 'Sep','on', 'Callback',...
				'calc_bonin_euler_pole(get(gco,''XData''), get(gco,''YData''));' );
		end
	elseif (IS_RECTANGLE && IamSpectrum)
		uimenu(cmenuHand, 'Label', 'Low Pass FFT filter', 'Callback', 'mirone(''GridToolsSectrum_CB'',guidata(gcbo), ''lpass'', gco)', 'Sep','on');
		uimenu(cmenuHand, 'Label', 'High Pass FFT filter','Callback', 'mirone(''GridToolsSectrum_CB'',guidata(gcbo), ''hpass'', gco)');
	end

	setLineWidth(uimenu(cmenuHand, 'Label', 'Line Width', 'Sep','on'), cb_LineWidth)
	setLineStyle(uimenu(cmenuHand, 'Label', 'Line Style'), {cb_solid cb_dashed cb_dotted cb_dashdot})
	if (IS_PATCH),		cb_color = uictx_color(h,'EdgeColor');	end      % there are 9 cb_color outputs
	setLineColor(uimenu(cmenuHand, 'Label', 'Line Color'), cb_color)
	set_stack_order(cmenuHand)      % Set options to change order in the stackpot

	if (strcmp(get(h,'Tag'),'FaultTrace'))      % For Okada modeling
		uimenu(cmenuHand, 'Label', 'Okada', 'Sep','on', 'Callback', {@okada_model,h,'okada'});    
		uimenu(cmenuHand, 'Label', 'Mansinha', 'Callback', {@okada_model,h,'mansinha'});    
	end

	if (IS_SEISPOLYG)                         % Seismicity options
		% gco gives the same handle as h 
		uimenu(cmenuHand, 'Label', 'Save events', 'Callback', 'save_seismicity(gcf,[],gco)', 'Sep','on');
		uimenu(cmenuHand, 'Label', 'Find clusters', 'Callback', 'find_clusters(gcf,gco)');
		itemHist = uimenu(cmenuHand, 'Label','Histograms');
		uimenu(itemHist, 'Label', 'Guttenberg & Richter', 'Callback', 'histos_seis(gco,''GR'')');
		uimenu(itemHist, 'Label', 'Cumulative number', 'Callback', 'histos_seis(gco,''CH'')');
		uimenu(itemHist, 'Label', 'Cumulative moment', 'Callback', 'histos_seis(gco,''CM'')');
		uimenu(itemHist, 'Label', 'Magnitude', 'Callback', 'histos_seis(gco,''MH'')');
		uimenu(itemHist, 'Label', 'Time', 'Callback', 'histos_seis(gco,''TH'')');
		uimenu(itemHist, 'Label', 'Display in Table', 'Callback', 'histos_seis(gcf,''HT'')','Sep','on');
		%uimenu(itemHist, 'Label', 'Hour of day', 'Callback', 'histos_seis(gco,''HH'')');
		itemTime = uimenu(cmenuHand, 'Label','Time series');
		uimenu(itemTime, 'Label', 'Time magnitude', 'Callback', 'histos_seis(gco,''TM'')');
		uimenu(itemTime, 'Label', 'Time depth', 'Callback', 'histos_seis(gco,''TD'')');
		uimenu(cmenuHand, 'Label', 'Mc and b estimate', 'Callback', 'histos_seis(gco,''BV'')');
		uimenu(cmenuHand, 'Label', 'Fit Omori law', 'Callback', 'histos_seis(gco,''OL'')');
		%uimenu(cmenuHand, 'Label', 'Skell', 'Callback', 'esqueleto_tmp(gco)','Sep','on');
	elseif (IS_SEISMICLINE)
		itemHist = uimenu(cmenuHand, 'Label','Set buffer zone', 'Sep','on');
		uimenu(itemHist, 'Label', '1 deg width', 'Callback', {@seismic_line,h,'def'});
		uimenu(itemHist, 'Label', 'other width', 'Callback', {@seismic_line,h,'buf'});
		uimenu(cmenuHand, 'Label', 'Project seismicity', 'Callback', {@seismic_line,h,'proj'}, 'Enable','off');
	end

% -----------------------------------------------------------------------------------------
function item_tools2 = set_ROI_uicontext(handles, cmenuHand, IamSpectrum)
% ...
	item_tools2 = uimenu(cmenuHand, 'Label', 'ROI Crop Tools','Sep','on');
	if (handles.validGrid)    % Option only available to recognized grids
		uimenu(item_tools2, 'Label', 'Crop Grid', 'Callback', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''CropaGrid_pure'')');
		uimenu(item_tools2, 'Label', 'Crop Image', 'Callback', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco)');
		uimenu(item_tools2, 'Label', 'Crop Image (with coords)', 'Callback', ...
			'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''CropaWithCoords'')');
		uimenu(item_tools2, 'Label', 'Set to constant', 'Callback', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''ROI_SetConst'')','Sep','on');
		item_tools3 = uimenu(item_tools2, 'Label', 'Stats');
		uimenu(item_tools3, 'Label', 'Mean', 'Callback', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''ROI_Mean'')');
		uimenu(item_tools3, 'Label', 'Median', 'Callback', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''ROI_Median'')');
		uimenu(item_tools3, 'Label', 'STD', 'Callback', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''ROI_STD'')');
		uimenu(item_tools2, 'Label', 'Clip grid', 'Callback', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''ROI_Clip'')');
		uimenu(item_tools2, 'Label', 'Median filter', 'Callback', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''ROI_MedianFilter'')');
		uimenu(item_tools2, 'Label', 'Spline smooth', 'Callback', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''ROI_SplineSmooth'')');
		uimenu(item_tools2, 'Label', 'Fill sinks', 'Callback', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''ROI_FillSinks_pitt'')');
		uimenu(item_tools2, 'Label', 'Slice peaks','Callback', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''ROI_FillSinks_peak'')');
		uimenu(item_tools2, 'Label', 'Histogram (grid)', 'Callback', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''CropaGrid_histo'')');
		hP = getappdata(handles.figure1, 'ParentFig');
		if (~isempty(hP) && ishandle(hP) && IamSpectrum)
			uimenu(item_tools2, 'Label', 'Low Pass FFT filter', 'Callback', 'mirone(''GridToolsSectrum_CB'',guidata(gcbo), ''lpass'', gco)');
			uimenu(item_tools2, 'Label', 'High Pass FFT filter','Callback', 'mirone(''GridToolsSectrum_CB'',guidata(gcbo), ''hpass'', gco)');
		end
	else			% We have an Image
		uimenu(item_tools2, 'Label', 'Crop Image', 'Callback', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco)');
		if (handles.image_type == 3)
			uimenu(item_tools2, 'Label', 'Crop Image (with coords)', 'Callback', ...
				'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''CropaWithCoords'')');
		end
		uimenu(item_tools2, 'Label', 'Histogram', 'Callback', 'image_histo(guidata(gcbo),gco)');
	end

% -----------------------------------------------------------------------------------------
function set_recTsu_uicontext(h)
% Set options particular to the NESTING nature of rectangles for TSUNAMI grid construction.

	handles = guidata(h(1));
	for (k = 1:numel(h))
		do_nest_sizes = false;
		cmenuHand = uicontextmenu('Parent',handles.figure1);
		set(h(k), 'UIContextMenu', cmenuHand)

		ud = get(h(k), 'Userdata');
		if (isempty(ud) || ud == 1)		% Root rectangle, set its Userdata to 1 to flag that fact (it's = 1 when restoring sess)
			if (handles.validGrid)		% Give this rectangle properties based on base level grid
				if (ud == 1)			% We are restoring from a session, so we know these already
					li = getappdata(h(k),'LineInfo');
					[x_inc, r] = strtok(li);	y_inc = strtok(r);
					x_inc = str2double(x_inc);	y_inc = str2double(y_inc);
				else
					resp = fix(abs(str2double(inputdlg({'Enter refinement factor'},'Refinement factor',[1 30],{'5'}))));
					if (isempty(resp) || isnan(resp) || resp == 0)		% OK, just make it a regular rectangle
						set(h(k), 'UIContextMenu', '')
						set_line_uicontext(h(k),'line')
						return
					end
					x_inc = handles.head(8) / resp;
					y_inc = handles.head(9) / resp;
				end
				setappdata(h(k),'LineInfo',[sprintf('%.16g %.16g',x_inc, y_inc), handles.head(7)]);
				do_nest_sizes = true;
			end
			uimenu(cmenuHand, 'Label', 'Rectangle limits (edit)', 'Callback', @rectangle_limits);
			ud = 1;		set(h(k),'UserData',ud)
		end
		if (handles.validGrid)			% Option only available to recognized grids
			uimenu(cmenuHand, 'Label', 'Crop Grid', 'Callback', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''CropaGrid_pure'')');
		end
% 		uimenu(cmenuHand, 'Label', 'Adjust to nesting dimensions', 'Sep', 'on', 'Callback', 'nesting_sizes(gcbo)');
		if (ud > 1 || handles.validGrid)
			uimenu(cmenuHand, 'Label', 'Show nesting info', 'Callback', 'nesting_sizes(gcbo,''Info'')');
			uimenu(cmenuHand, 'Label', 'Create blank grid', 'Callback', 'nesting_sizes(gcbo,''Grid'')');
		end
		if (k == numel(h))
			uimenu(cmenuHand, 'Label', 'New nested grid', 'Callback', 'nesting_sizes(gcbo,''New'')');
		end

		set_common_lineProps(h(k), cmenuHand, false)
		set(h(k),'Tag','NEST')
		setappdata(h(k),'RunCB',{'nesting_sizes', h(k)})	% Run this everytime rectangle is editted
		ui_edit_polygon(h(k))
		if (do_nest_sizes)
			nesting_sizes(h(k))			% Do the nesting size adjustment right away
		end
	end

% -----------------------------------------------------------------------------------------
function set_common_lineProps(h, cmenuHand, IS_PATCH)
% Common line properties set by some other functions

	cb_LineWidth = uictx_LineWidth(h);      % there are 5 cb_LineWidth outputs
	cb_solid = 'set(gco, ''LineStyle'', ''-''); refresh';
	cb_dashed = 'set(gco, ''LineStyle'', ''--''); refresh';
	cb_dotted = 'set(gco, ''LineStyle'', '':''); refresh';
	cb_dashdot = 'set(gco, ''LineStyle'', ''-.''); refresh';
	if (IS_PATCH),	cb_color = uictx_color(h,'EdgeColor');      % there are 9 cb_color outputs
	else,			cb_color = uictx_color(h);
	end
	setLineWidth(uimenu(cmenuHand, 'Label', 'Line Width', 'Sep','on'), cb_LineWidth)
	setLineStyle(uimenu(cmenuHand, 'Label', 'Line Style'), {cb_solid cb_dashed cb_dotted cb_dashdot})
	setLineColor(uimenu(cmenuHand, 'Label', 'Line Color'), cb_color)

% -----------------------------------------------------------------------------------------
function seismic_line(obj, evt, hL, opt)
% Plot seismicity projected along a polyline and enclosed inside a buffer zone of that pline
	handles = guidata(hL);
	hP = get(hL, 'UserData');	% handle to the buffer zone (if not exist, will be created later)
	if (opt(1) == 'b' || opt(1) == 'd')		% Create or expand a buffer zone
		if (opt(1) == 'b')
			resp = inputdlg('Width of buffer zone (deg)','Buffer width',[1 30],{'1.0'});
			if isempty(resp),	return,		end
			resp = str2double(resp{1}) / 2;	% Divide by 2 because this is actually the half-width`
		else
			resp = 0.5;			% The dafult value
		end
		[y, x] = buffer_j(get(hL, 'ydata'), get(hL, 'xdata'), resp, 'out', 13, 1);
		if (isempty(x)),	return,		end
		if (isempty(hP) || ~ishandle(hP))
			hP = patch('XData',x, 'YData',y, 'Parent',handles.axes1, 'EdgeColor',handles.DefLineColor, ...
				'FaceColor','none', 'LineWidth',handles.DefLineThick+1, 'Tag','SeismicBuffer');
			draw_funs(hP,'line_uicontext');	%uistack_j(hP,'bottom'),
			uistack_j(hL,'top')				% Put the 'Seismic' line on Top
			set(hL, 'UserData', hP)			% Save it there to ease its later retrival
		else
			set(hP,'XData',x, 'YData',y)	% Buffer already existed, just resize it
		end
		h = findobj(get(get(obj,'Parent'),'Parent') ,'Label','Project seismicity');
		set(h, 'Enable','on')
	else					% Project seismicity along the seismic line
		if (isempty(hP) || ~ishandle(hP))
			errordlg('Buffer zone is vanished, so what do you want me to do? Bye.','Error'),	return
		end
		hS = findobj(handles.axes1,'Tag','Earthquakes');
		if (isempty(hS))
			errordlg('Project what? The seismicity is gone. Bye Bye.','Chico Clever'),	return
		end

		x = get(hS,'XData');	y = get(hS,'YData');
		if (isa(x,'cell'))		% We don't want cells going into inpolygon
			x = [x{:}];			y = [y{:}];
		end
		IN = inpolygon(x,y, get(hP,'XData'),get(hP,'YData'));	% Find events inside the buffer zone
		x = x(IN);				y = y(IN);
		if (isempty(x))
			errordlg('Cou Cou!! There are no seisms inside this zone. Bye Bye.','Chico Clever'),	return
		end

		if (numel(hS) == 1)		% Must split between single and multiple (layered) sources
			evt_time = (getappdata(hS,'SeismicityTime'))';		% Get events time
			evt_dep  = (double(getappdata(hS,'SeismicityDepth')) / 10)';
			evt_mag  = (double(getappdata(hS,'SeismicityMag')) / 10)';
		else
			evt_time = zeros(numel(IN),1);		evt_dep = zeros(numel(IN),1);		evt_mag = zeros(numel(IN),1);
			k0 = 1;
			for (k = 1:numel(hS))
				t = getappdata(hS(k),'SeismicityTime');
				k1 = k0 + numel(t) - 1;
				evt_time(k0:k1) = t(:);
				t = double(getappdata(hS(k),'SeismicityDepth')) / 10;
				evt_dep(k0:k1) = t(:);
				t = double(getappdata(hS(k),'SeismicityMag')) / 10;
				evt_mag(k0:k1) = t(:);
				k0 = k1 + 1;
			end
		end
		evt_time = evt_time(IN);
		evt_dep = evt_dep(IN);
		evt_mag = evt_mag(IN);

		xL = get(hL, 'xdata');	yL = get(hL, 'ydata');

	    lat_i = yL(1:end-1);	lat_f = yL(2:end);
		lon_i = xL(1:end-1);	lon_f = xL(2:end);
		lineSegDist = vdist(lat_i,lon_i,lat_f,lon_f,handles.DefineEllipsoide([1 3])) / 1000;
		lineSegDist(isnan(lineSegDist)) = 0;				% Exact repeated points generate NaNs in solution
		lineSegDistAcum = [0; cumsum(lineSegDist(:))];
		lineSegDist(end+1) = 0;		lineSegDist = lineSegDist(:);	% Extend to ease the algo below

		f_name = [handles.path_tmp 'lixo.dat'];
		double2ascii(f_name,[xL(:) yL(:)],'%f\t%f');		% Save as file so we can use it mapproject
		ptFrac = c_mapproject([x(:) y(:)], ['-L' f_name '+']);	% Project and get fractional points along the line
		ptFrac = ptFrac.data(:,5) + 1;								% +1 because mapproject is 0 based
		intSeg = fix(ptFrac);
		rd = lineSegDistAcum(intSeg) + (ptFrac - intSeg) .* lineSegDist(intSeg);
		[rd, ind] = sort(rd);

		evt_time = evt_time(ind);	evt_dep = evt_dep(ind);
		figure;
		h = subplot(2,2,1);		plot(rd, evt_time, '.'),	set(h,'ylim',[min(evt_time) max(evt_time)],'xlim',[rd(1) rd(end)]);
		set(get(h,'Title'), 'String', 'Time vs distance (km)');
		h = subplot(2,2,2);		plot(rd, evt_dep, '.'),		set(h,'ylim',[min(evt_dep) max(evt_dep)+0.01],'xlim',[rd(1) rd(end)]);
		set(get(h,'Title'), 'String', 'Depth (km) vs distance (km)');
		inc = 25;	% Oroginal value. Keep it for backward compat.
		if (rd(end) < 250),	inc = rd(end) / 21;	end
		h = subplot(2,2,3);		histo_m('hist', rd, 0:inc:rd(end), [0 rd(end)]);
		set(get(h,'Title'), 'String', sprintf('Histogram with %.4g km bins', inc));
		h = subplot(2,2,4);		plot(rd, evt_mag, '.'),		set(h,'ylim',[min(evt_mag) max(evt_mag)],'xlim',[rd(1) rd(end)]);
		set(get(h,'Title'), 'String', 'Magnitude vs distance (km)');
	end

% -----------------------------------------------------------------------------------------
function hCopy = copy_line_object(obj, evt, hFig, hAxes)
% Make a copy of an line object and update some of its properties.
% If HCOPY output arg is specified, return without setting WindowButtonMotionFcn(). In this case one very
% likely need to transmit the original object in EVT, since this function was not called after a mouse event.

	handles = guidata(hFig);
	if (~isempty(evt) && handles.version7 < 8.4)
		oldH = evt;		% When calling this function by another that just wants a copy.
	else
		oldH = gco(hFig);
	end

	newH = copyobj(oldH,hAxes);
	h = findobj(get(newH,'uicontextmenu'),'label','Save line');
	if (~isempty(h))        % Replace the old line handle in the 'Save line' Callback by the just created one
		hFun = get(h,'Callback');
		hFun{2} = newH;
		set(h,'Callback',hFun)
	end
	if (isappdata(newH,'polygon_data'))
		rmappdata(newH,'polygon_data')		% Remove the parent's ui_edit_polygon appdata
	end
	if (isappdata(oldH,'cust_symb'))		% If copying an element with an associated GMT symbol copy it too
		cs_fname = getappdata(oldH, 'cust_symb');
		setappdata(newH, 'cust_symb', cs_fname)
	end

	state = uisuspend_j(hFig);				% Remember initial figure state
	x_lim = get(hAxes,'xlim');        y_lim = get(hAxes,'ylim');
	current_pt = get(hAxes, 'CurrentPoint');
	setappdata(newH,'old_pt',[current_pt(1,1) current_pt(1,2)])
	
	if (nargout)
		uirestore_j(state, 'nochildren');	% Restore the figure's initial state
		ui_edit_polygon(newH)				% Set the edition functions to the this handle
		hCopy = newH;
		return
	end

	set(hFig,'WindowButtonMotionFcn',{@wbm_MovePolygon,newH,[x_lim y_lim],hAxes},...
		'WindowButtonDownFcn',{@wbd_MovePolygon,newH,state}, 'Pointer','fleur');

% ---------
function wbm_MovePolygon(obj,evt,h,lim,hAxes)
	pt = get(hAxes, 'CurrentPoint');
	if (pt(1,1)<lim(1)) || (pt(1,1)>lim(2)) || (pt(1,2)<lim(3)) || (pt(1,2)>lim(4));   return; end
	old_pt = getappdata(h,'old_pt');
	xx = get(h,'XData');            yy = get(h,'YData');
	dx = pt(1,1) - old_pt(1);       dy = pt(1,2) - old_pt(2);
	xx = xx + dx;                   yy = yy + dy;
	setappdata(h,'old_pt',[pt(1,1) pt(1,2)])
	set(h, 'XData',xx, 'YData',yy);

% ---------
function wbd_MovePolygon(obj,evt,h,state)
	uirestore_j(state, 'nochildren');	% Restore the figure's initial state
	if (strcmp(get(h, 'Tag'), 'Arrow'))	% (Called by Copy object) Must also update the anchor points
		ud = get(h, 'UserData');
		x  = get(h, 'xdata');	y = get(h, 'ydata');
		ud.anchors = [(x(6)+x(7))/2 x(3); (y(6)+y(7))/2 y(3)]';
		set(h, 'UserData', ud);
	else
		ui_edit_polygon(h)				% Reset the edition functions with the correct handle
	end
% -----------------------------------------------------------------------------------------

% -----------------------------------------------------------------------------------------
function join_lines(obj,evt,hFig)
% Join lines that are NOT -- SEISPOLYGON, or MBTRACK, or FaultTrace or closed polygons

	hCurrLine = gco;
	hLines = get_polygon(hFig,'multi');				% Get the line handles
	if (isempty(hLines)),		return,		end
	hLines = setxor(hLines, hCurrLine);
	if (numel(hLines) == 0),	return,		end		% Nothing to join
	for (k = 1:numel(hLines))
		if (hCurrLine == hLines(k)),	continue,	end	% Must find a clever solution
		if (strcmp(get(hLines(k),'Type'),'patch')),		continue,	end
		[x, y, z, Z, was_closed] = join2lines([hCurrLine hLines(k)]);
		if (~was_closed),	delete(hLines(k)),	end		% Closed polygons are ignored 
		set(hCurrLine, 'XData',x, 'YData',y)
		if (~isempty(z)),	set(hCurrLine, 'ZData',z),	end
		if (~isempty(Z)),	setappdata(hCurrLine, 'ZData',Z),	end
	end

% ---------
function [x, y, z, Z, was_closed] = join2lines(hLines, TOL)
% Joint the two lines which have handles "hLines" by their closest connection points
% TOL is the max distance that the two lines can be apart and still be joint.
%	If not provided, defaults to Inf, but if it is and the two lines are too further
%	apart than X,Y & WAS_CLOSED are all set to empty. It's callers responsability to check it.

	if (nargin == 1),	TOL = inf;		end

	x1 = get(hLines(1),'XData');	x2 = get(hLines(2),'XData');
	y1 = get(hLines(1),'YData');	y2 = get(hLines(2),'YData');
	z1 = get(hLines(1),'ZData');	z2 = get(hLines(2),'ZData');
	Z1 = getappdata(hLines(1),'ZData');		Z2 = getappdata(hLines(2),'ZData');
	if (~isempty(z1) && isempty(z2)),		z2 = zeros(size(z1));
	elseif (~isempty(z2) && isempty(z1)),	z1 = zeros(size(z2));
	end
	if (~isempty(Z1) && isempty(Z2)),		Z2 = zeros(size(Z1));
	elseif (~isempty(Z2) && isempty(Z1)),	Z1 = zeros(size(Z2));
	end

	was_closed = false;
	if ((x2(1) == x2(end)) && (y2(1) == y2(end)))		% Ignore closed polygons
		x = x1;		y = y1;		z = [];		Z = [];
		was_closed = true;
		return
	end

	% Find how segments should be glued. That is find the closest extremities
	dif_x = [(x1(1) - x2(1)); (x1(1) - x2(end)); (x1(end) - x2(1)); (x1(end) - x2(end))];
	dif_y = [(y1(1) - y2(1)); (y1(1) - y2(end)); (y1(end) - y2(1)); (y1(end) - y2(end))];
	dist = sum([dif_x dif_y] .^2 ,2);	% Square of distances between the 4 extremities
	[mimi, I] = min(dist);				% We only care about the min location
	if (mimi > TOL)
		x = [];		y = [];		z = [];		Z = [];		was_closed = [];
		return
	end
	if (I == 1)				% Line 1 starts near the begining of line 2
		last = 1;
		if (x1(1) == x2(1) && y1(1) == y2(1)),		last = 2;		end
		x = [x2(end:-1:last) x1];	y = [y2(end:-1:last) y1];
		if (~isempty(z1)),	z = [z2(end:-1:last) z1];	end			% z is regular line 3rth dim, Z is the same stored in appdata
		if (~isempty(Z1)),	Z = [Z2(end:-1:last) Z1];	end
	elseif (I == 2)			% Line 2 ends near the begining of line 1 
		x = [x2 x1];				y = [y2 y1];
		if (~isempty(z1)),	z = [z2 z1];	end
		if (~isempty(Z1)),	Z = [Z2 Z1];	end
	elseif (I == 3)			% Line 1 ends near the begining of line 2
		x = [x1 x2];				y = [y1 y2];
		if (~isempty(z1)),	z = [z1 z2];	end
		if (~isempty(Z1)),	Z = [Z1 Z2];	end
	else					% Line 1 ends near the end of line 2
		last = 1;
		if (x1(end) == x2(end) && y1(end) == y2(end)),	last = 2;	end		% If points are equal, do not repeat it
		x = [x1 x2(end:-1:last)];	y = [y1 y2(end:-1:last)];
		if (~isempty(z1)),	z = [z1 z2(end:-1:last)];	end
		if (~isempty(Z1)),	Z = [Z1 Z2(end:-1:last)];	end
	end
	if (isempty(z1)),	z = [];	end
	if (isempty(Z1)),	Z = [];	end
% -----------------------------------------------------------------------------------------
	
% --------------------------------------------------------------------
function pt_quick_grd(hObj, evt, h, opt)
% Automatically calculate a grid from a PointZ data
	
	x = get(h, 'XData');	y = get(h, 'YData');	z = get(h, 'ZData');
	if (isempty(z))
		z = getappdata(h, 'ZData');
		if (isempty(z))		% Than sothing screwed in between
			errordlg('Smething screwed with storing the Z coords of this element. Have to quit.', 'Error')
			return
		end
	end
	z = double(z);
	if (isa(x, 'single'))
		x = double(x);		y = double(y);
	end
	x_min = min(x);			x_max = max(x);
	y_min = min(y);			y_max = max(y);

	if (opt == 'a')			% ------------- Estimate x|y_inc based on cheap heuristics ------
		dx = diff(x(1:min(300, numel(x))));
		dy = diff(y(1:min(300, numel(y))));
		dx(dx == 0) = [];		dy(dy == 0) = [];		% Remove influence of repeated points
		x_inc = abs(mean(dx));
		y_inc = abs(mean(dy));
		if (y_inc < x_inc/20),		y_inc = x_inc;		% Quite likely if data was previously gridded and dumped.
		elseif (x_inc < y_inc/20),	x_inc = y_inc;
		end

		nx = round((x_max - x_min) / x_inc) + 1;
		ny = round((y_max - y_min) / y_inc) + 1;
		if (nx > 1000 || ny > 1000)
			x_inc = (x_max - x_min) / min(nx,1000);
			y_inc = (y_max - y_min) / min(ny,1000);
			x_inc = (x_inc + y_inc) / 2;
			y_inc = x_inc;
			nx = round((x_max - x_min) / x_inc) + 1;
			ny = round((y_max - y_min) / y_inc) + 1;
		end
		if (x_inc < 1e-5 && y_inc < 1e-5)
			errordlg('Sorry, couldn''t figure out a sensible value off X & Y increments for the quick gridding.', 'Error')
			return
		end
	else
		resp = inputdlg({'Enter X increment:','Enter Y increment (if empty, equal x_inc):'},'Grid increments', 1);
		if (numel(resp) == 1 || (isempty(resp{1}) && isempty(resp{2})))		% User just Canceled
			return
		else
			% Allow only one answer
			if (isempty(resp{2}) && ~isempty(resp{1})),			resp{2} = resp{1};
			elseif (isempty(resp{1}) && ~isempty(resp{2})),		resp{1} = resp{2};
			end
			x_inc = str2double(resp{1});		y_inc = str2double(resp{2});
			if (isnan(x_inc) || x_inc <= 0 || isnan(y_inc) || y_inc <= 0)
				show_manguito,	return
			end
			nx = round((x_max - x_min) / x_inc) + 1;
			ny = round((y_max - y_min) / y_inc) + 1;
		end
	end

	opt_I = sprintf('-I%.8g/%.8g', x_inc, y_inc);
	opt_R = sprintf('-R%.12g/%.12g/%.12g/%.12g', x_min, x_min + (nx-1)*x_inc, y_min, y_min + (ny-1)*y_inc);
	[Z, head] = gmtmbgrid_m(x, y, z, opt_I, opt_R, '-Mz', '-C4');
	Z = single(Z);
	tmp.X = linspace(head(1), head(2), size(Z,2));
	tmp.Y = linspace(head(3), head(4), size(Z,1));
	tmp.head = head;
	tmp.name = 'Quick interpolated PointZ shape data';
	prjInfoStruc = aux_funs('getFigProjInfo', guidata(h));
	if (~isempty(prjInfoStruc.projWKT))
		tmp.srsWKT = prjInfoStruc.projWKT;
	elseif(~isempty(prjInfoStruc.proj4))
		tmp.srsWKT = ogrproj(prjInfoStruc.proj4);
	end
	mirone(Z, tmp)

% --------------------------------------------------------------------
function hh = loc_quiver(struc,varargin)
%QUIVER Quiver plot.
%   QUIVER(Struc,X,Y,U,V) plots velocity vectors as arrows with components (u,v)
%   at the points (x,y).  The matrices X,Y,U,V must all be the same size
%   and contain corresponding position and velocity components (X and Y
%   can also be vectors to specify a uniform grid).  QUIVER automatically
%   scales the arrows to fit within the grid.
%
%   QUIVER(Struc,X,Y,U,V,S) automatically scales the arrows to fit within the grid and
%   then stretches them by S. Use  S=0 to plot the arrows without the automatic scaling.
%
%	STRUC structure with this members
%		hQuiver - handles of a previously created arrow field (Default is [] )
%		spacingChanged - If different from zero previous arrows are deleted and reconstructed
%				  with the new input in varargin, but the handles remain valid (Default == 0). 
%		hAx - axes handle of the current figure. If empty a new fig is created
%		color - Optional field containing the line color. If absent, plot black lines.
%		thick - Optional field containing the line thickness. If absent, thick = 1.
%	To use the above default values, give and empty STRUC.
%
%   H = QUIVER(...) returns a vector of line handles.

	% Arrow head parameters
	alpha = 0.33;		% Size of arrow head relative to the length of the vector
	beta = 0.33;		% Width of the base of the arrow head relative to the length
	autoscale = 1;		% Autoscale if ~= 0 then scale by this.
	subsample = 1;		% Plot one every other grid node vector

	if (isempty(struc))
		hQuiver = [];		spacingChanged = [];		hAx = gca;		lc = 'k';	lThick = 1;
	else
		hQuiver = struc.hQuiver;		spacingChanged = struc.spacingChanged;		hAx = struc.hAx;
		if (isfield(struc, 'color')),	lc = struc.color;
		else,							lc = 'k';
		end
		if (isfield(struc, 'thick')),	lThick = struc.thick;
		else,							lThick = 1;
		end
	end

	nin = nargin - 1;

	% Check numeric input arguments
	if (nin < 4)					% quiver(u,v) or quiver(u,v,s)
		[msg,x,y,u,v] = xyzchk(varargin{1:2});
	else
		[msg,x,y,u,v] = xyzchk(varargin{1:4});
	end
	if ~isempty(msg), error(msg); end

	if (nin == 5)		% quiver(x,y,u,v,s)
		autoscale = varargin{nin};
	elseif  (nin == 6)
		autoscale = varargin{nin-1};
		subsample = abs(round(varargin{nin}));
	end

	% Scalar expand u,v
	if (numel(u) == 1),     u = u(ones(size(x))); end
	if (numel(v) == 1),     v = v(ones(size(u))); end

	if (subsample > 1)
		x = x(1:subsample:end,1:subsample:end);		y = y(1:subsample:end,1:subsample:end);
		u = u(1:subsample:end,1:subsample:end);		v = v(1:subsample:end,1:subsample:end);
	end

	if (autoscale)
		% Base autoscale value on average spacing in the x and y directions.
		% Estimate number of points in each direction as either the size of the
		% input arrays or the effective square spacing if x and y are vectors.
		if (min(size(x))==1),	n=sqrt(numel(x)); m=n;
		else,					[m,n]=size(x);
		end
		delx = diff([min(x(:)) max(x(:))])/n;
		dely = diff([min(y(:)) max(y(:))])/m;
		del = delx.^2 + dely.^2;
		if (del > 0)
			len = (u.^2 + v.^2)/del;
			maxlen = sqrt(max(len(:)));
		else
			maxlen = 0;
		end
		
		if maxlen > 0
			autoscale = autoscale*0.9 / maxlen;
		else
			autoscale = autoscale*0.9;
		end
		u = u*autoscale; v = v*autoscale;
	end

	% Make velocity vectors
	x = x(:).';		y = y(:).';
	u = u(:).';		v = v(:).';
	uu = [x; x+u; ones(size(u))*NaN];
	vv = [y; y+v; ones(size(u))*NaN];
	% Make arrow heads
	hu = [x+u-alpha*(u+beta*(v+eps)); x+u; x+u-alpha*(u-beta*(v+eps)); ones(size(u))*NaN];
	hv = [y+v-alpha*(v-beta*(u+eps)); y+v; y+v-alpha*(v+beta*(u+eps)); ones(size(v))*NaN];

	if (spacingChanged)
		try		delete(hQuiver),	hQuiver = [];	end		% Remove previous arrow field
	end

	if (isempty(hQuiver) || ~ishandle(hQuiver(1)))			% No arrows yet.
		h1 = line('XData',uu(:), 'YData',vv(:), 'Parent',hAx, 'Color',lc,'Linewidth',lThick);
		h2 = line('XData',hu(:), 'YData',hv(:), 'Parent',hAx, 'Color',lc,'Linewidth',lThick);
		if (nargout > 0),	hh = [h1;h2];	end
	else
		% We have the arrows and only want to change them
		set(hQuiver(1),'XData',uu(:), 'YData',vv(:))
		set(hQuiver(2),'XData',hu(:), 'YData',hv(:))
		if (nargout > 0),	hh = [hQuiver(1); hQuiver(2)];	end
	end

% -----------------------------------------------------------------------------------------
function set_country_uicontext(h)
% Minimalist patch uicontext to be used with countries patches due to the insane/ultrageous
% memory (and time) consumption taken by ML
	handles = guidata(h(1));
	for (i = 1:numel(h))
		cmenuHand = uicontextmenu('Parent',handles.figure1);
		set(h(i), 'UIContextMenu', cmenuHand);   
		uimenu(cmenuHand, 'Label', 'Save line', 'Callback', @save_line);
		uimenu(cmenuHand, 'Label', 'Delete', 'Callback', 'delete(gco)');
		if (handles.validGrid)
			item_ct = uimenu(cmenuHand, 'Label', 'ROI Crop Tools','Sep','on');
			uimenu(item_ct, 'Label', 'Crop Grid', 'Callback', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''CropaGrid_pure'')');
			uimenu(item_ct, 'Label', 'Set to constant', 'Callback', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''ROI_SetConst'')');
			uimenu(item_ct, 'Label', 'Mean', 'Callback', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''ROI_Mean'')');
			uimenu(item_ct, 'Label', 'Clip grid', 'Callback', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''ROI_Clip'')');
			uimenu(item_ct, 'Label', 'Median filter', 'Callback', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''ROI_MedianFilter'')');
			uimenu(item_ct, 'Label', 'Spline smooth', 'Callback', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''ROI_SplineSmooth'')');
			uimenu(item_ct, 'Label', 'Histogram', 'Callback', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''CropaGrid_histo'')');
		end
		cb_LineWidth = uictx_LineWidth(h(i));			% there are 5 cb_LineWidth outputs
		item_lw = uimenu(cmenuHand, 'Label', 'Line Width', 'Sep','on');
		uimenu(item_lw, 'Label', '1     pt', 'Callback', cb_LineWidth{1});
		uimenu(item_lw, 'Label', 'Other...', 'Callback', cb_LineWidth{5});
		item8 = uimenu(cmenuHand, 'Label','Fill Color', 'Sep','on');
		cb_color = uictx_color(h(i),'facecolor');		% there are 9 cb_color outputs
		uimenu(item8, 'Label', 'Other...', 'Callback', cb_color{9});
		uimenu(item8, 'Label', 'None', 'Callback', 'set(gco, ''FaceColor'', ''none'');refresh');
		uimenu(cmenuHand, 'Label', 'Transparency', 'Callback', @set_transparency);
		uimenu(cmenuHand, 'Label', 'Create Mask', 'Callback', 'poly2mask_fig(guidata(gcbo),gco)');
		if (handles.validGrid)
		end
		if (handles.image_type ~= 20)
			uimenu(cmenuHand, 'Label', 'Region-Of-Interest', 'Sep','on', 'Callback', ...
				'mirone(''DrawClosedPolygon_CB'',guidata(gcbo),gco)');
		end
	end

% -----------------------------------------------------------------------------------------
function okada_model(obj,eventdata,h,opt)
	if (nargin == 3),   opt = 'okada';		end
	hh = findobj('Tag','FaultTrace');		% Check if we have more than one (multi-segment?)faults
	if (isempty(hh)),   errordlg('This is just a line, NOT a fault trace. Can''t you see the difference?','Error'); return; end
	h_fig = get(0,'CurrentFigure');		handles = guidata(h_fig);

	% Guess minimum length segment that could be due to a bad line drawing
	if (handles.geog)
		min_len = 0.05;
	else
		imgLims = getappdata(handles.axes1,'ThisImageLims');
		if (abs(imgLims(2) - imgLims(1)) < 5000),   min_len = 5;    % Assume that the grid is in km
		else,                                        min_len = 5000; % Assume meters
		end
	end
	if (numel(hh) > 1)
		az = cell(1,numel(hh));
		for k=1:numel(hh)
			xx = get(hh(k),'XData');    yy = get(hh(k),'YData');
			dx = diff(xx);  dy = diff(yy);              dr = sqrt(dx.*dx + dy.*dy);
			ind = find(dr < min_len);
			if (~isempty(ind) && length(xx) > 2)		% Remove too short segments
				xx(ind) = [];   yy(ind) = [];
				set(hh(k),'XData',xx,'YData',yy)
			end
			azim = show_lineAzims([],[],hh(k));
			az{k} = azim.az;
		end
		h = hh;
	else
		xx = get(h,'XData');    yy = get(h,'YData');
		dx = diff(xx);			dy = diff(yy);			dr = sqrt(dx.*dx + dy.*dy);
		ind = find(dr < min_len);
		if (~isempty(ind) && length(xx) > 2)			% Remove too short segments
			xx(ind) = [];   yy(ind) = [];
			set(h,'XData',xx,'YData',yy)
		end
		azim = show_lineAzims([],[],h);
		az = azim.az;
	end

	if (strcmp(opt,'okada')),			deform_okada(handles,h,az);
	elseif (strcmp(opt,'mansinha')),	deform_mansinha(handles,h,az);
	end
	% Feigl's example
	% u1 = -22;    u2 = 514;     u3 = 0;
	% W = 3.1;    depth = 3.4;    dip = 28;

% -----------------------------------------------------------------------------------------
function set_SRTM_rect_uicontext(h,opt)
	% h is a handle to a line object (that can be closed)
	handles = guidata(h(1));	cmenuHand = uicontextmenu('Parent',handles.figure1);
	set(h, 'UIContextMenu', cmenuHand);
	ui_edit_polygon(h)    % Set edition functions
	uimenu(cmenuHand, 'Label', 'Delete', 'Callback', 'delete(gco)');
	cb_Fill_surface = 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''FillGaps'',''surface'');delete(gco)';
	cb_Fill_cubic = 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''FillGaps'',''cubic'');delete(gco)';
	cb_Fill_linear = 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''FillGaps'',''linear'');delete(gco)';
	cb_Fill_sea   = 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''FillGaps'',''sea'');delete(gco)';
	uimenu(cmenuHand, 'Label', 'Fill gaps (surface)', 'Callback', cb_Fill_surface);
	uimenu(cmenuHand, 'Label', 'Fill gaps (cubic)', 'Callback', cb_Fill_cubic);
	uimenu(cmenuHand, 'Label', 'Fill gaps (linear)', 'Callback', cb_Fill_linear);
	uimenu(cmenuHand, 'Label', 'Fill gaps (sea)', 'Callback', cb_Fill_sea);

% -----------------------------------------------------------------------------------------
function set_ContourLines_uicontext(h,h_label)
% h is the handle to the contour value. Each contour is given this uicontext
	handles = guidata(h(1));	cmenuHand = uicontextmenu('Parent',handles.figure1);
	set(h, 'UIContextMenu', cmenuHand);
	% cb1     = 'mirone(''DrawEditLine_CB'',gcbo,[],guidata(gcbo))';
	ui_edit_polygon(h)				% Set edition functions
	cb_rac = {@remove_symbolClass,h};		% It will also remove the labels because they have the same tag.
	cb_LineWidth = uictx_LineWidth(h);		% there are 5 cb_LineWidth outputs
	cb18 = 'set(gco, ''LineStyle'', ''-''); refresh';   cb19 = 'set(gco, ''LineStyle'', ''--''); refresh';
	cb20 = 'set(gco, ''LineStyle'', '':''); refresh';   cb21 = 'set(gco, ''LineStyle'', ''-.''); refresh';
	cb_color = uictx_color(h);				% there are 9 cb_color outputs

	uimenu(cmenuHand, 'Label', 'Delete contour', 'Callback',{@remove_singleContour,h});
	uimenu(cmenuHand, 'Label', 'Delete all contours', 'Callback', cb_rac);
	% item1 = uimenu(cmenuHand, 'Label', 'Edit contour (left-click on it)', 'Callback', cb1);
	uimenu(cmenuHand, 'Label', 'Join contours', 'Callback', {@join_lines,handles.figure1});
	uimenu(cmenuHand, 'Label', 'Save contour', 'Callback', {@save_formated,h});
	uimenu(cmenuHand, 'Label', 'Contour length', 'Callback', {@show_LineLength,[]});
	uimenu(cmenuHand, 'Label', 'Area under contour', 'Callback', @show_Area);
	item_lw = uimenu(cmenuHand, 'Label', 'Contour Line Width', 'Sep','on');
	setLineWidth(item_lw,cb_LineWidth)
	item_ls = uimenu(cmenuHand, 'Label', 'Contour Line Style');
	setLineStyle(item_ls,{cb18 cb19 cb20 cb21})
	item_lc = uimenu(cmenuHand, 'Label', 'Contour Line Color');
	setLineColor(item_lc,cb_color)
	cb_CLineWidth = uictx_Class_LineWidth(h);           % there are 5 cb_CLineWidth outputs
	item8 = uimenu(cmenuHand, 'Label', 'All Contours Line Width', 'Sep','on');
	uimenu(item8, 'Label', '1       pt', 'Callback', cb_CLineWidth{1});
	uimenu(item8, 'Label', '2       pt', 'Callback', cb_CLineWidth{2});
	uimenu(item8, 'Label', '3       pt', 'Callback', cb_CLineWidth{3});
	uimenu(item8, 'Label', 'Other...', 'Callback', cb_CLineWidth{5});
	cb_CLineStyle = uictx_Class_LineStyle(h);        % there are 4 cb_CLineStyle outputs
	item9 = uimenu(cmenuHand, 'Label', 'All Contours Line Style');
	uimenu(item9, 'Label', 'solid', 'Callback', cb_CLineStyle{1});
	uimenu(item9, 'Label', 'dashed', 'Callback', cb_CLineStyle{2});
	uimenu(item9, 'Label', 'dotted', 'Callback', cb_CLineStyle{3});
	uimenu(item9, 'Label', 'dash-dotted', 'Callback', cb_CLineStyle{4});
	cb_CLineColor = uictx_Class_LineColor(h);              % there are 9 cb_CLineColor outputs
	item_lc = uimenu(cmenuHand, 'Label', 'All Contours Line Color');
	setLineColor(item_lc,cb_CLineColor)

% -----------------------------------------------------------------------------------------
function setCoastLineUictx(h)
% h is a handle to a line object (a CoastLine one)
	tag = get(h,'Tag');
	if (strcmp(tag,'CoastLineNetCDF')),			label = 'Delete coastlines';
	elseif (strcmp(tag,'PoliticalBoundaries')),	label = 'Delete boundaries';
	elseif (strcmp(tag,'Rivers')),				label = 'Delete rivers';
	end
	handles = guidata(h);
	cmenuHand = uicontextmenu('Parent',handles.figure1);
	set(h, 'UIContextMenu', cmenuHand)
	cb_LineWidth = uictx_LineWidth(h);		% there are 5 cb_LineWidth outputs
	cb13 = 'set(gco, ''LineStyle'', ''-''); refresh';   cb14 = 'set(gco, ''LineStyle'', ''--''); refresh';
	cb15 = 'set(gco, ''LineStyle'', '':''); refresh';   cb16 = 'set(gco, ''LineStyle'', ''-.''); refresh';
	cb_color = uictx_color(h);				% there are 9 cb_color outputs
	
	uimenu(cmenuHand, 'Label', label, 'Callback', 'delete(gco)');
	uimenu(cmenuHand, 'Label', 'Save coastline', 'Callback', {@save_formated,h});

	if (strcmp(tag,'CoastLineNetCDF'))		% Options to apply ocean/land masking
		item = uimenu(cmenuHand, 'Label', 'Mask', 'Sep','on');
		uimenu(item, 'Label', 'Land',  'Callback', {@apply_grdlandMask,h, 'L'})
		uimenu(item, 'Label', 'Ocean', 'Callback', {@apply_grdlandMask,h, 'O'})
	end

	setLineWidth(uimenu(cmenuHand, 'Label', 'Line Width', 'Sep','on'), cb_LineWidth)
	item_ls = uimenu(cmenuHand, 'Label', 'Line Style');
	setLineStyle(item_ls,{cb13 cb14 cb15 cb16})
	setLineColor(uimenu(cmenuHand, 'Label', 'Line Color'), cb_color)
	ui_edit_polygon(h)

% -----------------------------------------------------------------------------------------
function apply_grdlandMask(hObj, evt, h, opt)
% Blank current grid/image with the default's grdlandmask setting for this coastline resolution
	handles = guidata(h);
	if (~handles.is_projected)
		opt_R = sprintf('-R%.12g/%.12g/%.12g/%.12g',handles.head(1), handles.head(2), handles.head(3), handles.head(4));
		opt_I = sprintf('-I%0.12g/%0.12g',handles.head(8), handles.head(9));
	else			% Now this adds quite a lot of processing to the job
		wkt = aux_funs('get_proj_string', handles, 'wkt');
		if (isempty(wkt))
			errordlg('Figure is projected, but failed to get proj info.', 'Error'),		return
		end
		data = [handles.head(1) handles.head(3); handles.head(1) handles.head(4); ...	% Use the 4 corners
		        handles.head(2) handles.head(4); handles.head(2) handles.head(3)];		
		[out, msg, opt_R] = proj2proj_pts(handles, data, 'srcWKT',wkt, 'dstProj4','+proj=longlat');
		if (~isempty(msg)),		errordlg(msg, 'Error'),		return,		end
		inc = handles.head(8:9) / (handles.EarthRad * 1000) * 180 / pi;
		opt_I = sprintf('-I%0.12g/%0.12g',inc(1), inc(2));
	end
	opt_F = ' ';	opt_e = ' ';
	if (handles.head(7)),		opt_F = '-F';	end
	if (handles.IamCompiled),	opt_e = '-e';	end
	opt_D = sprintf('-D%s', getappdata(h,'resolution'));	% Get the resolution as stored in line's appdata
	opt_N = '-N0/255';			% Continent masking
	set(handles.figure1,'pointer','watch'),		pause(0.01)
	[mask, new_hdr] = c_grdlandmask(opt_R, opt_I, opt_D, opt_N, opt_F, opt_e, '-A0/0/1', '-V');
	mask = uint8(mask);
	if (handles.is_projected)	% Then MASK is in geogs and we need to project it.
		hdrStruct.ULx = new_hdr(1) - new_hdr(8) / 2 * (~new_hdr(7));	% Goto pixel reg
		hdrStruct.ULy = new_hdr(4) + new_hdr(9) / 2 * (~new_hdr(7));
		hdrStruct.Xinc = new_hdr(8);		hdrStruct.Yinc = new_hdr(9);
		hdrStruct.ResampleAlg = 'bilinear';
		hdrStruct.t_res = handles.head(8:9);
		hdrStruct.SrcProjSRS = '+proj=longlat';
		hdrStruct.DstProjSRS = wkt;
		[mask, att] = gdalwarp_mex(mask, hdrStruct);	% The output is normally larger then base grid
		rect_crop = [handles.head(1) handles.head(3) diff(handles.head(1:2)) diff(handles.head(3:4))];
		[mask,r_c] = cropimg(att.GMT_hdr(1:2), att.GMT_hdr(3:4), mask, rect_crop, 'out_grid');
	end

	if (opt == 'O')				% Ocean masking.
		set(handles.hImg, 'AlphaData', mask)
		mask = ~logical(mask);
	else
		mask = logical(mask);
		new_mask = alloc_mex(size(mask,1), size(mask,2), 'uint8', 255);
		new_mask(mask) = 0;
		set(handles.hImg, 'AlphaData', new_mask)
	end

	img = get(handles.hImg, 'CData');	% Mask the image as well
	if (ndims(img) == 3)
		bg_color = uint8(handles.bg_color * 255);
		tmp = img(:,:,1);		tmp(mask) = bg_color(1);	img(:,:,1) = tmp;
		tmp = img(:,:,2);		tmp(mask) = bg_color(2);	img(:,:,2) = tmp;
		tmp = img(:,:,3);		tmp(mask) = bg_color(3);	img(:,:,3) = tmp;
	else
		img(mask) = 0;
		pal = get(handles.figure1, 'Colormap');
		pal(1,:) = handles.bg_color;
		set(handles.figure1, 'Colormap', pal);
	end
	set(handles.hImg, 'CData', img)
	set(handles.figure1,'pointer','arrow');

	if (handles.validGrid)
		[X,Y,Z] = load_grd(handles);	Z(mask) = NaN;
		setappdata(handles.figure1,'dem_z',Z);
		handles.have_nans = grdutils(Z,'-N');
		zz = grdutils(Z,'-L');			handles.head(5:6) = [zz(1) zz(2)];
		handles.firstIllum = true;
		set(handles.haveNaNs,'Vis','on')
		guidata(handles.figure1, handles)
	end
% -----------------------------------------------------------------------------------------

% -----------------------------------------------------------------------------------------
function set_PB_uicontext(h,data)
% h is a handle to the lines of the PB_All (P. Bird Plate Boundaries) object

for i = 1:7     % Loop over all Plate Boundaries Types
	h_cur = [];
	switch i
		case 1,			h_cur = h.OSR;  data_cur = data.OSR;    % class = 'OSR'
		case 2,			h_cur = h.OTF;  data_cur = data.OTF;    % class = 'OTF'
		case 3,			h_cur = h.CRB;  data_cur = data.CRB;    % class = 'CRB'
		case 4,			h_cur = h.CTF;  data_cur = data.CTF;    % class = 'CTF'
		case 5,			h_cur = h.CCB;  data_cur = data.CCB;    % class = 'CCB'
		case 6,			h_cur = h.OCB;  data_cur = data.OCB;    % class = 'OCB'
		case 7,			h_cur = h.SUB;  data_cur = data.SUB;    % class = 'SUB'
	end
	if (isempty(h_cur)),	continue,	end
	cmenuHand = uicontextmenu;
	set(h_cur, 'UIContextMenu', cmenuHand);
	cb_LineWidth = uictx_Class_LineWidth(h_cur);    % there are 5 cb_PB_LineWidth outputs
	cb_color = uictx_Class_LineColor(h_cur);        % there are 9 cb_PB_color outputs
	uimenu(cmenuHand, 'Label', 'Segment info', 'Callback', {@PB_All_Info,h_cur,data_cur});
	uimenu(cmenuHand, 'Label', 'Delete class', 'Callback', 'delete(findobj(''Tag'',''PB_All''))', 'Sep','on');
	uimenu(cmenuHand, 'Label', 'Segment length', 'Callback', {@show_LineLength,[]});
	item3 = uimenu(cmenuHand, 'Label', 'Line Width', 'Sep','on');
	uimenu(item3, 'Label', '2       pt', 'Callback', cb_LineWidth{2});
	uimenu(item3, 'Label', '3       pt', 'Callback', cb_LineWidth{3});
	uimenu(item3, 'Label', '4       pt', 'Callback', cb_LineWidth{4});
	uimenu(item3, 'Label', 'Other...', 'Callback', cb_LineWidth{5});
	item_lc = uimenu(cmenuHand, 'Label', 'Color');
	setLineColor(item_lc,cb_color)
end

% -----------------------------------------------------------------------------------------
function set_isochrons_uicontext(h, data)
% h are handles to the lines of isochrons (or other lines with a info)
	if (isempty(h)),	return,		end
	tag = get(h,'Tag');
	if (iscell(tag)),   tag = tag{1};   end

	handles = guidata(get(h(1),'Parent'));				% Get Mirone handles
	cmenuHand = uicontextmenu('Parent',handles.figure1, 'Callback', @store_clicked_pt);
	set(h, 'UIContextMenu', cmenuHand);
	cb_LineWidth = uictx_LineWidth(h);		% there are 5 cb_LineWidth outputs
	cb_color = uictx_color(h);				% there are 9 cb_color outputs
	cbls1 = 'set(gco, ''LineStyle'', ''-''); refresh';   cbls2 = 'set(gco, ''LineStyle'', ''--''); refresh';
	cbls3 = 'set(gco, ''LineStyle'', '':''); refresh';   cbls4 = 'set(gco, ''LineStyle'', ''-.''); refresh';
	if (~isempty(data) && ~all(isempty(cat(2,data{:}))) )
		uimenu(cmenuHand, 'Label', [tag ' info'], 'Callback', {@Isochrons_Info,data});
		uimenu(cmenuHand, 'Label', ['Delete this ' tag ' line'], 'Callback', {@del_line,h}, 'Sep','on');
	else
		uimenu(cmenuHand, 'Label', ['Delete this ' tag ' line'], 'Callback', {@del_line,h});
	end
	uimenu(cmenuHand, 'Label', ['Save this '  tag ' line'],  'Callback', @save_line);
	uimenu(cmenuHand, 'Label', 'Join lines', 'Callback', {@join_lines,handles.figure1});
	uimenu(cmenuHand, 'Label', 'Line azimuths', 'Callback', @show_lineAzims);
	uimenu(cmenuHand, 'Label', 'Line length', 'Callback', {@show_LineLength,[],'nikles'});
	LINE_ISCLOSED = 0;
	for (i = 1:numel(h))
		x = get(h(i),'XData');      y = get(h(i),'YData');
		if (numel(x) > 2 && (x(1) == x(end)) && (y(1) == y(end)))		% See if we have at least one closed line
			LINE_ISCLOSED = 1;		break
		end
	end
	% If at least one is closed, activate the Area option
	if (LINE_ISCLOSED)
		uimenu(cmenuHand, 'Label', 'Area under polygon', 'Callback', @show_Area, 'Sep','on');
		uimenu(cmenuHand, 'Label', 'Create Mask', 'Callback', 'poly2mask_fig(guidata(gcbo),gco)');
		set_ROI_uicontext(handles, cmenuHand, false);
	end
	if (handles.image_type ~= 20 && (ndims(get(handles.hImg,'CData')) == 2 || handles.validGrid))
		if (handles.nLayers > 1)
			cbTrack = 'setappdata(gcf,''TrackThisLine'',gco); mirone(''ExtractProfile_CB'',guidata(gcbo),''3D'')';
			uimenu(cmenuHand, 'Label', '3D interpolation', 'Callback', cbTrack);
		end
		cbTrack = 'setappdata(gcf,''TrackThisLine'',gco); mirone(''ExtractProfile_CB'',guidata(gcbo),''point'')';
		uimenu(cmenuHand, 'Label', 'Point interpolation', 'Callback', cbTrack, 'Sep','on');
		cbTrack = 'setappdata(gcf,''TrackThisLine'',gco); mirone(''ExtractProfile_CB'',guidata(gcbo))';
		uimenu(cmenuHand, 'Label', 'Extract profile', 'Callback', cbTrack);
	end
	setLineWidth(uimenu(cmenuHand, 'Label', 'Line Width', 'Sep','on'), cb_LineWidth)
	setLineStyle(uimenu(cmenuHand, 'Label', 'Line Style'), {cbls1 cbls2 cbls3 cbls4})
	setLineColor(uimenu(cmenuHand, 'Label', 'Color'), cb_color)
	set_stack_order(cmenuHand)		% Set options to change order in the stackpot	

	% --------- Now set the class properties
	item_allThis = uimenu(cmenuHand, 'Label', ['All ' tag '...'], 'Sep','on');
	uimenu(item_allThis, 'Label', 'Delete all', 'Callback', {@remove_symbolClass,h});
	uimenu(item_allThis, 'Label', 'Save all', 'Callback',   {@save_line,h});

	cb_ClassColor = uictx_Class_LineColor(h);        % there are 9 cb_color outputs. SHIT, IF h IS MULTIPLE WE CAN'T SET INDIVIDUAL COLORS
	setLineColor(uimenu(item_allThis, 'Label', 'Color'), cb_ClassColor)
	cb_ClassLineWidth = uictx_Class_LineWidth(h);    % there are 5 cb_ClassLineWidth outputs
	item_Class_lw = uimenu(item_allThis, 'Label', 'Line Width');
	uimenu(item_Class_lw, 'Label', '1       pt', 'Callback', cb_ClassLineWidth{1});
	uimenu(item_Class_lw, 'Label', '2       pt', 'Callback', cb_ClassLineWidth{2});
	uimenu(item_Class_lw, 'Label', '3       pt', 'Callback', cb_ClassLineWidth{3});
	uimenu(item_Class_lw, 'Label', '4       pt', 'Callback', cb_ClassLineWidth{4});
	uimenu(item_Class_lw, 'Label', 'Other...', 'Callback', cb_ClassLineWidth{5});
	cb_ClassLineStyle = uictx_Class_LineStyle(h);    % there are 4 cb_ClassLineStyle outputs
	item_Class_lt = uimenu(item_allThis, 'Label', 'Line Style');
	setLineStyle(item_Class_lt,{cb_ClassLineStyle{1} cb_ClassLineStyle{2} cb_ClassLineStyle{3} cb_ClassLineStyle{4}})
	item = uimenu(cmenuHand, 'Label', 'Isochron pole operations', 'Sep','on');
	uimenu(item, 'Label', 'Update poles in headers', 'Callback', {@pole2neighbor, [], 'update_poles'});
	item2 = uimenu(item, 'Label', 'Compute pole to neighbor', 'Sep','on');
	uimenu(item2, 'Label', 'Best-Fit (all in this plate)', 'Callback', {@pole2neighbor, [], 'anglefit'});
	uimenu(item2, 'Label', 'Best-Fit (only me)', 'Callback', {@pole2neighbor, [], 'anglefit', 1});
	uimenu(item2, 'Label', 'Best-Fit (only me -> iterate)', 'Callback', {@pole2neighbor, [], 'anglefit', 10});
	uimenu(item2, 'Label', 'Show Results', 'Callback', {@pole2neighbor, [], 'showresults'});
	uimenu(item, 'Label', 'Plates stage poles', 'Callback', {@pole2neighbor, [], 'plate_stages'});
	uimenu(item, 'Label', 'Reconstruct at this age', 'Callback', {@pole2neighbor, [], 'reconst'}, 'Sep','on');
	uimenu(item, 'Label', 'Age Grid', 'Callback', {@pole2neighbor, [], 'agegrid'}, 'Sep','on');
	uimenu(item, 'Label', 'Make Age-script', 'Callback', @make_age_script);	
	uimenu(cmenuHand, 'Label', 'Euler rotation', 'Sep','on', 'Callback', 'euler_stuff(gcf,gco)');
	uimenu(cmenuHand, 'Label', 'Rotate all in plate', 'Callback', 'pole2neighbor([],[], gco, ''rot_all'')');
	for (i=1:length(h)),		ui_edit_polygon(h(i)),		end		% Set edition functions

% -----------------------------------------------------------------------------------------
function store_clicked_pt(hObject, evt)
% Store the clicked point on the base cmenuHand appdata. This is usefull when one want to know
% the current point at the time of first clicking on figure. The axes current point prop only
% tells the current point of last clicking, which normally was somewhere on a uictx option.
	handles = guidata(hObject);
	setappdata(hObject, 'clicked_pt', get(handles.axes1, 'CurrentPoint'))
	
% -----------------------------------------------------------------------------------------
function set_gmtfile_uicontext(h, data)
% h is a handle to the line of a gmtfile
	tag = get(h,'Tag');
	if (iscell(tag)),   tag = tag{1};   end

	handles = guidata(h(1));	cmenuHand = uicontextmenu('Parent',handles.figure1);
	set(h, 'UIContextMenu', cmenuHand);
	cb_LineWidth = uictx_LineWidth(h);       % there are 5 cb_LineWidth outputs
	cb_color = uictx_color(h);               % there are 9 cb_color outputs
	cbls1 = 'set(gco, ''LineStyle'', ''-''); refresh';   cbls2 = 'set(gco, ''LineStyle'', ''--''); refresh';
	cbls3 = 'set(gco, ''LineStyle'', '':''); refresh';   cbls4 = 'set(gco, ''LineStyle'', ''-.''); refresh';
	uimenu(cmenuHand, 'Label', [tag ' info'], 'Callback', {@gmtfile_Info,h,data});
	uimenu(cmenuHand, 'Label', ['Delete this ' tag ' line'], 'Callback', 'delete(gco)', 'Sep','on');
	uimenu(cmenuHand, 'Label', ['Save this ' tag ' line'], 'Callback', @save_line);
	uimenu(cmenuHand, 'Label', 'Open with gmtedit', 'Callback', {@call_gmtedit,h});
 	uimenu(cmenuHand, 'Label', 'Extract Mag Chunk', 'Callback', {@call_gmtedit,h,'nikles'});
	uimenu(cmenuHand, 'Label', 'Create Mask', 'Callback', 'poly2mask_fig(guidata(gcbo),gco)');
	deal_opts('mgg_coe', cmenuHand);
	%uimenu(cmenuHand, 'Label', 'Try to relocate', 'Callback', {@tryRelocate,h});
	setLineWidth(uimenu(cmenuHand, 'Label', 'Line Width', 'Sep','on'), cb_LineWidth)
	setLineStyle(uimenu(cmenuHand, 'Label', 'Line Style'), {cbls1 cbls2 cbls3 cbls4})
	setLineColor(uimenu(cmenuHand, 'Label', 'Color'), cb_color)

% -----------------------------------------------------------------------------------------
% function tryRelocate(obj,evt,h)
% 	handles = guidata(obj);
% 	[X,Y,Z] = load_grd(handles);
% 	if (isempty(Z))		return,		end
% 	hg = gmtedit(getappdata(h,'FullName'));
% 	set(hg, 'Vis','off')								% Hide it
% 	handGmtedit = guidata(hg);
% 	y_m = get(handGmtedit.h_mm,'YData');	% Get the Mag values
% 	delete(hg)
% 	x = get(h, 'XData');		y = get(h, 'YData');
% 	new_x = x;					new_y = y;
% 	dr = 0.005;
% 	for (j = 1:numel(x))
% 		if (isnan(x(j)))		continue,		end
% 		yy_c =  y(j) + [-5:5] * dr;
% 		xx_c = repmat(x(j), numel(yy_c), 1);
% 		zc = grdtrack_m(Z, handles.head, [xx_c yy_c'], '-Z', '-Q');
% 		xx_r =  x(j) + [-5:5] * dr;
% 		yy_r = repmat(y(j), numel(xx_r), 1);
% 		zr = grdtrack_m(Z, handles.head, [xx_r' yy_r], '-Z', '-Q');
% 		dcol = abs(zc - y_m(j));
% 		drow = abs(zr - y_m(j));
% 		[this_min_c, n] = min(dcol);
% 		[this_min_r, m] = min(drow);
% 		if (this_min_c < this_min_r)
% 			new_x(j) = x(j);		new_y(j) = yy_c(n);
% 		else
% 			new_x(j) = xx_r(m);		new_y(j) = y(j);
% 		end
% 	end
% 	set(h, 'XData', new_x, 'YData', new_y);

% -----------------------------------------------------------------------------------------
function call_gmtedit(obj, evt, h, opt)
	handles = guidata(h);
	if (handles.is_projected)
		proj = aux_funs('get_proj_string', handles.figure1, 'proj');
	end
	if (nargin == 4)		% Call helper window to extract a chunk of the mag anom profile
		hFig = get(get(h,'Parent'),'Parent');
		[xp,yp] = getline_j(hFig);
		if (numel(xp) < 2),		return,		end
		if (handles.is_projected)
			xy = proj2proj_pts([], [xp(:) yp(:)], 'srcProj4', proj, 'dstProj4', '+proj=longlat');
			xp = xy(:,1);	yp = xy(:,2);
			bak = handles.geog;		handles.geog = 1;		% Need this to cheat a later call to mag_synthetic
			guidata(handles.figure1, handles)
		end
		hGL = line('XData', xp, 'YData', yp,'Color','y','Parent',get(h,'Parent'),'LineWidth',3,'Tag','polyline');
		guidelineAzim = azimuth_geo(yp(1), xp(1), yp(end), xp(end));
		mag_synthetic(hFig, h, hGL, guidelineAzim)
		if (handles.is_projected)
			handles.geog = bak;		guidata(handles.figure1, handles)
		end
		return
	end
	pt = get(get(h,'Parent'), 'CurrentPoint');
	if (handles.is_projected)
		pt = proj2proj_pts([], pt(1,1:2), 'srcProj4', proj, 'dstProj4', '+proj=longlat');
	end
	vars = getappdata(h,'VarsName');		opt_V = '  ';	% To be ignored opt_V needs to have at least 2 chars
	if (~isempty(vars))
		opt_V = ['-V' vars{1} ','  vars{2} ',' vars{3}];	% Need to encode the Vars info in a single string
		if (strcmp(opt_V,'-V,,')),		opt_V = '  ';	end	% When vars is actually a 3 empties cell
	end
	gmtedit(getappdata(h,'FullName'), sprintf('-P%.6f/%.6f',pt(1,1:2)), opt_V, get(h,'Parent'));

% -----------------------------------------------------------------------------------------
function cb = uictx_setMarker(h,prop)
% Set uicontext colors in a PB object class hose handles are contained in h
% PROP is either "Marker" or "MarkerSize". OPT is either the symbol or it's size
	if (strcmp(prop,'Marker'))
		cb{1} = {@other_Marker,h,prop,'+'};		cb{2} = {@other_Marker,h,prop,'o'};
		cb{3} = {@other_Marker,h,prop,'*'};		cb{4} = {@other_Marker,h,prop,'.'};
		cb{5} = {@other_Marker,h,prop,'x'};		cb{6} = {@other_Marker,h,prop,'s'};
		cb{7} = {@other_Marker,h,prop,'d'};		cb{8} = {@other_Marker,h,prop,'^'};
		cb{9} = {@other_Marker,h,prop,'v'};		cb{10} = {@other_Marker,h,prop,'>'};
		cb{11} = {@other_Marker,h,prop,'<'};	cb{12} = {@other_Marker,h,prop,'p'};
		cb{13} = {@other_Marker,h,prop,'h'};
	elseif (strcmp(prop,'MarkerSize'))
		cb{1} = {@other_Marker,h,prop,7};		cb{2} = {@other_Marker,h,prop,8};
		cb{3} = {@other_Marker,h,prop,9};		cb{4} = {@other_Marker,h,prop,10};
		cb{5} = {@other_Marker,h,prop,12};		cb{6} = {@other_Marker,h,prop,14};
		cb{7} = {@other_Marker,h,prop,16};		cb{8} = {@other_SymbSize,h};
	end

	function other_Marker(obj,eventdata,h,prop,opt)
		set(h,prop,opt);    refresh

% -----------------------------------------------------------------------------------------
function set_font_size(h, opt)
% Set font size for axes lables when H = gca and nargin == 3 or any other property that fits
% in the H (handle), OPT (property).
	resp  = inputdlg({'Enter new font size (pt)'}, 'Font size', [1 30]);
	resp = str2double(resp);
	if isnan(resp),		return,		end
	if (nargin == 1),	opt = 'FontSize';	end
	set(h,opt,resp);

% -----------------------------------------------------------------------------------------
function cb = uictx_Class_LineWidth(h)
% Set uicontext LineWidths in a PB object class hose handles are contained in h
	cb{1} = {@other_Class_LineWidth,h,1};		cb{2} = {@other_Class_LineWidth,h,2};
	cb{3} = {@other_Class_LineWidth,h,3};		cb{4} = {@other_Class_LineWidth,h,4};
	cb{5} = {@other_Class_LineWidth,h,[]};

function other_Class_LineWidth(obj,eventdata,h,opt)
% If individual Lines were previously removed (by "Remove Line") h has invalid
% handles, so make sure all handles are valid
	h=h(ishandle(h));
	if ~isempty(opt)   
		set(h,'LineWidth',opt);        refresh;
	else
		resp  = inputdlg({'Enter new line width (pt)'}, 'Line width', [1 30]);
		resp = str2double(resp);
		if isnan(resp),		return,		end
		set(h,'LineWidth',str2double(resp));        refresh
	end
% -----------------------------------------------------------------------------------------

% -----------------------------------------------------------------------------------------
function cb = uictx_Class_LineColor(h,prop)
% Set uicontext colors in a PB object class whose handles are contained in h.
% PROP, when given, is either "MarkerFaceColor" or "MarkerEdgeColor"
	if (nargin == 1),   prop = [];  end
	cb{1} = {@other_Class_LineColor,h,'k',prop};       cb{2} = {@other_Class_LineColor,h,'w',prop};
	cb{3} = {@other_Class_LineColor,h,'r',prop};       cb{4} = {@other_Class_LineColor,h,'g',prop};
	cb{5} = {@other_Class_LineColor,h,'b',prop};       cb{6} = {@other_Class_LineColor,h,'y',prop};
	cb{7} = {@other_Class_LineColor,h,'c',prop};       cb{8} = {@other_Class_LineColor,h,'m',prop};
	cb{9} = {@other_Class_LineColor,h,[],prop};

function other_Class_LineColor(obj,eventdata,h,cor,prop)
% If individual Lines were previously removed (by "Remove Line") h has invalid
% handles, so make sure all handles are valid
	h=h(ishandle(h));
	if ~isempty(cor)   
		if isempty(prop)	% line
			set(h,'color',cor);   refresh;
		else				% marker
			set(h,prop,cor);   refresh;
		end
	else
		c = uisetcolor;
		if (length(c) > 1)			% That is, if a color was selected
			if isempty(prop)
				set(h,'color',c);	refresh;
			else
				set(h,prop,c);		refresh;
			end
		end
	end
% -----------------------------------------------------------------------------------------

% -----------------------------------------------------------------------------------------
function cb = uictx_Class_LineStyle(h)
	% Set uicontext LineStyles in a class object hose handles are contained in h
	cb{1} = {@other_Class_LineStyle,h,'-'};      cb{2} = {@other_Class_LineStyle,h,'--'};
	cb{3} = {@other_Class_LineStyle,h,':'};      cb{4} = {@other_Class_LineStyle,h,'-.'};
	
	function other_Class_LineStyle(obj,eventdata,h,opt)
	% If individual Lines were previously removed (by "Remove Line") h has invalid
	% handles, so make sure all handles are valid
	h = h(ishandle(h));
	set(h,'LineStyle',opt);        refresh;
% -----------------------------------------------------------------------------------------

% -----------------------------------------------------------------------------------------
function set_greatCircle_uicontext(h)
% h is a handle to a great circle arc (in geog coords) object
	handles = guidata(h(1));	cmenuHand = uicontextmenu('Parent',handles.figure1);
	set(h, 'UIContextMenu', cmenuHand);
	cb_solid  = 'set(gco, ''LineStyle'', ''-''); refresh';   cb_dashed      = 'set(gco, ''LineStyle'', ''--''); refresh';
	cb_dotted = 'set(gco, ''LineStyle'', '':''); refresh';   cb_dash_dotted = 'set(gco, ''LineStyle'', ''-.''); refresh';
	
	uimenu(cmenuHand, 'Label', 'Delete', 'Callback', 'delete(gco)');
	uimenu(cmenuHand, 'Label', 'Save line', 'Callback', {@save_formated,h});
	uimenu(cmenuHand, 'Label', 'Line length', 'Callback', {@show_LineLength,[],'total'});
	uimenu(cmenuHand, 'Label', 'Line azimuth','Callback', @show_lineAzims)
	setLineWidth(uimenu(cmenuHand, 'Label', 'Line Width', 'Sep','on'), uictx_LineWidth(h))
	setLineStyle(uimenu(cmenuHand, 'Label', 'Line Style'),{cb_solid cb_dashed cb_dotted cb_dash_dotted})
	setLineColor(uimenu(cmenuHand, 'Label', 'Line Color'), uictx_color(h))

% -----------------------------------------------------------------------------------------
function set_circleGeo_uicontext(h)
% h is a handle to a circle (in geog coords) object
% NOTE: on 1-1-04 I finished a function called uicirclegeo that draws circles and provides
% controls to change various circle parameters. Because it makes extensive use of the lines
% userdata, the move_circle function of this file cannot be used, for it also changes userdata.
	tag = get(h,'Tag');
	handles = guidata(h(1));	cmenuHand = uicontextmenu('Parent',handles.figure1);
	set(h, 'UIContextMenu', cmenuHand);
	cb_solid  = 'set(gco, ''LineStyle'', ''-''); refresh';   cb_dashed      = 'set(gco, ''LineStyle'', ''--''); refresh';
	cb_dotted = 'set(gco, ''LineStyle'', '':''); refresh';   cb_dash_dotted = 'set(gco, ''LineStyle'', ''-.''); refresh';
	cb_roi = 'mirone(''DrawClosedPolygon_CB'',guidata(gcbo),gco)';
	
	uimenu(cmenuHand, 'Label', 'Delete', 'Callback', 'delete(gco)')
	uimenu(cmenuHand, 'Label', 'Save circle', 'Callback', {@save_formated,h})
	uimenu(cmenuHand, 'Label', 'Line length', 'Callback', {@show_LineLength,[]})
	if ~strcmp(tag,'CircleEuler')       % "Just" a regular geographical circle
		uimenu(cmenuHand, 'Label', 'Make me Donut', 'Sep','on', 'Callback', {@donutify, h(1), 'geog'})
		uimenu(cmenuHand, 'Label', 'Region-Of-Interest', 'Sep','on', 'Callback', cb_roi)
	else
		uimenu(cmenuHand, 'Label', 'Compute velocity', 'Sep','on', 'Callback', {@report_EulerVel,h})
		uimenu(cmenuHand, 'Label', 'Plot tangent arrow', 'Callback', {@report_EulerVel, h, 'tg'})
		uimenu(cmenuHand, 'Label', 'Plot normal arrow',  'Callback', {@report_EulerVel, h, 'no'})
	end
	setLineWidth(uimenu(cmenuHand, 'Label', 'Line Width', 'Sep','on'), uictx_LineWidth(h))
	setLineStyle(uimenu(cmenuHand, 'Label', 'Line Style'), {cb_solid cb_dashed cb_dotted cb_dash_dotted})
	setLineColor(uimenu(cmenuHand, 'Label', 'Line Color'), uictx_color(h))

% -----------------------------------------------------------------------------------------
function donutify(obj, evt, h, opt)
% Make a donut out of a circle (well, actually 2 circles). To be used on radial averages
	LLR = getappdata(h,'LonLatRad');
	prompt = {['Inner circle radius (< ' sprintf('%.6g)', LLR(3))]};
	rad = abs(str2double(inputdlg(prompt, 'Donut inner Rad', [1 35],{'0.0'})) );
	if (isempty(rad) || isnan(rad) || rad >= LLR(3)),	return,		end
	handles = guidata(h);		LS = get(h, 'LineStyle');	LW = get(h, 'LineWidth');	cC = get(h, 'color');
	x = get(h, 'XData');		y = get(h, 'YData');
	if (opt(1) == 'g')
		[latc, lonc] = circ_geo(LLR(2), LLR(1), rad, []);
		if (handles.geog == 2),		lonc = lonc + 360;		end		% Longitudes in the [0 360] interval
		tipo = 1;
	else		% Cartesian
		xx = linspace(-pi,pi,360);		yy = sin(xx);		xx = cos(xx);
		lonc = LLR(1) + rad * xx;		latc = LLR(2) + rad * yy;
		tipo = 0;
	end
	latc = latc(end:-1:1);	lonc = lonc(end:-1:1);		% Revert order as it will be in the inner side
	x = [x lonc x(1)];		y = [y latc y(1)];
	hP = patch('XData',x, 'YData',y, 'Parent',handles.axes1, 'LineWidth',LW, 'FaceColor','none',...
		'EdgeColor',cC, 'LineStyle',LS, 'Tag', 'Donut');
	setappdata(hP,  'donut', [LLR rad tipo])			% Store circle's center, outer, inner radius and geogacity
	cmenuHand = uicontextmenu('Parent',handles.figure1);
	set(hP, 'UIContextMenu', cmenuHand);
	uimenu(cmenuHand, 'Label', 'Delete', 'Callback', 'delete(gco)')
	uimenu(cmenuHand, 'Label', 'Radial average', 'Callback', 'grid_profiler(gcf, gco)', 'Sep', 'on')

% -----------------------------------------------------------------------------------------
function set_circleCart_uicontext(h)
% h is a handle to a circle (in cartesian coords) object
	handles = guidata(h(1));	cmenuHand = uicontextmenu('Parent',handles.figure1);
	set(h, 'UIContextMenu', cmenuHand);
	cb_solid  = 'set(gco, ''LineStyle'', ''-''); refresh';   cb_dashed      = 'set(gco, ''LineStyle'', ''--''); refresh';
	cb_dotted = 'set(gco, ''LineStyle'', '':''); refresh';   cb_dash_dotted = 'set(gco, ''LineStyle'', ''-.''); refresh';
	cb_roi = 'mirone(''DrawClosedPolygon_CB'',guidata(gcbo),gco)';
	
	uimenu(cmenuHand, 'Label', 'Delete', 'Callback', 'delete(gco)');
	uimenu(cmenuHand, 'Label', 'Save circle', 'Callback', {@save_formated,h});
	uimenu(cmenuHand, 'Label', 'Circle perimeter', 'Callback', {@show_LineLength,[]});
	uimenu(cmenuHand, 'Label', 'Move (interactive)', 'Callback', {@move_circle,h});
	uimenu(uimenu(cmenuHand, 'Label', 'Change'), 'Label', 'By coordinates', 'Callback', {@change_CircCenter1,h});
	uimenu(cmenuHand, 'Label', 'Make me Donut', 'Sep','on', 'Callback', {@donutify, h(1), 'cart'})
	uimenu(cmenuHand, 'Label', 'Region-Of-Interest', 'Sep','on', 'Callback', cb_roi);
	hp = getappdata(handles.figure1, 'ParentFig');
	if (~isempty(hp) && ishandle(hp) && ~isempty(strfind(get(handles.figure1,'Name'), 'spectrum')))
		uimenu(cmenuHand, 'Label', 'Low Pass FFT filter', 'Callback', 'mirone(''GridToolsSectrum_CB'',guidata(gcbo), ''lpass'', gco)');
		uimenu(cmenuHand, 'Label', 'High Pass FFT filter','Callback', 'mirone(''GridToolsSectrum_CB'',guidata(gcbo), ''hpass'', gco)');
	end
	setLineWidth(uimenu(cmenuHand, 'Label', 'Line Width', 'Sep','on'), uictx_LineWidth(h))
	setLineStyle(uimenu(cmenuHand, 'Label', 'Line Style'), {cb_solid cb_dashed cb_dotted cb_dash_dotted})
	setLineColor(uimenu(cmenuHand, 'Label', 'Line Color'), uictx_color(h))

% -----------------------------------------------------------------------------------------
function report_EulerVel(obj, evt, h, opt)
% 1- Report the volocity and azimuth of a movement around an Euler Pole
% 2- If OPT is used it must be a string with either 'tg' or 'no' that stands
%	 for plotting a tangent or normal arror to the circle at the click point.

	s = get(h,'Userdata');
	hAxes = get(h, 'Parent');
	pt = get(hAxes, 'CurrentPoint');
	[vel, azim] = compute_EulerVel(pt(1,2),pt(1,1),s.clat,s.clon,s.omega);

	if (nargin == 3)
		msg = [ sprintf('Pole name:  %s\n', s.plates), ...
				sprintf('Pole lon  = %3.3f\n',s.clon), ...
				sprintf('Pole lat   = %2.3f\n',s.clat), ...
				sprintf('Pole rate = %.3f\n',s.omega), ...
				sprintf('At point: \n'), ...
				sprintf('Lon = %3.3f\n',pt(1,1)), ...
				sprintf('Lat  = %2.3f\n',pt(1,2)), ...
				sprintf('Speed (cm/yr) =   %2.2f\n',vel), ...
				sprintf('Azimuth (CW from North) = %3.1f',azim)];
		msgbox(msg,'Euler velocity')
	else
		resp = inputdlg({'Enter scale in points per speed (cm/yr). Default 30 pt = 1 cm/yr'}, ...
			'Scale pt/velocity', [1 50], {'30'});
		resp = sscanf(resp{1},'%f');
		if isnan(resp),		return,		end
		
		hFig = get(hAxes, 'Parent');		handles = guidata(hFig);
		[hscale, vscale] = vectorFirstButtonDown(hFig, handles.axes1);
		if (~strcmp(opt,'tg')),		azim = azim + 90;		end
		mag = resp * ((hscale + vscale) / 2) * 111110;		% Length in meters, which is what vreckon wants

		[lat2,lon2] = vreckon(pt(1,2), pt(1,1), mag, [azim azim], 1);		% Get arrow tip point
		[xt, yt] = make_arrow([pt(1,1) lon2; pt(1,2) lat2], hscale, vscale, 10);
	
		hVec = patch('XData',xt, 'YData', yt, 'FaceColor',handles.DefLineColor,'EdgeColor', ...
			handles.DefLineColor,'LineWidth',0.5,'Tag','Arrow');
		ud.arrow_xy = [xt(:) yt(:)];		
		ud.vFac = 1.3;			ud.headLength = 10;
		ud.aspectRatio = 3/2;	ud.length = mag;
		ud.hscale = hscale;		ud.vscale = vscale;
		ud.anchors = [pt(1,1) lon2; pt(1,2) lat2]';		% Transpose because we store by columns
		ud.mag = mag;			ud.azim = azim;
		set(hVec, 'UserData', ud)
		set_vector_uicontext(hVec)
	end

% -----------------------------------------------------------------------------------------
function [vel, azim] = compute_EulerVel(alat,alon,plat,plon,omega, opt)
% alat & alon are the point coords. plat, plon & omega are the pole parameters (All in degrees)
% OPT, if given (anything will do), selects output as X,Y velocity components

	D2R = pi / 180;
	earth_rad = 6371e3;    % Earth radius in km
	plat = plat*D2R;      plon = plon*D2R;      omega = omega*D2R;
	alat = alat*D2R;      alon = alon*D2R;

	% Convert to geocentric
	ecc = 0.0818191908426215;		% WGS84
	alat = atan2((1-ecc^2) .* sin(alat), cos(alat));

	x = cos(plat).*sin(plon).*sin(alat) - cos(alat).*sin(alon).*sin(plat);    % East vel
	y = cos(alat).*cos(alon).*sin(plat) - cos(plat).*cos(plon).*sin(alat);    % North vel
	z = cos(plat).*cos(alat).*sin(alon-plon);
	vlon = -sin(alon).*x + cos(alon).*y;
	vlat = -sin(alat).*cos(alon).*x-sin(alat).*sin(alon).*y + cos(alat).*z;
	azim = 90 - atan2(vlat,vlon) / D2R;

	ind = (azim < 0);
	if (any(ind)),		azim(ind) = azim(ind) + 360;	end		% Give allways the result in the 0-360 range

	x = sin(alat).*sin(plat) + cos(alat).*cos(plat).*cos(plon-alon);
	delta = acos(x);
	vel = omega/1e+4 * earth_rad .* sin(delta);			% to give velocity in cm/Ma
	
	if (nargin == 6)
		azim = azim * D2R;			% Put it back in radians
		vel_ = vel .* sin(azim);	% X	(azim is angle to North)
		azim = vel .* cos(azim);	% Y
		vel = vel_;
	end

% -----------------------------------------------------------------------------------------
function set_vector_uicontext(h)
% h is a handle to a vector object
	h = h(1);
	handles = guidata(h);	cmenuHand = uicontextmenu('Parent',handles.figure1);
	set(h, 'UIContextMenu', cmenuHand);
	uimenu(cmenuHand, 'Label', 'Delete', 'Callback', 'delete(gco)');
	uimenu(cmenuHand, 'Label', 'Save line', 'Callback', {@save_formated, h});
	uimenu(cmenuHand, 'Label', 'Copy', 'Callback', {@copy_line_object, handles.figure1, handles.axes1});
 	uimenu(cmenuHand, 'Label', 'Edit (redraw)', 'Callback', @arrowRedraw);
	uimenu(cmenuHand, 'Label', 'Reshape (head)','Callback', {@arrow_shape,h});
	if (handles.geog)		% No solution yet for cartesian arrows
		uimenu(cmenuHand, 'Label', 'Copy (mirror)', 'Callback', @mirror_arrow);
	end
	x = get(h, 'xdata');	y = get(h, 'ydata');
	L = [x(3) x(6); y(3) y(6)]';	% To compute arrow length (head included)
	uimenu(cmenuHand, 'Label', 'Arrow length (with head)', 'Callback', {@show_LineLength, L}, 'Sep','on')
	L = [x(7) x(8); y(7) y(8)]';	% To compute arrow azimuth
	uimenu(cmenuHand, 'Label', 'Arrow azimuth','Callback', {@show_lineAzims, L})
	setLineColor(uimenu(cmenuHand, 'Label', 'Line Color','Sep','on'), uictx_color(h,'EdgeColor'))
	item2 = uimenu(cmenuHand, 'Label','Fill Color');
	setLineColor( item2, uictx_color(h, 'facecolor') )
	uimenu(item2, 'Label', 'None', 'Sep','on', 'Callback', 'set(gco, ''FaceColor'', ''none'');refresh');
	uimenu(cmenuHand, 'Label', 'Transparency', 'Callback', @set_transparency);
	set(h, 'ButtonDownFcn', 'move_obj(1)')

% -----------------------------------------------------------------------------------------
function arrowRedraw(obj, evt)
	% Edit a previously existing arrow
	hVec = gco;
	hAxes = get(hVec,'parent');
	hFig = get(hAxes,'parent');
	state = uisuspend_j(hFig);		% Remember initial figure state
	set(hFig,'Pointer', 'crosshair');
	hVec(2) = line('XData', [], 'YData', [],'LineWidth',0.5,'Tag','Arrow');
	ud = get(hVec(1),'UserData');
	vectorFirstButtonDown(hFig, hAxes, hVec, state, ud.anchors, ud.headLength, ud.vFac, ud.aspectRatio)

% -----------------------------------------------------------------------------------------
function mirror_arrow(obj, evt)
	% Make a copy of the arrow whose handle is H, but mirrored. That is pointing in the oposite sense
	h = gco;
	handles = guidata(h);
	ud = get(h, 'UserData');
	x1 = (ud.arrow_xy(end-2,1) + ud.arrow_xy(end-1,1)) / 2;
	y1 = (ud.arrow_xy(end-2,2) + ud.arrow_xy(end-1,2)) / 2;
	[y2, x2] = vreckon(y1, x1, ud.mag, [ud.azim ud.azim]+180, 1);
	[xt, yt] = make_arrow([x1 x2; y1 y2], ud.hscale, ud.vscale, ud.headLength, ud.vFac);
	hVec = patch('XData',xt, 'YData', yt,'FaceColor',handles.DefLineColor,'EdgeColor', ...
		handles.DefLineColor,'LineWidth',0.5,'Tag','Arrow');
	ud.arrow_xy = [xt(:) yt(:)];		ud.azim = rem(ud.azim + 180, 360);
	ud.anchor = ud.anchor(end:-1:1,:);	% Other ud fields don't need updating
	set(hVec, 'UserData', ud)
	set_vector_uicontext(hVec)

% -----------------------------------------------------------------------------------------
function show_swhatRatio(obj,evt,h)
    msgbox(sprintf('Swath Ratio for this track is: %g',getappdata(h,'swathRatio')),'')

% -----------------------------------------------------------------------------------------
% function stack_profiles(obj,evt)
% % ...
% 	h = gco;
% 	prompt = {'Cross profile length (km):','Distance between cross profiles:'};
% 	resp = inputdlg(prompt, 'Stacked profiles', [1 30], {'',''});
% 	if (isempty(resp)),		return,		end
% 	len = str2double(resp{1});	sep = str2double(resp{2});
% 	if (isnan(len) || isnan(sep))
% 		errordlg('#@&%$#!&%$$#&?', 'Error'),	return
% 	end
% 	
% 	handles = guidata(h);
% 	[X, Y, Z] = load_grd(handles);
% 	G = gmt('wrapgrid', Z, handles.head);
% 	x = get(h, 'XData');	y = get(h, 'YData');
% 	ds = min(handles.head(8:9)) * pi / 180 * handles.DefineEllipsoide(1) / 1000;	% In approx km
% 	opt_C = sprintf('-C%sk/%g/%s', resp{1}, ds, resp{2});
% 	[profiles, stack] = gmtmex(['grdtrack -G ' opt_C ' -Sm+s+a --FORMAT_GEO_OUT=+D'], [x(:) y(:)], G);
% 	% Show the stacked profiles
% 	figure(2);	hold on
% 	for (k = 1:length(profiles))
% 		p = polyfit(profiles(k).data(:,3),profiles(k).data(:,5),1);		y_fit = polyval(p,profiles(k).data(:,3));
% 		plot(profiles(k).data(:,3), profiles(k).data(:,5)-y_fit)
% 	end
% 	p = polyfit(stack.data(:,1), stack.data(:,2),1);		y_fit = polyval(p,stack.data(:,1));
% 	plot(stack.data(:,1), stack.data(:,2)-y_fit, 'LineWidth', 3)

% -----------------------------------------------------------------------------------------
function chessify(obj, evt, h)
% Create a chess-board grid
	prompt = {'Enter DX patch size:','Enter DY patch size:','One Z value','Other Z value'};
	resp = inputdlg(prompt, 'Chessboard size and values', [1 30], {'','','1','-1'});
	if (isempty(resp)),		return,		end
	handles = guidata(h);
	dx = str2double(resp{1});	dy = str2double(resp{2});
	up_down = [str2double(resp{3}) 0 str2double(resp{4})];
	x = get(h,'XData');			y = get(h,'YData');
	dx_r = max(x) - min(x);		dy_r = max(y) - min(y);
	nCols = ceil(dx_r / dx);	nRows = ceil(dy_r / dy);
	for (row = 1:nRows)
		for (col = 1:nCols)
			xr = x(1) + [(col - 1) col] * dx;
			yr = y(1) + [(row - 1) row] * dy;
			rec = [xr(1) yr(1); xr(1) yr(2); xr(2) yr(2); xr(2) yr(1); xr(1) yr(1)];
			mirone('ImageCrop_CB', handles, rec, 'SetConst', up_down((-1)^(row+col-1) + 2))
		end
	end

% -----------------------------------------------------------------------------------------
function show_Area(obj, evt, h)
% Compute area under line and insult the user if the line is not closed
% NOTE that H is optional. Use only when want to make sure that this fun
% uses that handle (does not work with copyied objects)

	if (nargin == 3)
		if (size(h,1) >= 2 && size(h,2) == 2)
			x = h(:,1);				y = h(:,2);
			handles = guidata(get(0,'CurrentFigure'));
		elseif (ishandle(h))
			x = get(h,'XData');		y = get(h,'YData');
			handles = guidata(h);
		end
	elseif (nargin == 2 || isempty(h) || length(h) > 1)
		h = gco;
		x = get(h,'XData');			y = get(h,'YData');
		handles = guidata(h);
	else
		handles = guidata(get(0,'CurrentFigure'));
	end

	% Contour lines for example have NaNs and not at the same x,y positions (???)
	ix = isnan(x);
	x(ix) = [];				y(ix) = [];
	iy = isnan(y);
	x(iy) = [];				y(iy) = [];
	if ~( (x(1) == x(end)) && (y(1) == y(end)) )
		msg{1} = 'This is not a closed line. Therefore the result is probably ...';
	else
		msg{1} = '';
	end
	if (handles.geog)
		area = area_geo(y,x,handles.DefineEllipsoide);    % Area is reported on the default ellipsoide
		str_units = 'm^2';		frmt = '%.2f';
		if (handles.DefineMeasureUnit(1) == 'k' && area > 1e4)
			area = area * 1e-6;		str_units = 'km^2';		frmt = '%.3f';
		end
		msg{2} = sprintf(['Area = ' frmt ' ' str_units], area);
		msgbox(msg,'Area')
	else
		area = polyarea(x,y);   % Area is reported in map user unites
		msg{2} = ['Area = ' sprintf('%g',area) ' map units ^2'];
		msgbox(msg,'Area')
	end

% -----------------------------------------------------------------------------------------
function ll = show_LineLength_XY(obj, evt, h)
% Print or return the length per X & Y component.
% If output, return a structure ll.len_x and ll.len_y, with the X & Y length of each line segment
% This function is meant to use on lines plotted into a profile ('ecran()') figure where X and
% Y components have different units (y = f(x)) and therefore only per component info has meaning.

	if (nargin == 2 || isempty(h)),		h = gco;	end
	x = get(h,'XData');		y = get(h,'YData');
	dx = diff(x);			dy = diff(y);
	if (nargout == 0)
		if (numel(dx) > 1)
			msg = cell(1, numel(dx) + 1);
			for (i = 1:numel(dx))
				msg{i} = sprintf('Length%d_X = %g;  Length%d_Y = %g   map units',i, dx, i, dy);
			end
			msg{i+1} = sprintf('Total (map units): length_X = %g;\nlength_y = %g',sum(dx), sum(dy));
		else
			msg = {sprintf('Map Units\nlength_X = %g;\nlength_Y = %g',dx, dy)};
		end
		msgbox(msg,'Line(s) length')
	else
		ll.len_x = dx;		ll.len_y = dy;
	end

% -----------------------------------------------------------------------------------------
function ll = show_LineLength(obj, evt, h, opt)
% Line length (perimeter if it is a closed polyline). If output argument, return a structure
% ll.len, ll.seg_len and ll.type, where "len" is total line length, "seg_len" is the length of
% each segment (for when there are more than one) and "type" is either 'geog' or 'cart'.
% For polylines ll.len contains only the total length.
% 22-09-04  Added OPT option. If it exists, report only total length (for nargout == 0)
%	2012    OPT can also be used to transmit the coordinates type in which report the result.
%           E.G. 'n', 'k', 'm' or 'u'
% 22-10-05  H is now only to be used if we want to specificaly use that handle. Otherwise use []
%           to fish it with gco (MUST use this form to work with copied objects)
% 16-08-07  H can contain a Mx2 column vector with the line vertices.

	n_args = nargin;
	if (n_args == 2 || isempty(h)),		h = gco;	end
	if (n_args <= 3)
		opt = [];
	elseif (numel(opt) == 1 && isa(opt, 'char'))
		measureUnit = opt;
	end

	if (n_args == 3 || (numel(h) > 1 && ~ishandle(h(1))) )	% Test only h(1) to reduce the risk of h(?) = 0 and ... fcked
		if (size(h,1) >= 2 && size(h,2) == 2)
			x = h(:,1);     y = h(:,2);
			handles = guidata(get(0,'CurrentFigure'));
		elseif (ishandle(h))
			x = get(h,'XData');    y = get(h,'YData');
			handles = guidata(h);
		end
	elseif ((n_args == 2 || n_args == 4 || length(h) > 1) && ~isempty(h) && ishandle(h(1)))
        x = get(h,'XData');    y = get(h,'YData');
		handles = guidata(h);
	else
		errordlg('Unknown case in show_LineLength()','error'),	return
	end

	if (isempty(opt) || strcmp(opt, 'nikles') || strcmp(opt, 'total'))	% That call with nikles was shity bad idea
		measureUnit = handles.DefineMeasureUnit(1);
	end

	% Contour lines for example have NaNs and not at the same x,y positions (???)
	ix = isnan(x);      x(ix) = [];     y(ix) = [];
	iy = isnan(y);      x(iy) = [];     y(iy) = [];
	if (handles.geog && measureUnit ~= 'u')
		lat_i = y(1:end-1);   lat_f = y(2:end);     clear y;
		lon_i = x(1:end-1);   lon_f = x(2:end);     clear x;
		tmp = vdist(lat_i,lon_i,lat_f,lon_f,handles.DefineEllipsoide([1 3]));

		switch measureUnit
			case 'n'        % Nautical miles
				scale = 1852;   str_unit = ' NM';
			case 'k'        % Kilometers
				scale = 1000;   str_unit = ' kilometers';
			case {'m','u'}  % Meters or user unites
				scale = 1;   str_unit = ' meters(?)';
		end
		total_len = sum(tmp) / scale;
		if (nargout == 0 && isempty(opt))
			len_i = tmp / scale;
			if (numel(tmp) <= 20)
				msg = cell(1, numel(tmp) + 1);
				for i = 1:numel(tmp)
					msg{i} = [sprintf('Length%d',i) '  =  ' sprintf('%.5f',len_i(i)) str_unit];
				end
				msg{i+1} = sprintf('Total length = %.5f %s',total_len, str_unit);
			else
				msg = {sprintf('Total length = %.5f %s',total_len, str_unit)};
			end
			msgbox(msg,'Line(s) length')
		elseif (nargout == 0 && ~isempty(opt))
			msgbox(sprintf('Total length = %.5f %s',total_len, str_unit),'Line length')
		else        % Should we also out output also the partial lengths?
			ll.len = total_len;   ll.seg_len = tmp / scale;	ll.type = 'geog';
		end
	else
		dx = diff(x);		dy = diff(y);
		len_i = sqrt(dx.*dx + dy.*dy);
		total_len = sum(len_i);
		if (nargout == 0 && isempty(opt))
			if (numel(dx) <= 200)
				msg = cell(1, numel(dx) + 1);
				for i = 1:numel(dx)
					msg{i} = [sprintf('Length%d',i) '  =  ' sprintf('%.5f',len_i(i)) ' map units'];
				end
				msg{i+1} = sprintf('Total length = %.5f  map units',total_len);
			else
				msg = {sprintf('Total length = %.5f  map units',total_len)};
			end
			msgbox(msg,'Line(s) length')
		elseif (nargout == 0 && ~isempty(opt))
			msgbox([sprintf('Total length = %.5f',total_len) ' map units'],'Line length')
		else		% The same question as in the geog case
			ll.len = total_len;		ll.seg_len = len_i;		ll.type = 'cart';
		end
	end

% -----------------------------------------------------------------------------------------
function show_AllTrackLength(obj,eventdata)
% Compute the length of all MB tracks present in the figure
	ALLlineHand = findobj(get(gca,'Child'),'Type','line');
	len = 0;
	for i = 1:length(ALLlineHand)
		tag = get(ALLlineHand(i),'Tag');
		if ~isempty(strfind(tag,'MBtrack'))       % case of a MBtrack line
			tmp = show_LineLength(obj,eventdata,ALLlineHand(i));        
			len = len + tmp.len;
		end
	end
	if (len > 0)
        msgbox(sprintf('Total tracks length = %g NM',len))
	end

% -----------------------------------------------------------------------------------------
function azim = show_lineAzims(obj, evt, h)
% Works either in geog or cart coordinates. Otherwise the result is a non-sense
% If an output argument is requested , return a structure azim.az and azim.type, where "az" is line
% azimuth and "type" is either 'geog' or 'cart'.
% 22-10-05  H is now only to be used if we want to specifically use that handle. Otherwise use
%           either [] or don't pass the H argument to fish it with gco (MUST use this form to work with copied objects)
% 16-08-07  H can contain a Mx2 column vector with the line vertices.
% 28-08-08  If > 10 azimuths & ~nargout send the result to "ecran"

	if (nargin == 3)
		if (size(h,1) >= 2 && size(h,2) == 2)   
			x = h(:,1);     y = h(:,2);
			handles = guidata(get(0,'CurrentFigure'));
		elseif (ishandle(h))
			x = get(h,'XData');    y = get(h,'YData');
			handles = guidata(h);
		end
	elseif (nargin == 2 || isempty(h) || numel(h) > 1)
		h = gco;  
		az = getappdata(h, 'Azim');
		if (~isempty(az))		% Report azim of great-circle arc and return
			msgbox(sprintf('Azimuth = %.2f', az),'Azimuth'),	return
		end
		x = get(h,'XData');    y = get(h,'YData');
		handles = guidata(h);
	end
	
	if (handles.geog)
		% I should do this instead
		%[lix, az] = vdist(y(1:end-1), x(1:end-1), y(2:end), x(2:end),handles.DefineEllipsoide([1 3]));
        az = azimuth_geo(y(1:end-1), x(1:end-1), y(2:end), x(2:end));
        azim.type = 'geog';                 % Even if it is never used
	else
        dx = diff(x);   dy = diff(y);
        angs = atan2(dy,dx) * 180/pi;       % and convert to degrees
        hFig = get(0,'CurrentFigure');      hAxes = get(hFig,'CurrentAxes');
        if(strcmp(get(hAxes,'YDir'),'reverse')),    angs = -angs;   end
        az = (90 - angs);                   % convert to azim (cw from north)
        ind = find(az < 0);
        az(ind) = 360 + az(ind);
        azim.type = 'cart';                 % Even if it is never used
	end

	if (nargout == 0)
		if (numel(az) <= 10)
			msg = cell(numel(az),1);
			for (i = 1:numel(az))
				msg{i} = ['Azimuth' sprintf('%d',i) '  =  ' sprintf('%3.1f',az(i)) '  degrees'];
			end
			msg{end+1} = '';
			id = (az > 270);    az(id) = az(id) - 360;
			az_mean = mean(az);
			msg{end+1} = sprintf('Mean azimuth = %.1f  degrees',az_mean);
			msgbox(msg,'Line(s) Azimuth')
		else
			ecran(handles, x(2:end), y(2:end), az, 'Polyline azimuths (deg)')
		end
	else
		azim.az = az;
	end
	refresh	

% -----------------------------------------------------------------------------------------
function [xx, yy] = smooth_line(obj, evt, h)
% Smooth data by 2-D spline interpolation
% Core code from Duane Hanselman (BSD License) http://www.mathworks.com/matlabcentral/fileexchange/38862

	x = get(h, 'XData');	y = get(h, 'YData');

	if (size(x,1) > 1)
		x = x(:)';		y = y(:)';
	end
	isopen = ~((x(1) == x(end)) && (y(1) == y(end)));
	if ~isopen  % wrap data so closed curve is smooth at joint
		x = [x(end-1) x x(2)];
		y = [y(end-1) y y(2)];
	end

	t = [0 cumsum(sqrt(diff(x).^2 + diff(y).^2))];	% get path length to create independent variable

	% place interpolation points in between those in t    
	n  = max(2,ceil(20/sqrt(numel(t))));
	ti = repmat(t(1:end-1),n,1);
	d  = repmat((0:n-1)'/n,1,length(x)-1);
	dt = repmat(diff(t),n,1);
	ti = ti + d .* dt;
	ti = [ti(:); t(end)];							% independent variable interpolation points

	% computer new contour points from spline fit
	xi = akimaspline(t,x,ti);
	yi = akimaspline(t,y,ti);
	if ~isopen   % take out redundant data if curve was closed
		xi = xi(n+1:end-n);
		yi = yi(n+1:end-n);
	end

	set(h, 'XData',xi, 'YData',yi)
	if (nargout)
		xx = xi;		yy = yi;
	end

% -----------------------------------------------------------------------------------------
function set_bar_uicontext(h)
	% Set uicontexts for the bars in a multibeam track. h is the handle to the bar objects
	handles = guidata(h(1));	cmenuHand = uicontextmenu('Parent',handles.figure1);
	set(h, 'UIContextMenu', cmenuHand);
	cb_LineWidth = uictx_LineWidth(h);      % there are 5 cb_LineWidth outputs, but I only use 5 here
	cb7 = 'set(gco, ''LineWidth'', 10); refresh';
	cb10 = 'set(gco, ''LineStyle'', ''-''); refresh';   cb11 = 'set(gco, ''LineStyle'', ''--''); refresh';
	cb12 = 'set(gco, ''LineStyle'', '':''); refresh';   cb13 = 'set(gco, ''LineStyle'', ''-.''); refresh';
	cb_color = uictx_color(h);      % there are 9 cb_color outputs
	item1 = uimenu(cmenuHand, 'Label', 'Line Width');
	uimenu(item1, 'Label', '1       pt', 'Callback', cb_LineWidth{1});
	uimenu(item1, 'Label', '2       pt', 'Callback', cb_LineWidth{2});
	uimenu(item1, 'Label', '3       pt', 'Callback', cb_LineWidth{3});
	uimenu(item1, 'Label', '4       pt', 'Callback', cb_LineWidth{4});
	uimenu(item1, 'Label', '5       pt', 'Callback', cb7);
	uimenu(item1, 'Label', 'Other...', 'Callback', cb_LineWidth{5});
	item_ls = uimenu(cmenuHand, 'Label', 'Line Style');
	setLineStyle(item_ls,{cb10 cb11 cb12 cb13})
	item_lc = uimenu(cmenuHand, 'Label', 'Line Color');
	setLineColor(item_lc,cb_color)

% -----------------------------------------------------------------------------------------
function cb = uictx_color(h, opt)
% Set uicontext colors in object whose handle is gco (or h for "other color")
% If opt is not given opt = 'Color' is assumed
	if (nargin == 1),   opt = [];   end
	if (~isempty(opt) && ischar(opt))
		c_type = opt;
	else
		c_type = 'Color';
	end
	cb{1} = ['set(gco,''' c_type ''',''k'');refresh'];       cb{2} = ['set(gco,''' c_type ''',''w'');refresh'];
	cb{3} = ['set(gco,''' c_type ''',''r'');refresh'];       cb{4} = ['set(gco,''' c_type ''',''g'');refresh'];
	cb{5} = ['set(gco,''' c_type ''',''b'');refresh'];       cb{6} = ['set(gco,''' c_type ''',''y'');refresh'];
	cb{7} = ['set(gco,''' c_type ''',''c'');refresh'];       cb{8} = ['set(gco,''' c_type ''',''m'');refresh'];
	cb{9} = {@other_color,h,opt};
	if (strcmp(opt, 'facecolor') && get(h, 'FaceAlpha') <= 0.005)	% Because of the R2015 shit with patches
		ind = strfind(cb{1}, ')');
		for (k = 1:8)
			cb{k} = [cb{k}(1:ind-1) ',''FaceAlpha'',1);refresh'];	% Have to set opacity so color can be seen
		end
	end

	function other_color(obj, evt, h, opt)
	if (nargin == 3),   opt = [];   end
	c = uisetcolor;
	if (isequal(c, 0)),		return,		end			% User gave up
	if (~isempty(opt) && ischar(opt))
		set(h, opt, c);
		if (strcmp(opt, 'facecolor') && get(h, 'FaceAlpha') <= 0.005)
			set(h, 'FaceAlpha', 1)
		end
		refresh
	else
		set(h,'color',c);   refresh
	end
% -----------------------------------------------------------------------------------------

% -----------------------------------------------------------------------------------------
function cb = uictx_LineWidth(h)
% Set uicontext colors in object hose handle is gco (or h for "other color")
	cb{1} = 'set(gco, ''LineWidth'', 1);refresh';		cb{2} = 'set(gco, ''LineWidth'', 2);refresh';
	cb{3} = 'set(gco, ''LineWidth'', 3);refresh';		cb{4} = 'set(gco, ''LineWidth'', 4);refresh';
	cb{5} = {@other_LineWidth,h};

function other_LineWidth(obj,evt,h)
	lw = get(h, 'LineWidth');
	resp  = inputdlg({'Enter new line width (pt)'}, 'Line width', [1 30], {num2str(lw)});
	if isempty(resp),	return,		end
	set(h,'LineWidth',str2double(resp));        refresh
% -----------------------------------------------------------------------------------------

% -----------------------------------------------------------------------------------------
function hVec = DrawVector
	hFig = get(0,'CurrentFigure');		handles = guidata(hFig);
	hVec(1) = patch('XData',[], 'YData', [],'FaceColor',handles.DefLineColor,'EdgeColor', ...
	                handles.DefLineColor,'LineWidth',0.5,'Tag','Arrow');
	hVec(2) = line('XData', [], 'YData', [],'Color',handles.DefLineColor,'LineWidth',0.5,'Tag','Arrow');
	state = uisuspend_j(hFig);			% Remember initial figure state
	set(hFig,'Pointer', 'crosshair');
	w = waitforbuttonpress;
	if (w == 0)							% A mouse click
		vectorFirstButtonDown(hFig, handles.axes1, hVec, state)
	else
		set(hFig,'Pointer', 'arrow');	hVec = [];
	end

function [hs, vs] = vectorFirstButtonDown(hFig, hAxes, h, state, anchor, ah, vFac, aspect)
% Outputs are used by the report_EulerVel (when PLOTing) function
% ANCHOR, if used, contains the arrow anchor point coords (used by the edit mode)
% When output is requested, minimum args is ([], HAXES).

	if (nargin >= 5),		pt = anchor;
	else,					pt = get(hAxes, 'CurrentPoint');
	end
	if (nargin ~= 8)
		ah = 12;	vFac = 1.3;		aspect = 3/2;
	end
 	axLims = getappdata(hAxes,'ThisImageLims');
	% create a conversion from data to points for the current axis
	oldUnits = get(hAxes,'Units');			set(hAxes,'Units','points');
	Pos = get(hAxes,'Position');			set(hAxes,'Units',oldUnits);
	vscale = 1/Pos(4) * diff(axLims(3:4));	hscale = 1/Pos(3) * diff(axLims(1:2));
	vscale = (vscale + hscale) / 2;			hscale = vscale;	% For not having a head direction dependency
	DAR = get(hAxes, 'DataAspectRatio');
	if (DAR(1) == 1 && DAR(1) ~= DAR(2))	% To account for the "Scale geog images at mean lat" effect
		vscale = vscale * DAR(2);		hscale = hscale * DAR(1);
	end
	if (~nargout)
		setappdata(h(1),'anchor', [pt(1,1) pt(1,2)])	% To eventual use in/if edit mode
		set(hFig,'WindowButtonMotionFcn',{@wbm_vector,[pt(1,1) pt(1,2)],h,hAxes,hscale,vscale, ah, vFac, aspect}, ...
			'WindowButtonDownFcn',{@wbd_vector, h, state, hscale, vscale});
	else
		hs = hscale;	vs = vscale;
	end

function wbm_vector(obj, evt, origin, h, hAxes, hscale, vscale, ah, vFac, aspect)
	pt = get(hAxes, 'CurrentPoint');
	x  = [origin(1) pt(1,1)];   y = [origin(2) pt(1,2)];
 	set(h(2),'XData',x, 'YData',y)
	[xt, yt] = make_arrow(h(2) , hscale, vscale, ah, vFac, aspect);
 	set(h(1),'XData',xt, 'YData',yt);

function wbd_vector(obj, evt, h, state, hscale, vscale)
% End interactive edition
	uirestore_j(state, 'nochildren');		% Restore the figure's initial state
	x  = get(h(2), 'XData');	y  = get(h(2), 'YData');
	xx = get(h(1), 'XData');	yy = get(h(1), 'YData');
	ud_old = get(h(1), 'UserData');
	ud.arrow_xy = [xx(:) yy(:)];
	ud.anchors = [x(:) y(:)];				% Extremities, usefull for arrow reconstruction
	ud.hscale = hscale;			ud.vscale = vscale;
	ll = show_lineAzims([], evt, h(1));		% Compute azimuth (geog or cartesian)
	ud.azim = ll.az(end);
	ll = show_LineLength([], evt, [x(:) y(:)], 'm');
	ud.mag = ll.len;
	if (isempty(ud_old))					% First time usage. Need to add these to UD
		ud.vFac = 1.3;				ud.headLength = 12;
		ud.aspectRatio = 3/2;				% Header lenght/width aspect ratio
	else
		ud.vFac = ud_old.vFac;		ud.headLength = ud_old.headLength;
		ud.aspectRatio = ud_old.aspectRatio;
	end
	set(h(1), 'UserData', ud)
	delete(h(2));							% We don't need this (support) line anymore
	set_vector_uicontext(h(1))
% -----------------------------------------------------------------------------------------

% -----------------------------------------------------------------------------------------
function h_gcirc = DrawGreatCircle
    hFig = get(0,'CurrentFigure');		handles = guidata(hFig);
    h_gcirc = line('XData', [], 'YData', [],'Color',handles.DefLineColor,'LineWidth',handles.DefLineThick);
	state = uisuspend_j(hFig);		% Remember initial figure state
	set(hFig,'Pointer', 'crosshair');	% to avoid the compiler BUG
	w = waitforbuttonpress;
	if w == 0							% A mouse click
        gcircFirstButtonDown(hFig,h_gcirc,state)
		set_greatCircle_uicontext(h_gcirc)
	else
        set(hFig,'Pointer', 'arrow');
        h_gcirc = [];
	end
%---------------
function gcircFirstButtonDown(hFig,h,state)
	hAxes = get(hFig,'CurrentAxes');	pt = get(hAxes, 'CurrentPoint');
	set(hFig,'WindowButtonMotionFcn',{@wbm_gcircle,[pt(1,1) pt(1,2)],h,hFig,hAxes},'WindowButtonDownFcn',{@wbd_gcircle,h,state});
%---------------
function wbm_gcircle(obj,eventdata,first_pt,h,hFig,hAxes)
	pt = get(hAxes, 'CurrentPoint');
	[x,y] = gcirc(first_pt(1),first_pt(2),pt(1,1),pt(1,2));
	% Find the eventual Date line discontinuity and insert a NaN on it
	% ind = find(abs(diff(x)) > 100);		% 100 is good enough
	% if (~isempty(ind))
	%     if (length(ind) == 2)
	%         x = [x(1:ind(1)) NaN x(ind(1)+1:ind(2)) NaN x(ind(2)+1:end)];
	%         y = [y(1:ind(1)) NaN y(ind(1)+1:ind(2)) NaN y(ind(2)+1:end)];
	%     elseif (length(ind) == 1)
	%         x = [x(1:ind) NaN x(ind+1:end)];		y = [y(1:ind) NaN y(ind+1:end)];
	%     end
	% end
	set(h, 'XData', x, 'YData', y);
%---------------
function wbd_gcircle(obj,eventdata,h,state)
	x = get(h, 'XData');			y = get(h, 'YData');
	az = azimuth_geo(y(1), x(1), y(end), x(end));
	setappdata(h,'Azim',az);		set(h,'Tag','GreatCircle')
	uirestore_j(state, 'nochildren');	% Restore the figure's initial state
% -----------------------------------------------------------------------------------------

% -----------------------------------------------------------------------------------------
function h_circ = DrawCartesianCircle
% THIS IS ONLY USED NOW WITH CARTESIAN CIRCLES
% Given one more compiler BUG, (WindowButtonDownFcn cannot be redefined)
% I found the following workaround.
	hFig = get(0,'CurrentFigure');          handles = guidata(hFig);
	h_circ = line('XData', [], 'YData', [],'Color',handles.DefLineColor,'LineWidth',handles.DefLineThick);
	%set(hFig,'WindowButtonDownFcn',{@circFirstButtonDown,h_circ}, 'Pointer', 'crosshair');
	state = uisuspend_j(hFig);		% Remember initial figure state
	set(hFig,'Pointer', 'crosshair');	% to avoid the compiler BUG
	w = waitforbuttonpress;
	if w == 0       % A mouse click
        circFirstButtonDown(h_circ,state)
		set_circleCart_uicontext(h_circ)
	else
        set(get(0,'CurrentFigure'),'Pointer', 'arrow');
        h_circ = [];
	end

%---------------
function circFirstButtonDown(h,state)
    x = linspace(-pi,pi,360);
	setappdata(h,'X',cos(x));       setappdata(h,'Y',sin(x))    % Save unit circle coords
    hFig = get(0,'CurrentFigure');  hAxes = get(hFig,'CurrentAxes');
	pt = get(hAxes, 'CurrentPoint');
	set(hFig,'WindowButtonMotionFcn',{@wbm_circle,[pt(1,1) pt(1,2)],h,hAxes}, ...
        'WindowButtonDownFcn',{@wbd_circle,h,state});

%---------------
function wbm_circle(obj,eventdata,center,h,hAxes)
	pt = get(hAxes, 'CurrentPoint');
	rad = sqrt( (pt(1,1)-center(1))^2 + (pt(1,2)-center(2))^2);
	%[y,x] = circ_geo(center(2),center(1),rad);
	x = getappdata(h,'X');          y = getappdata(h,'Y');
	x = center(1) + rad * x;        y = center(2) + rad * y;
	set(h, 'XData', x, 'YData', y,'Userdata',[center(1) center(2) rad]);

%---------------
function wbd_circle(obj,eventdata,h,state)
	lon_lat_rad = get(h,'UserData');	setappdata(h,'LonLatRad',lon_lat_rad)   % save this in appdata
	set(h,'Tag','circleCart')
    rmappdata(h,'X');			rmappdata(h,'Y');
	set(h,'UserData',[])				% Clean UserData so that it doesn't risk to crash write_fleder
	set(h, 'ButtonDownFcn', 'move_obj(1)')
	uirestore_j(state, 'nochildren');	% Restore the figure's initial state
% -----------------------------------------------------------------------------------------

% -----------------------------------------------------------------------------------------
function h_circ = draw_circleEulerPole(lon,lat)
% Draw a circle (or arc of a circle) about the Euler Pole (or any other origin)
% See notes above for the reason why waitforbuttonpress is used.
	h_circ = line('XData', [], 'YData', []);
	hFig = get(0,'CurrentFigure');		hAxes = get(hFig,'CurrentAxes');
	set(hFig,'Pointer', 'crosshair');
	w = waitforbuttonpress;
	if w == 0       % A mouse click
        set(hFig,'WindowButtonMotionFcn',{@wbm_circle,[lon lat],h_circ,hAxes},'WindowButtonDownFcn',{@wbd_circle,h_circ});
		set_circleGeo_uicontext(h_circ)
	else
        set(hFig,'Pointer', 'arrow');
        h_circ = [];
	end

% -----------------------------------------------------------------------------------------
function move_circle(obj,eventdata,h)
% ONLY FOR CARTESIAN CIRCLES.
	hFig = get(0,'CurrentFigure');		hAxes = get(hFig,'CurrentAxes');
	state = uisuspend_j(hFig);		% Remember initial figure state
	np = numel(get(h,'XData'));			x = linspace(-pi,pi,np);
	setappdata(h,'X',cos(x));			setappdata(h,'Y',sin(x))	% Save unit circle coords
	center = getappdata(h,'LonLatRad');
	set(hFig,'WindowButtonMotionFcn',{@wbm_MoveCircle,h,center,hAxes},...   
		'WindowButtonDownFcn',{@wbd_MoveCircle,h,state,hAxes},'Pointer', 'crosshair');

function wbm_MoveCircle(obj,eventdata,h,center,hAxes)
	pt = get(hAxes, 'CurrentPoint');
	x = getappdata(h,'X');				y = getappdata(h,'Y');
	x = pt(1,1) + center(3)*x;			y = pt(1,2) + center(3)*y;
	set(h, 'XData', x, 'YData', y, 'Userdata',[pt(1,1) pt(1,2) center(3)]);

function wbd_MoveCircle(obj,eventdata,h,state,hAxes)
	% check if x,y is inside of axis
	pt = get(hAxes, 'CurrentPoint');	x = pt(1,1);	y = pt(1,2);
	x_lim = get(hAxes,'xlim');			y_lim = get(hAxes,'ylim');
	if (x < x_lim(1)) || (x > x_lim(2)) || (y < y_lim(1)) || (y > y_lim(2)),	return,	end
	lon_lat_rad = get(h,'UserData');	setappdata(h,'LonLatRad',lon_lat_rad)	% save this in appdata
	set(h,'UserData',[])			% Clean UserData so that it doesn't risk to crash write_fleder
    rmappdata(h,'X');				rmappdata(h,'Y');
	uirestore_j(state, 'nochildren');	% Restore the figure's initial state

% -----------------------------------------------------------------------------------------
function change_CircCenter1(obj,eventdata,h)
% Change the Circle's center by asking it's coordinates
% ONLY FOR CARTESIAN CIRCLES.
	lon_lat_rad = getappdata(h,'LonLatRad');
	prompt = {'Enter new lon (or x)' ,'Enter new lat (or y)', 'Enter new radius'};
	num_lines= [1 30; 1 30; 1 30];
	def = {num2str(lon_lat_rad(1)) num2str(lon_lat_rad(2)) num2str(lon_lat_rad(3))};
	resp  = inputdlg(prompt, 'Change circle', num_lines, def);
	if isempty(resp),		return,		end
	np = numel(get(h,'XData'));		x = linspace(-pi,pi,np);
	y = sin(x);						x = cos(x);			% unit circle coords
	x = str2double(resp{1}) + str2double(resp{3}) * x;
	y = str2double(resp{2}) + str2double(resp{3}) * y;
	set(h, 'XData', x, 'YData', y);
	setappdata(h,'LonLatRad',[str2double(resp{1}) str2double(resp{2}) str2double(resp{3})])
% -----------------------------------------------------------------------------------------

% -----------------------------------------------------------------------------------------
function rectangle_limits(obj,eventdata)
% Change the Rectangle's limits by asking it's corner coordinates
	h = gco;
	x = get(h,'XData');     y = get(h,'YData');

	region = bg_region('with_limits',[x(1) x(3) y(1) y(3)]);
	if isempty(region),    return;  end     % User gave up
	x_min = region(1);      x_max = region(2);
	y_min = region(3);      y_max = region(4);

	set(h, 'XData', [x_min,x_min,x_max,x_max,x_min], 'YData', [y_min,y_max,y_max,y_min,y_min]);
% -----------------------------------------------------------------------------------------

% -----------------------------------------------------------------------------------------
function rectangle_register_img(obj, event)
% Prompt user for rectangle corner coordinates and use them to register the image
	do_flip = false;

	h = gco;
	handles = guidata(get(h,'Parent'));
	rect_x = get(h,'XData');   rect_y = get(h,'YData');		% Get rectangle limits

	region = bg_region('emptyREGIST');
	if isempty(region),		return,		end		% User gave up
	x_min = region(1);		x_max = region(2);
	y_min = region(3);		y_max = region(4);
	handles.geog = aux_funs('guessGeog',region(1:4));		% Trust more in the test here
	ax = handles.axes1;

	x(1) = rect_x(1);     x(2) = rect_x(2);     x(3) = rect_x(3);
	y(1) = rect_y(1);     y(2) = rect_y(2);
	img = get(handles.hImg,'CData');
	% Transform the ractangle limits into row-col limits
	limits = getappdata(handles.axes1,'ThisImageLims');
	limits(1:2) = limits(1:2) + [0.5 -0.5]*handles.head(8);	% Use the grid registration limits.
	limits(3:4) = limits(3:4) + [0.5 -0.5]*handles.head(9);
	r_c = cropimg(limits(1:2), limits(3:4), img, [x(1) y(1) (x(3)-x(2)) (y(2)-y(1))], 'out_precise');
	% Find if we are dealing with a image with origin at upper left (i.e. with y positive down)
	if(strcmp(get(ax,'YDir'),'reverse'))
		img = flipdim(img,1);
		% We have to invert the row count to account for the new origin in lower left corner
		imgHeight = size(img,1);
		tmp = r_c(1);
		r_c(1) = imgHeight - r_c(2) + 1;
		r_c(2) = imgHeight - tmp + 1;
		do_flip = true;
	end
	% Compute and apply the affine transformation
	base  = [x_min y_min; x_min y_max; x_max y_max; x_max y_min];
	input = [r_c(3) r_c(1); r_c(3) r_c(2); r_c(4) r_c(2); r_c(4) r_c(1)];

	trans = AffineTransform(input,base);
	x_pt = [1; size(img,2)];    y_pt = [1; size(img,1)];    % For more X points, change accordingly
	X1 = [x_pt y_pt ones(size(x_pt,1),1)];
	U1 = X1 * trans;
	new_xlim = U1(:,1)';        new_ylim = U1(:,2)';

	% Get coords of all line objects and recreate them later. Have to do it before deleting hImg
	ALLlineHand = findobj(handles.axes1,'Type','line');
	line_x = cell(1, numel(ALLlineHand));			line_y = cell(1, numel(ALLlineHand));
	for (k = 1:numel(ALLlineHand))
		x = get(ALLlineHand(k), 'XData');			y = get(ALLlineHand(k), 'YData');
		if (do_flip),		y = imgHeight - y + 1;	end
		X1 = [x(:) y(:) ones(numel(x),1)];
		U1 = X1 * trans;					% Reproject them
		line_x{k} = U1(:,1);	line_y{k} = U1(:,2);
	end

	% Rebuild the image with the new limits. After many atempts I found that kill and redraw is the safer way.
	m = size(img,1);		n = size(img,2);
	[new_xlim,new_ylim] = aux_funs('adjust_lims',new_xlim,new_ylim,m,n);
	delete(handles.hImg);
	handles.hImg = image(new_xlim,new_ylim,img,'Parent',handles.axes1);
	handles.head(1:4) = [new_xlim new_ylim];
	resizetrue(handles, [], 'xy');
	if (new_ylim(2) < new_ylim(1)),		set(ax,'YDir','reverse'),	end		% resizetrue doesn't know how to do this
	if (new_xlim(2) < new_xlim(1)),		set(ax,'XDir','reverse'),	end
	setappdata(ax,'ThisImageLims',[get(ax,'XLim') get(ax,'YLim')])
	handles.old_size = get(handles.figure1,'Pos');      % Save fig size to prevent maximizing
	handles.origFig = img;

	% Redraw the rectangle that meanwhile has gone to the ether togheter with gca.
	lt = handles.DefLineThick;				lc = handles.DefLineColor;
	for (k = 1:numel(ALLlineHand))
		ALLlineHand(k) = line('XData',line_x{k},'YData',line_y{k},'Color',lc,'LineWidth',lt,'Parent',handles.axes1);
	end
	if (handles.image_type == 2)					% Lets pretend that we have a GeoTIFF image
		handles.image_type = 3;
		Hdr.LL_prj_xmin = new_xlim(1);		Hdr.LR_prj_xmax = new_xlim(2);
		Hdr.LL_prj_ymin = new_ylim(1);		Hdr.UR_prj_ymax = new_ylim(2);
		Hdr.projection = 'linear';			Hdr.datum = 'unknown';
		set(handles.figure1,'UserData',Hdr);		% Minimalist Hdr to allow saving as a GeoTIFF image
	end
	x_inc = (new_xlim(2)-new_xlim(1)) / (size(img,2) - 1);
	y_inc = (new_ylim(2)-new_ylim(1)) / (size(img,1) - 1);
	handles.head(8:9) = [x_inc y_inc];   

	handles.fileName = [];							% Not loadable in session
	if (handles.validGrid)
		new_xlim = linspace(new_xlim(1),new_xlim(2),size(img,2));
		new_ylim = linspace(new_ylim(1),new_ylim(2),size(img,1));
		setappdata(handles.figure1,'dem_x',new_xlim);  	setappdata(handles.figure1,'dem_y',new_ylim);
	end

	if (handles.geog)
		set(handles.toGE,'Enable','on')
		mirone('SetAxesNumericType',handles,[]);	% Set axes uicontextmenus
	end
	guidata(handles.figure1, handles);
	for (k = 1:numel(ALLlineHand))
		draw_funs(ALLlineHand(k),'line_uicontext')	% Set lines's uicontextmenu
	end

% -----------------------------------------------------------------------------------------
function trans = AffineTransform(uv,xy)
% For an affine transformation:
%                     [ A D 0 ]
% [u v 1] = [x y 1] * [ B E 0 ]
%                     [ C F 1 ]
% There are 6 unknowns: A,B,C,D,E,F
% Another way to write this is:
%                   [ A D ]
% [u v] = [x y 1] * [ B E ]
%                   [ C F ]
% Rewriting the above matrix equation:
% U = X * T, where T = reshape([A B C D E F],3,2)
%
% With 3 or more correspondence points we can solve for T,
% T = X\U which gives us the first 2 columns of T, and
% we know the third column must be [0 0 1]'.

	K = 3;      M = size(xy,1);     X = [xy ones(M,1)];
	U = uv;		% just solve for the first two columns of T

	% We know that X * T = U
	if rank(X) >= K
		Tinv = X \ U;
	else
		errordlg('At least 3 non-collinear points needed to infer affine transform.','Error');
	end

	Tinv(:,3) = [0 0 1]';       % add third column
	trans = inv(Tinv);
	trans(:,3) = [0 0 1]';

% -----------------------------------------------------------------------------------------
function copy_text_object(obj, eventdata)
	hAx = gca;
	handMir = guidata(hAx);
	if (handMir.version7 < 8.4)
		copyobj(gco, hAx);
	else
		copyobj(gco, hAx, 'legacy');	% R2015 Fcker
	end
    move_text([],[])

% -----------------------------------------------------------------------------------------
function move_text(obj,evt)
	h = gco;
    hFig = get(0,'CurrentFigure');  hAxes = get(hFig,'CurrentAxes');
	state = uisuspend_j(hFig);     % Remember initial figure state
	set(hFig,'WindowButtonMotionFcn',{@wbm_txt,h,hAxes},'WindowButtonDownFcn',{@wbd_txt,h,state,hAxes});
	refresh
function wbm_txt(obj,eventdata,h,hAxes)
	pt = get(hAxes, 'CurrentPoint');
	pos = get(h,'Position');    pos(1) = pt(1,1);   pos(2) = pt(1,2);   set(h,'Position',pos);
	refresh
function wbd_txt(obj,eventdata,h,state,hAxes)
% check if x,y is inside of axis
	pt = get(hAxes, 'CurrentPoint');  x = pt(1,1);    y = pt(1,2);
	x_lim = get(hAxes,'xlim');      y_lim = get(hAxes,'ylim');
	if (x < x_lim(1)) || (x > x_lim(2)) || (y < y_lim(1)) || (y > y_lim(2));	return,	end
	refresh
	uirestore_j(state, 'nochildren');	% Restore the figure's initial state
% -----------------------------------------------------------------------------------------

% -----------------------------------------------------------------------------------------
function rotate_text(obj,eventdata)
	resp  = inputdlg({'Enter angle of rotation'}, '', [1 30]);
	if isempty(resp),	return,		end
	set(gco,'Rotation',str2double(resp))
	refresh

% -----------------------------------------------------------------------------------------
function text_FontSize(obj,eventdata)
	h = gco;
	ft = uisetfont(h,'Change Font');
	if (~isstruct(ft) && ft == 0),	return,		end
	set(h,'FontName',ft.FontName,'FontUnits',ft.FontUnits,'FontSize',ft.FontSize, ...
		'FontWeight',ft.FontWeight,'FontAngle',ft.FontAngle)
	refresh

% -----------------------------------------------------------------------------------------
function export_text(obj,eventdata)
	h = gco;		handles = guidata(h);
	pos = get(h,'Position');    font = get(h,'FontName');      size = get(h,'FontSize');
	str = get(h,'String');      angle = get(h,'Rotation');

	str1 = {'*.txt;*.TXT', 'Text file (*.txt,*.TXT)'; '*.*', 'All Files (*.*)'};
	[FileName,PathName] = put_or_get_file(handles,str1,'Select Text File name','put','.txt');
	if (isequal(FileName,0)),	refresh,	return,		end
	fname = [PathName FileName];
	fid = fopen(fname, 'w');
	if (fid < 0),   errordlg(['Can''t open file:  ' fname],'Error');    return;     end
	fprintf(fid,'%g\t%g\t%g\t%g\t%s ML  %s\n',pos(1), pos(2),size,angle,font,str);
	fclose(fid);

% -----------------------------------------------------------------------------------------
function set_symbol_uicontext(h,data)
% Set uicontexts for the symbols. h is the handle to the marker (line in fact) object
% This funtion is a bit messy because it serves for setting uicontexes of individual
% symbols, points and of "volcano", "hotspot" & "ODP" class symbols. 
	if (isempty(h)),	return,		end
	tag = get(h,'Tag');
	if (isa(tag, 'cell')),	tag = tag{1};	end
	if (numel(h) == 1 && length(get(h,'XData')) > 1)
		more_than_one = true;		% Flags that h points to a multi-vertice object
	else
		more_than_one = false;
	end

	handles = guidata(h(1));
	cmenuHand = uicontextmenu('Parent',handles.figure1);	set(h, 'UIContextMenu', cmenuHand);
	separator = false;
	this_not = true;		% for class symbols "this_not = 1". Used for not seting some options inapropriate to class symbols
	seismicity_options = false;
	tide_options = false;
	sep = 'off';

	if strcmp(tag,'hotspot')		% DATA must be a structure containing name & age fields
		uimenu(cmenuHand, 'Label', 'Hotspot info', 'Callback', {@hotspot_info,h,data.name,data.age,[]});
		uimenu(cmenuHand, 'Label', 'Plot name', 'Callback', {@hotspot_info,h,data.name,data.age,'text'});
		separator = true;
	elseif strcmp(tag,'volcano')    % DATA must be a structure containing name, description & dating fields
		uimenu(cmenuHand, 'Label', 'Volcano info', 'Callback', {@volcano_info,h,data.name,data.desc,data.dating});
		separator = true;
	elseif strcmp(tag,'meteor')		% DATA must be a structure containing name, diameter & dating fields
		uimenu(cmenuHand, 'Label', 'Impact info', 'Callback', {@meteor_info,h,data.name,data.diameter,data.dating,data.exposed,data.btype});
		separator = true;
	elseif strcmp(tag,'mar_online')	% DATA must be a structure containing name, code station & country fields
		uimenu(cmenuHand, 'Label', 'Download Mareg (2 days)', 'Callback', {@mareg_online,h,data});
		uimenu(cmenuHand, 'Label', 'Download Mareg (Calendar)', 'Callback', {@mareg_online,h,data,'cal'});
		separator = true;
	elseif strcmp(tag,'hydro')		% DATA must be a cell array with 5 cols contining description of each Vent
		uimenu(cmenuHand, 'Label', 'Hydrotermal info', 'Callback', {@hydro_info,h,data});
		separator = true;	
	elseif (strcmp(tag,'DSDP') || strcmp(tag,'ODP') || strcmp(tag,'IODP'))	% DATA is a struct with leg, site, z & penetration fields
		uimenu(cmenuHand, 'Label', [tag ' info'], 'Callback', {@ODP_info,h,data.leg,data.site,data.z,data.penetration});
		separator = true;
	elseif strcmp(tag,'City_major') || strcmp(tag,'City_other')
		this_not = true;
	elseif strcmp(tag,'Earthquakes')	% DATA is empty because I didn't store any info (they are too many)
		seismicity_options = isappdata(h,'SeismicityTime');
	elseif strcmp(tag,'Pointpolyline')	% DATA is empty because it doesn't have any associated info
		s_info = getappdata(h, 'LineInfo');		% But there might be one in appdata
		if (~isempty(s_info))
			uimenu(cmenuHand, 'Label', 'Symbol info', 'Callback', {@Isochrons_Info,s_info});
			sep = 'on';
		end
		this_not = false;
	elseif strcmp(tag,'TTT')			% DATA is empty
		this_not = false;
	elseif strcmp(tag,'TideStation')	% DATA is empty
		tide_options = true;
		separator = false;
	elseif strcmp(tag,'Anchor')			% The Anchor spot of the "Tiles Tool"
		uimenu(cmenuHand, 'Label', 'Move (precise)', 'Callback', {@change_SymbPos,h});
		return
	else
		this_not = false;
	end

	if (~this_not)						% non class symbols can be moved
		ui_edit_polygon(h)				% Set edition functions
		if strcmp(tag,'ARGO')
			uimenu(cmenuHand, 'Label', 'Plot data', 'Callback', {@ARGO_profile,h});
		else
%			uimenu(cmenuHand, 'Label', 'LOS', 'Callback', 'line_of_sight(gco)');
			uimenu(cmenuHand, 'Label', 'Move (precise)', 'Callback', {@change_SymbPos,h}, 'Sep', sep);
		end
	end

	if separator
		if (~more_than_one)				% Single symbol
			uimenu(cmenuHand, 'Label', 'Remove', 'Callback', 'delete(gco)', 'Sep','on');
		else							% Multiple symbols
			uimenu(cmenuHand, 'Label', 'Remove this', 'Callback', {@remove_one_from_many,h}, 'Sep','on');
		end
	else
		if (~more_than_one)				% Single symbol
			uimenu(cmenuHand, 'Label', 'Remove', 'Callback', 'delete(gco)');
		else							% Multiple symbols
			uimenu(cmenuHand, 'Label', 'Remove this', 'Callback', {@remove_one_from_many,h});
			uimenu(cmenuHand, 'Label', 'Delete all', 'Callback', {@del_line,h});
		end
	end

	if (this_not)						% individual symbols don't belong to a class
		uimenu(cmenuHand, 'Label', 'Remove class', 'Callback', {@remove_symbolClass,h});
	end

	if strcmp(tag,'LinkedSymb')			% Symbols plotted by the "Linked" option in either 'ecran' or 'gmtedit'
		uimenu(cmenuHand, 'Label', 'Remove class', 'Callback', {@remove_symbolClass,h});
		uimenu(cmenuHand, 'Label', 'Save all linked pts', 'Callback', {@save_line,h});
		this_not = true;				% Set this to avoid going inside next IF-test
	end

	if (~this_not)						% class symbols don't export
		uimenu(cmenuHand, 'Label', 'Save', 'Callback', {@export_symbol,h});
		if (handles.image_type ~= 20)
			if (strcmp(tag,'Pointpolyline'))	% Allow pure grdtrack interpolation
				cbTrack = 'setappdata(gcf,''TrackThisLine'',gco); mirone(''ExtractProfile_CB'',guidata(gcbo),''point'')';
				uimenu(cmenuHand, 'Label', 'Point interpolation', 'Callback', cbTrack, 'Sep','on');
			end
			if (handles.nLayers > 1)
				cbTrack = 'setappdata(gcf,''TrackThisLine'',gco); mirone(''ExtractProfile_CB'',guidata(gcbo),''3D'')';
				uimenu(cmenuHand, 'Label', '3D interpolation', 'Callback', cbTrack, 'Sep','on');
			end
		end
	end

	% -------------- See if we may offer the Quick interpolation option ------------------
	z = get(h, 'ZData');
	if (isempty(z)),	z = getappdata(h, 'ZData');		end
	if (~isempty(z))
		item = uimenu(cmenuHand, 'Label', 'Quick grid', 'Sep', 'on');
		uimenu(item, 'Label', 'auto', 'Callback', {@pt_quick_grd,h,'a'});
		uimenu(item, 'Label', 'set increments', 'Callback', {@pt_quick_grd,h,'q'});
	end

	if (seismicity_options)
		uimenu(cmenuHand, 'Label', 'Save events', 'Callback', 'save_seismicity(gcf,gco)', 'Sep','on');
		uimenu(cmenuHand, 'Label', 'Seismicity movie', 'Callback', 'animate_seismicity(gcf,gco)');
		uimenu(cmenuHand, 'Label', 'Draw polygon', 'Callback', ...
		       'mirone(''DrawClosedPolygon_CB'',guidata(gcbo),''SeismicPolyg'')');
		uimenu(cmenuHand, 'Label', 'Draw seismic line', 'Callback', ...
		       'mirone(''DrawLine_CB'',guidata(gcbo),''SeismicLine'')');
		itemHist = uimenu(cmenuHand, 'Label','Histograms');
		uimenu(itemHist, 'Label', 'Guttenberg & Richter', 'Callback', 'histos_seis(gco,''GR'')');
		uimenu(itemHist, 'Label', 'Cumulative number', 'Callback', 'histos_seis(gco,''CH'')');
		uimenu(itemHist, 'Label', 'Cumulative moment', 'Callback', 'histos_seis(gco,''CM'')');
		uimenu(itemHist, 'Label', 'Magnitude', 'Callback', 'histos_seis(gco,''MH'')');
		uimenu(itemHist, 'Label', 'Time', 'Callback', 'histos_seis(gco,''TH'')');
		uimenu(itemHist, 'Label', 'Display in Table', 'Callback', 'histos_seis(gcf,''HT'')','Sep','on');
		%uimenu(itemHist, 'Label', 'Hour of day', 'Callback', 'histos_seis(gco,''HM'')');
		itemTime = uimenu(cmenuHand, 'Label','Time series');
		uimenu(itemTime, 'Label', 'Time magnitude', 'Callback', 'histos_seis(gco,''TM'')');
		uimenu(itemTime, 'Label', 'Time depth', 'Callback', 'histos_seis(gco,''TD'')');
		uimenu(cmenuHand, 'Label', 'Mc and b estimate', 'Callback', 'histos_seis(gcf,''BV'')');
		uimenu(cmenuHand, 'Label', 'Fit Omori law', 'Callback', 'histos_seis(gcf,''OL'')');
	end

	if (tide_options)
		uimenu(cmenuHand, 'Label', 'Plot tides', 'Callback', {@tidesStuff,h,'plot'}, 'Sep','on');
		uimenu(cmenuHand, 'Label', 'Station Info', 'Callback', {@tidesStuff,h,'info'});
		%uimenu(cmenuHand, 'Label', 'Tide Calendar', 'Callback', {@tidesStuff,h,'calendar'});
	end

	itemSymb = uimenu(cmenuHand, 'Label', 'Symbol', 'Sep','on');
	cb_mark = uictx_setMarker(h,'Marker');              % there are 13 uictx_setMarker outputs
	uimenu(itemSymb, 'Label', 'plus sign', 'Callback', cb_mark{1});
	uimenu(itemSymb, 'Label', 'circle', 'Callback', cb_mark{2});
	uimenu(itemSymb, 'Label', 'asterisk', 'Callback', cb_mark{3});
	uimenu(itemSymb, 'Label', 'point', 'Callback', cb_mark{4});
	uimenu(itemSymb, 'Label', 'cross', 'Callback', cb_mark{5});
	uimenu(itemSymb, 'Label', 'square', 'Callback', cb_mark{6});
	uimenu(itemSymb, 'Label', 'diamond', 'Callback', cb_mark{7});
	uimenu(itemSymb, 'Label', 'upward triangle', 'Callback', cb_mark{8});
	uimenu(itemSymb, 'Label', 'downward triangle', 'Callback', cb_mark{9});
	uimenu(itemSymb, 'Label', 'right triangle', 'Callback', cb_mark{10});
	uimenu(itemSymb, 'Label', 'left triangle', 'Callback', cb_mark{11});
	uimenu(itemSymb, 'Label', 'five-point star', 'Callback', cb_mark{12});
	uimenu(itemSymb, 'Label', 'six-point star', 'Callback', cb_mark{13});
	itemSize = uimenu(cmenuHand, 'Label', 'Size');
	cb_markSize = uictx_setMarker(h,'MarkerSize');              % there are 8 uictx_setMarker outputs
	uimenu(itemSize, 'Label', '7       pt', 'Callback', cb_markSize{1});
	uimenu(itemSize, 'Label', '8       pt', 'Callback', cb_markSize{2});
	uimenu(itemSize, 'Label', '9       pt', 'Callback', cb_markSize{3});
	uimenu(itemSize, 'Label', '10     pt', 'Callback', cb_markSize{4});
	uimenu(itemSize, 'Label', '12     pt', 'Callback', cb_markSize{5});
	uimenu(itemSize, 'Label', '14     pt', 'Callback', cb_markSize{6});
	uimenu(itemSize, 'Label', '16     pt', 'Callback', cb_markSize{7});
	uimenu(itemSize, 'Label', 'other...', 'Callback', cb_markSize{8});
	setLineColor(uimenu(cmenuHand, 'Label', 'Fill Color'), uictx_Class_LineColor(h,'MarkerFaceColor'))
	setLineColor(uimenu(cmenuHand, 'Label', 'Edge Color'), uictx_Class_LineColor(h,'MarkerEdgeColor'))

% -----------------------------------------------------------------------------------------
function change_SymbPos(obj,eventdata,h)
% Change the Symbol position by asking it's coordinates

tag = get(h,'Tag');
if (strcmp(tag,'Pointpolyline') || strcmp(tag,'Maregraph') || strcmp(tag,'TTT'))
	pt = get(gca,'CurrentPoint');
	xp = get(h,'XData');    yp = get(h,'YData');
	% Find out which symb was selected
	dif_x = xp - pt(1,1);   dif_y = yp - pt(1,2);
	dist = sqrt(dif_x.^2 + dif_y.^2);   clear dif_x dif_y;
	[B,IX] = sort(dist);    i = IX(1);  clear dist IX;
	xx = xp(i);             yy = yp(i);
	is_single = 0;
else                % Individual symbol
	xx = get(h,'XData');        yy = get(h,'YData');
	i = 1;
	is_single = 1;
end

% Show the coordinates with same format as the axes label
labelType = getappdata(gca,'LabelFormatType');
if (~isempty(labelType))
	switch labelType
		case 'DegMin'
			x_str = degree2dms(str2double( ddewhite(sprintf('%8f',xx)) ),'DDMM',0,'str');   % x_str is a structure with string fields
			y_str = degree2dms(str2double( ddewhite(sprintf('%8f',yy)) ),'DDMM',0,'str');
			xx = [x_str.dd ':' x_str.mm];
			yy = [y_str.dd ':' y_str.mm];
		case 'DegMinDec'
			x_str = degree2dms(str2double( ddewhite(sprintf('%8f',xx)) ),'DDMM.x',2,'str');
			y_str = degree2dms(str2double( ddewhite(sprintf('%8f',yy)) ),'DDMM.x',2,'str');
			xx = [x_str.dd ':' x_str.mm];
			yy = [y_str.dd ':' y_str.mm];
		case 'DegMinSec'
			x_str = degree2dms(str2double( ddewhite(sprintf('%8f',xx)) ),'DDMMSS',0,'str');
			y_str = degree2dms(str2double( ddewhite(sprintf('%8f',yy)) ),'DDMMSS',0,'str');
			xx = [x_str.dd ':' x_str.mm ':' x_str.ss];
			yy = [y_str.dd ':' y_str.mm ':' y_str.ss];
		case 'DegMinSecDec'
			x_str = degree2dms(str2double( ddewhite(sprintf('%8f',xx)) ),'DDMMSS.x',1,'str');
			y_str = degree2dms(str2double( ddewhite(sprintf('%8f',yy)) ),'DDMMSS.x',1,'str');
			xx = [x_str.dd ':' x_str.mm ':' x_str.ss];
			yy = [y_str.dd ':' y_str.mm ':' y_str.ss];
		otherwise
			xx = sprintf('%.9g',xx);		yy = sprintf('%.9g',yy);
	end
else   
	xx = sprintf('%.9g',xx);		yy = sprintf('%.9g',yy);
end

prompt = {'Enter new lon (or x)' ,'Enter new lat (or y)'};
resp  = inputdlg(prompt,'Move symbol',[1 30; 1 30],{xx yy});
if isempty(resp),		return,		end

val_x = test_dms(resp{1});			% See if coords were given in dd:mm or dd:mm:ss format
val_y = test_dms(resp{2});
x = 0;     y = 0;
for (k = 1:length(val_x)),  x = x + sign(str2double(val_x{1}))*abs(str2double(val_x{k})) / (60^(k-1));    end
for (k = 1:length(val_y)),  y = y + sign(str2double(val_y{1}))*abs(str2double(val_y{k})) / (60^(k-1));    end

if (is_single)		% Individual symbol   
	xp = x;			yp = y;
else				% Picked symbol from a list
	xp(i) = x;		yp(i) = y;
end
set(h, 'XData', xp, 'YData', yp);

% -----------------------------------------------------------------------------------------
function ARGO_profile(obj,evt,hFig)
% Download and display ARGO data associated with the selected buoy
	pt = get(gca,'CurrentPoint');
	xp = get(hFig,'XData');    yp = get(hFig,'YData');
	% Find out which marker was picked
	dif_x = xp - pt(1,1);   dif_y = yp - pt(1,2);
	dist = sqrt(dif_x.^2 + dif_y.^2);   clear dif_x dif_y;
	[B,IX] = sort(dist);    i = IX(1);  clear dist IX;
	ncfiles = getappdata(hFig, 'ncfiles');
	thisFile = ncfiles{i};
	[lat,lon,date,P,T,S,pn] = argo_floats('argodata',thisFile);
	if (all(isnan(lat))),	return,		end
	if (size(P{1},1) == 2)						% Data may come in a misery shape
		p = P{1}(2,:);		p(isnan(p)) = [];		P = [p P{1}(1,:)];
		t = T{1}(2,:);		t(isnan(t)) = [];		T = [t T{1}(1,:)];
		s = S{1}(2,:);		s(isnan(s)) = [];		S = [s S{1}(1,:)];
		clear p t s
	else
		P = P{1}(1,:);		T = T{1}(1,:);			S = S{1}(1,:);
	end
	if (all(isnan(T)))
		warndlg('This instrument has no data.','Warning'),	return
	end
	s.figSize = [350 500];	s.xlabel = 'Salinity';	s.ylabel = 'Temperature';	s.title = strtok(datestr(date));
	%s.fhandle = {'scatter', {S,T,30,P,'filled'}};
	hFig = ecran({s});
	hS = scatter(S,T,30,P,'filled');
	handEcran = guidata(hFig);
	handEcran.hLine = hS;
	setappdata(hFig, 'location', [lon lat])		% Save these for eventual later usage in Ecran
	setappdata(hFig, 'Pressure', P),	setappdata(hFig, 'Temp', T),	setappdata(hFig, 'Sal', S)
	h = findobj(hFig, 'Label','Analysis');
	hC = get(h, 'Children');
	set(hC(end), 'Label','Sound Velocity Profile', 'Tag','AnalysisSVP');	% Recycle this one
	delete(hC(2:end-1))										% First is add_uictx(). a not visible one
	delete(findobj(hFig, 'Label','Misc'))
	delete(findobj(hFig, 'Tag','isocs_but')),	delete(findobj(hFig, 'Tag','rectang_but')),	delete(findobj(hFig, 'Tag','DynSlope'))
	guidata(hFig, handEcran)
	set(hFig, 'Name', ['TS diagram (' s.title ')'])
	%tightfig(hFig);

	% The PT fig
	s.xlabel = 'Temperature';	s.ylabel = 'Pressure (decibar)';
	hFig = ecran({s});
	hS = scatter(T,P,30,'filled');
	delete(findobj(hFig, 'Label','Misc')),		delete(findobj(hFig, 'Label','Analysis'))
	delete(findobj(hFig, 'Tag','isocs_but')),	delete(findobj(hFig, 'Tag','rectang_but')),	delete(findobj(hFig, 'Tag','DynSlope'))
	handEcran = guidata(hFig);		handEcran.hLine = hS;	guidata(hFig, handEcran)
 	hAx = findobj(hFig, 'Type', 'axes', 'Tag', 'axes1');
	set(hFig, 'Name', ['PT diagram (' s.title ')']),	set(hAx, 'YDir', 'reverse')
	%tightfig(hFig);
	pos = get(hFig, 'Pos');		pos(2) = pos(2) - (2-1)*30;		set(hFig, 'Pos', pos);

	% The PS fig
	hFig = ecran({s});
	hS = scatter(S,P,30,'filled');
	delete(findobj(hFig, 'Label','Misc')),		delete(findobj(hFig, 'Label','Analysis'))
	delete(findobj(hFig, 'Tag','isocs_but')),	delete(findobj(hFig, 'Tag','rectang_but')),	delete(findobj(hFig, 'Tag','DynSlope'))
	handEcran = guidata(hFig);		handEcran.hLine = hS;	guidata(hFig, handEcran)
	hAx = findobj(hFig, 'Type', 'axes', 'Tag', 'axes1');
 	set(get(hAx, 'XLabel'), 'Str','Salinity'),		set(get(hAx, 'YLabel'), 'Str','Pressure (decibar)')
	set(hFig, 'Name', ['PS diagram (' s.title ')']),	set(hAx, 'YDir', 'reverse')
	%tightfig(hFig)
	pos = get(hFig, 'Pos');		pos(2) = pos(2) - (3-1)*30;		set(hFig, 'Pos', pos);
	
% -----------------------------------------------------------------------------------------
function remove_one_from_many(obj,evt,h)
%Delete one symbol that belongs to a class (in fact a vertex of a polyline)
	pt = get(gca,'CurrentPoint');
	xp = get(h,'XData');    yp = get(h,'YData');
	l = numel(xp);
	if (iscell(xp))     % These stupids might be cell arrays, so they need to be converted to vectors   
		x = zeros(1,l);    y = zeros(1,l);   
		for i=1:l
			x(i) = xp{i};   y(i) = yp{i};
		end   
		xp = x;     yp = y;
	end
	% Find out which symb was selected
	dif_x = xp - pt(1,1);   dif_y = yp - pt(1,2);
	dist = sqrt(dif_x.^2 + dif_y.^2);   clear dif_x dif_y;
	[B,IX] = sort(dist);    i = IX(1);  clear dist IX;
	xp(i) = [];     yp(i) = [];
	set(h, 'XData',xp, 'YData',yp,'LineStyle','none');
	zz = get(h, 'UserData');
	if (~isempty(zz))
		zz(i) = [];
		set(h, 'UserData', zz)
	end

% -----------------------------------------------------------------------------------------
function other_SymbSize(obj,eventdata,h)
	resp  = inputdlg({'Enter new size (pt)'}, 'Symbol Size', [1 30]);
	if isempty(resp),		return,		end
	set(h,'MarkerSize',str2double(resp));        refresh

% -----------------------------------------------------------------------------------------
function res = check_IsRectangle(h)
% Check if h is a handle to a rectangle. This is used to verify if a rectangle has
% not been deformed (edited). We need it because some options are only available
% to operate with rectangles (e.g. Crop, Transplant Image, etc...)
	x = get(h,'XData');   y = get(h,'YData');
	if ~( (x(1) == x(end)) && (y(1) == y(end)) && length(x) == 5 && ...
			(x(1) == x(2)) && (x(3) == x(4)) && (y(1) == y(4)) && (y(2) == y(3)) )
		res = 0;
	else
		res = 1;
	end

% -----------------------------------------------------------------------------------------
function remove_symbolClass(obj,eventdata,h)
% Delete all symbol that belong to the class of "h". We do this by fishing it's tag.
% If individual symbols were previously removed (by "Remove symbol") h has invalid
% handles, so make sure all handles are valid
	h = h(ishandle(h));
	tag = get(h,'Tag');
	if iscell(tag)          % When several symbol of class "tag" exists
		h_all = findobj(gca,'Tag',tag{1});
	else                    % Only one symbol of class "tag" exists
		h_all = findobj(gca,'Tag',tag);
	end
	delete(h_all)

% -----------------------------------------------------------------------------------------
function make_age_script(obj, evt)
% ...
	pato = 'C:\a1\mgd77\AGU12\t\';
	hAllIsocs = findobj('Tag', get(gco,'Tag'));
	fidJob = fopen([pato 'generate_age_pts.bat'],'wt');			% The batch file
	fprintf(fidJob, '@echo off\n\nset fname=ages_stages.dat\n');
	fprintf(fidJob, 'set opt_D=-D0.03\n');
	fprintf(fidJob, 'set opt_S=-S0.03\n');
	fprintf(fidJob, 'set opt_P=-P90\n\n');
	fprintf(fidJob, 'del /Q %%fname%%\n');
	for (i = 1:numel(hAllIsocs))
		LineInfo = getappdata(hAllIsocs(i),'LineInfo');
		if (isempty(LineInfo)),		continue,	end
		%ind = strfind(LineInfo, 'STG1"');
		ind = strfind(LineInfo, 'STG3_');
		if (isempty(ind)),		continue,	end
		ind2 = strfind(LineInfo(ind:end),'"') + ind - 1;	% So that ind reports to the begining string too
		isoc_name = strtok(LineInfo);
		pole_name = sprintf('%spolo_%d.stg', pato, i);
		A = sscanf(LineInfo(ind2(1)+1:ind2(2)-1), '%f %f %f %f %f');
		stg = sprintf('%.9g %.9g %.9g %.9g %.9g', A(1),A(2),A(3)-A(4),0,A(5));
		fid = fopen(pole_name,'wt');					%  The stage pole file
		fprintf(fid, '%s\n', stg);
		fclose(fid);
		x = get(hAllIsocs(i), 'XData');		y = get(hAllIsocs(i), 'YData');
		isoc_name = sprintf('%sisoca_%d_%s.dat', pato, i, isoc_name);	% reuse the variable 'isoc_name'
		fid = fopen(isoc_name,'wt');					%  The isochrone file
		fprintf(fid, '> %s\n', LineInfo);
		fprintf(fid, '%.4f\t%.4f\n', [x(:)'; y(:)']);
		fclose(fid);
		fprintf(fidJob, 'telha %s -E%s -A %%opt_D%% %%opt_S%% -O%g %%opt_P%% >> %%fname%%\n', isoc_name, pole_name, A(4));

		% The Ridge is a special case because it's not duplicated and has two stage poles.
		if (strcmp(LineInfo(1:3), ' 0 '))				% Found the Ridge
			pole_name = sprintf('%spolo__%d.stg', pato, i);
			fid = fopen(pole_name,'wt');				%  The stage pole file
			stg = sprintf('%.9g %.9g %.9g %.9g %.9g', A(1),A(2),A(3)-A(4),0,-A(5));
			fprintf(fid, '%s\n', stg);
			fclose(fid);
			fprintf(fidJob, 'telha %s -E%s -A %%opt_D%% %%opt_S%% -O%g %%opt_P%% >> %%fname%%\n', isoc_name, pole_name, A(4));
		end
	end
	fclose(fidJob);

% -----------------------------------------------------------------------------------------
function remove_singleContour(obj, evt, h)
% Delete an individual contour and its eventual label(s)
	labHand = getappdata(h,'LabelHands');
	if (~isempty(labHand))
		try     delete(labHand);   end
	end
	delete(h)

% -----------------------------------------------------------------------------------------
function save_line(obj, evt, h)
% Save either individual as well as class lines. The latter uses the ">" symbol to separate segments

	if (nargin == 3),		h = h(ishandle(h));
	elseif (nargin == 2),	h = gco;			% Save only one line
	end
	handles = guidata(gcbo);
	name = '.dat';				% Default to this extension but if true isochron we try to construct a good name
	last_dir = '';
	if (strcmp(get(h, 'Tag'),'isochron'))	% Special case
		name = create_isoc_name(h);
		if (isempty(name)),		name = '.dat';	end		% Something did not work well
	else
		pato = getappdata(h, 'filePath');	% Now that we have this the other branch seems completely superfelous.
		if (~isempty(pato))
			last_dir = handles.last_dir;
			handles.last_dir = pato;		% Temporary change 
			fn = getappdata(h, 'fileName');
			if (~isempty(fn) && numel(h) == 1),		name = fn;	end
		end
	end
	str1 = {'*.dat;*.DAT', 'Line file (*.dat,*.DAT)'; '*.*', 'All Files (*.*)'};
	[FileName,PathName,handles] = put_or_get_file(handles,str1,'Select Line File name','put', name);
	if (~isempty(last_dir)),	handles.last_dir = last_dir;	end		% Reset the right dir
	guidata(handles.figure1,handles)
	if isequal(FileName,0),		return,		end
	x = get(h,'XData');			y = get(h,'YData');

	fname = [PathName FileName];
	[PATH,FNAME,EXT] = fileparts([PathName FileName]);
	if isempty(EXT),    fname = [PathName FNAME '.dat'];    end

	if (strcmp(get(h, 'Tag'),'GMT_DBpolyline'))
		save_GMT_DB_asc(h, fname)
		return
	end

	prj4 = '';	cvt2geogs = false;
	if (handles.is_projected)
		prj4 = aux_funs('get_proj_string',handles.figure1, 'proj');
		cvt2geogs = (getappdata(handles.figure1,'DispInGeogs') == 1);	% True when user selected to display in geogs
	end

	fid = fopen(fname, 'w');
	if (fid < 0),	errordlg(['Can''t open file:  ' fname],'Error'),	return,		end

	if (strcmp(get(h, 'Tag'),'LinkedSymb'))		% Those are rare cases (gmtedit & ecran), need to check it first
		hAll = findobj(handles.axes1,'-depth',1,'type','line','Tag','LinkedSymb');
		for (k = 1:numel(hAll))
			x = get(hAll(k),'XData');		y = get(hAll(k),'YData');	str = getappdata(hAll(k),'box');
			if (isempty(str)),	fprintf(fid,'%.6f\t%.6f\n',[x(:)'; y(:)']);
			else,				fprintf(fid,'%.6f\t%.6f\t%s\n',[x(:)'; y(:)'], str);
			end
		end
	elseif (~iscell(x))
		LineInfo = getappdata(h,'LineInfo');
		LineInfo_bak = '';
		if (isa(LineInfo, 'cell'))
			LineInfo_bak = LineInfo;	LineInfo = LineInfo{1};
		end
		if (~isempty(prj4))
			if (cvt2geogs)		% OK, so in this case we need to convert from the prj4
				xy = proj2proj_pts([], [x(:) y(:)], 'srcProj4', prj4, 'dstProj4', '+proj=longlat');
				x = xy(:,1);	y = xy(:,2);
				prj4 = '+proj=longlat +ellps=WGS84';
			end
			LineInfo = upd_LineInfo(LineInfo, ['"' prj4 '"']);
		end
		if (~isempty(LineInfo_bak))			% If LineInfo was a cellarray, restore it.
			LineInfo_bak{1} = LineInfo;		LineInfo = LineInfo_bak;
		end

		if (strcmp(LineInfo, '>')),	LineInfo = '';	end		% Sometimes we get a trailling one
		if (~isempty(LineInfo))
			if (isa(LineInfo, 'cell'))
				fprintf(fid,'> %s',LineInfo{1});
				for (k = 2:numel(LineInfo))
					fprintf(fid,'\t%s',LineInfo{k});
				end
				fprintf(fid,'\n');
			else
				fprintf(fid,'> %s\n',LineInfo);
			end
		end
		z = getappdata(h,'ZData');
		if (isempty(z) || numel(z) ~= numel(x))
			fprintf(fid,'%.6f\t%.6f\n',[x(:)'; y(:)']);
			if (~isempty(z) && numel(z) ~= numel(x))
				warndlg('Cannot save the Z column due to a bug in line editing.', 'Warning')
			end
		else
			fprintf(fid,'%.6f\t%.6f\t%.6f\n',[x(:)'; y(:)'; z(:)']);
		end
	else
		for (i = 1:numel(h))
			LineInfo = getappdata(h(i),'LineInfo');

			if (cvt2geogs)		% OK, so in this case we need to convert from the prj4
				if (i == 1)
					LineInfo = upd_LineInfo(LineInfo, '"+proj=longlat +ellps=WGS84"');
				end
				xy = proj2proj_pts([], [x{i}(:) y{i}(:)], 'srcProj4', prj4, 'dstProj4', '+proj=longlat');
				x{i} = xy(:,1);	y{i} = xy(:,2);
			elseif (i == 1 && ~isempty(prj4))
				LineInfo = sprintf('%s "%s"', LineInfo, prj4);		% Works even if LineInfo == ''
			end

			if (~isempty(LineInfo)),	str = ['> ' LineInfo];
			else,						str = '>';
			end
			fprintf(fid,'%s\n',str);

			z = getappdata(h(i),'ZData');
			if (isempty(z) || numel(z) ~= numel(x{i}))
				fprintf(fid,'%.6f\t%.6f\n',[x{i}(:)'; y{i}(:)']);
				if (~isempty(z) && numel(z) ~= numel(x{i}))
					warndlg('Cannot save the Z column due to a bug in line editing.', 'Warning')
				end
			else
				fprintf(fid,'%.6f\t%.6f\t%.6f\n',[x{i}(:)'; y{i}(:)'; z(:)']);
			end
		end
	end
	fclose(fid);

function name = create_isoc_name(h)
% Generate the likely name that I use to use for the isochrons
	try
		LineInfo = getappdata(h,'LineInfo');
		[isoc, plates] = strtok(LineInfo);	plates = ddewhite(plates);
		
		ind = strfind(plates, '/');
		if (isempty(ind)),		name = '';		return,		end		% No pair plates so no name fishing
		name1 = plates(1:ind(1)-1);
		[t,r] = strtok(plates(ind(1)+1:end));	r = ddewhite(r);
		if (~isempty(r))
			if (~strncmp(r, 'FIN', 3))		% Second name is made of two words
				name2 = [t ' ' r];
			else
				name2 = t;
			end
		else
			name2 = t;
		end
		
		ind = strfind(name1, ' ');
		if (isempty(ind)),		n1 = name1(1:2);
		else,					n1 = [name1(1) name1(ind(1)+1)];
		end
		ind = strfind(name2, ' ');
		if (isempty(ind)),		n2 = name2(1:2);
		else,					n2 = [name2(1) name2(ind(1)+1)];
		end
		name = sprintf('c%s_%s_%s.dat', isoc, n1, n2);
	catch
		name = [];
	end

function save_GMT_DB_asc(h, fname)
% Go through all GMT_DB polygons present in figure and save only those who have been edited
	fid = fopen(fname, 'w');
	if (fid < 0),	errordlg(['Can''t open file:  ' fname],'Error'),	return,		end
	h = findobj('Tag','GMT_DBpolyline');
	for (k = 1:numel(h))
		if (isempty(getappdata(h(k), 'edited'))),	continue,	end		% Skip because it was not modified
		GSHHS_str = getappdata(h(k),'GSHHS_str');
		if (k == 1 && ~isempty(GSHHS_str))		% Write back the magic string that allows us to recognize these type of files
			fprintf(fid,'# $Id: draw_funs.m 11441 2019-07-24 14:35:47Z j $\n#\n%s\n#\n', GSHHS_str);
		end
		hdr = getappdata(h(k), 'LineInfo');
		x = get(h(k), 'XData');			y = get(h(k), 'YData');
		indNaN = find(isnan(x));
		if (~isempty(indNaN)),	x(indNaN) = [];		y(indNaN) = [];		end	% We don't want NaNs in this story
		ind1 = strfind(hdr,'N = ');		ind2 = strfind(hdr,'G = ');
		hdr = sprintf('> %s%d %s', hdr(1:ind1+3), numel(x), hdr(ind2:end));
		fprintf(fid,'%s\n',hdr);
		fprintf(fid,'%.6f\t%.6f\n',[x; y]);
	end
	fclose(fid);
% -----------------------------------------------------------------------------------------

% -----------------------------------------------------------------------------------------
function export_symbol(obj, evt, h, opt)
	h = h(ishandle(h));
	tag = get(h,'Tag');
	xx = get(h,'XData');    yy = get(h,'YData');
	if (numel(xx) > 1 && ~strcmp(tag,'Pointpolyline') && ~strcmp(tag,'Maregraph'))     % (don't remember why)
		% Points and Maregraphs may be many but don't belong to a class
		msgbox('Only individual symbols may be exported and this one seams to belong to a class of symbols. Exiting','Warning')
		return
	end
	zz = get(h,'UserData');
	if (isempty(zz)),		doSave_formated(xx, yy)
	else,					doSave_formated(xx, yy, zz)
	end

% -----------------------------------------------------------------------------------------
function save_formated(obj, evt, h, opt)
% Save x,y[,z] vars into a file but taking into account the 'LabelFormatType'
% If OPT is given than it must contain either a:
%	- Mx3 array with the x,y,z data to be saved
%	- A struct with fields 'x', 'y', 'z', 'time_z' issued from the 3D interpolation
%	  of a multi-layered grid. See mirone.m ExtractProfile_CB() for details

	if (nargin < 3)
		errordlg('save_formated: called with a wrong number of arguments.','ERROR'),	return
	elseif (nargin == 3)
		if (~isa(h,'cell')),	h = gco;		% Fish the handle so that it works with copyied objs
		else,					h = h{1};		% Really use this handle.
		end
		xx = get(h,'XData');    yy = get(h,'YData');	zz = get(h, 'ZData');
		doSave_formated(xx, yy, zz, getappdata(h,'LineInfo'))
	else
		if (~isa(opt,'struct'))					% The Mx3 array case
			if (size(opt,2) ~= 3)
				errordlg('save_formated: variable must contain a Mx3 array.','ERROR')
				return
			end
			doSave_formated(opt(:,1), opt(:,2), opt(:,3))
		else									% The 3D interpolation case (this needs further doc)
			handles = guidata(get(0,'CurrentFigure'));
			[FileName,PathName] = put_or_get_file(handles,{'*.dat;*.DAT','ASCII file'; '*.*', 'All Files (*.*)'},'Output file','put');
			if isequal(FileName,0),		return,		end
			fname = [PathName FileName];
			wmode = 'wt';
			if (~ispc),		wmode = 'w';	end	% ASCII read mode type
			fid = fopen(fname, wmode);
			t = opt.time_z;						% Layers's times (or whatever unit)
			fprintf(fid,'# Interpolated file: %s\n', handles.grdname);
			fprintf(fid, ['#  \t', repmat('%g(X)\t', [1,size(opt.z,2)]) '\n'], opt.x(:));
			fprintf(fid, ['#  \t', repmat('%g(Y)\t', [1,size(opt.z,2)]) '\n'], opt.y(:));
			fprintf(fid, '>XY\n');
			fprintf(fid, ['%.2f\t' repmat('%f\t',[1,size(opt.z,2)]) '\n'], [t(:) double(opt.z)]');
			fclose(fid);
		end
	end

% -----------------------------------------------------------------------------------------
function doSave_formated(xx, yy, opt_z, LineInfo)
% Save x,y[,z] vars into a file but taking into account the 'LabelFormatType'
% OPT_Z is what the name says, optional

	if (nargin < 4),	LineInfo = '';	end

	hFig = get(0,'CurrentFigure');
	handles = guidata(hFig);
	str1 = {'*.dat;*.DAT', 'Symbol file (*.dat,*.DAT)'; '*.*', 'All Files (*.*)'};
	[FileName,PathName] = put_or_get_file(handles,str1,'Select File name','put','.dat');
	if isequal(FileName,0),		return,		end
	f_name = [PathName FileName];
	
	% Save data with a format determined by axes format
	labelType = getappdata(handles.axes1,'LabelFormatType');	% find the axes label format
	if isempty(labelType),		labelType = ' ';		end		% untempered matlab axes labels
	switch labelType
		case {' ','DegDec','NotGeog'}
			[xy, LineInfo] = lineinfo_proj(handles, xx, yy, LineInfo);
			if (isempty(xy)),	xy = [xx(:) yy(:)];		end		% Means it was not projected
			fmt = '%f\t%f';
		case 'DegMin'
			out_x = degree2dms(xx,'DDMM',0,'numeric');        out_y = degree2dms(yy,'DDMM',0,'numeric');
			xy = [out_x.dd(:) out_x.mm(:) out_y.dd(:) out_y.mm(:)];
			fmt = '%4d %02d\t%4d %02d';
		case 'DegMinDec'        % I'm writing the minutes with a precision of 2 decimals
			out_x = degree2dms(xx,'DDMM.x',2,'numeric');      out_y = degree2dms(yy,'DDMM.x',2,'numeric');
			xy = [out_x.dd(:) out_x.mm(:) out_y.dd(:) out_y.mm(:)];
			fmt = '%4d %02.2f\t%4d %02.2f';
		case 'DegMinSec'
			out_x = degree2dms(xx,'DDMMSS',0,'numeric');      out_y = degree2dms(yy,'DDMMSS',0,'numeric');
			xy = [out_x.dd(:) out_x.mm(:) out_x.ss(:) out_y.dd(:) out_y.mm(:) out_y.ss(:)];
			fmt = '%4d %02d %02d\t%4d %02d %02d';
		case 'DegMinSecDec'     % I'm writing the seconds with a precision of 2 decimals
			out_x = degree2dms(xx,'DDMMSS',2,'numeric');      out_y = degree2dms(yy,'DDMMSS',2,'numeric');
			xy = [out_x.dd(:) out_x.mm(:) out_x.ss(:) out_y.dd(:) out_y.mm(:) out_y.ss(:)];
			fmt = '%4d %02d %02.2f\t%4d %02d %02.2f';
		case 'Date'
			xy = [xx(:) yy(:)];
			fmt = '%.14g\t%f';
	end
	
	if (nargin == 3 && ~isempty(opt_z))      
		xy = [xy opt_z(:)];    fmt = [fmt '\t%f'];
	end
	if (~isempty(LineInfo))		% Add a line info multiseg first line
		fmt = {LineInfo, fmt};
	end
	double2ascii(f_name,xy,fmt,'maybeMultis');

% -----------------------------------------------------------------------------------------
function [xy, LineInfo] = lineinfo_proj(handles, x, y, LineInfo)
% Conditionally project x,y back to geogs and also update or create LineInfo with prj4 info
	xy = [];
	if (handles.is_projected)
		prj4 = aux_funs('get_proj_string',handles.figure1, 'proj');
		cvt2geogs = (getappdata(handles.figure1,'DispInGeogs') == 1);	% User selected to display in geogs
		if (~isempty(prj4) && cvt2geogs)		% Convert to geogs
			xy = proj2proj_pts([], [x(:) y(:)], 'srcProj4', prj4, 'dstProj4', '+proj=longlat');
			if (~isempty(LineInfo))
				LineInfo = upd_LineInfo(LineInfo, '"+proj=longlat +ellps=WGS84"');
			end
		elseif (~isempty(prj4))		% prj4 should not be, but...
			LineInfo = upd_LineInfo(LineInfo, prj4);
		end
	end
	
	function LineInfo = upd_LineInfo(LineInfo, new_prj4)
	% 1. Update the +proj4 string if it exists in LineInfo
	% 2. Append the new_prj4 string if LineInfo is not empty
	% 3. Create a LineInfo with the new_prj4 if LineInfo is empty
		ind = strfind(LineInfo, '"+proj');
		if (~isempty(ind))		% A +proj=... already exists. Replace it with the longlat
			ind2 = strfind(LineInfo(ind(1)+1:end), '"');	% This guy MUST exist
			if (isempty(ind2))
				error('+proj string in header MUST be wraped with double quotes')
			end
			LineInfo = strrep(LineInfo, LineInfo(ind(1):ind(1)+ind2(1)), [' ' new_prj4]);
		elseif (~isempty(LineInfo))
			LineInfo = [LineInfo [' ' new_prj4]];
		else
			if (new_prj4(1) == '"'),	LineInfo = ['> ' new_prj4];
			else,						LineInfo = ['> "' new_prj4 '"'];
			end
		end

% -----------------------------------------------------------------------------------------
function cb = uictx_SymbColor(h,prop)
% Set uicontext colors in object hose handle is gco (or h for "other color")
	cb{1} = ['set(gco, ''' prop ''', ''k'');refresh'];  cb{2} = ['set(gco, ''' prop ''', ''w'');refresh'];
	cb{3} = ['set(gco, ''' prop ''', ''r'');refresh'];  cb{4} = ['set(gco, ''' prop ''', ''g'');refresh'];
	cb{5} = ['set(gco, ''' prop ''', ''b'');refresh'];  cb{6} = ['set(gco, ''' prop ''', ''y'');refresh'];
	cb{7} = ['set(gco, ''' prop ''', ''c'');refresh'];  cb{8} = ['set(gco, ''' prop ''', ''m'');refresh'];
	cb{9} = {@other_SymbColor,h,prop};
% -----------------------------------------------------------------------------------------
function other_SymbColor(obj,eventdata,h,prop)
	c = uisetcolor;
	if length(c) > 1            % That is, if a color was selected
		set(h,prop,c)
	end

% -----------------------------------------------------------------------------------------
function hotspot_info(obj,eventdata,h,name,age,opt)
	i = get(gco,'Userdata');
	if isempty(opt)
		msgbox( sprintf(['Hotspot name: ' name{i} '\n' 'Hotspot age:   ' sprintf('%g',age(i)) ' Ma'] ),'Fogspot info')
	else
		name = strrep(name,'_',' ');            % Replace '_' by ' '
		textHand = text(get(h(i),'XData'),get(h(i),'YData'),0,name{i});
		draw_funs(textHand,'DrawText')          % Set text's uicontextmenu
	end

% -----------------------------------------------------------------------------------------
function tidesStuff(obj,eventdata,h,opt)
	pt = get(gca,'CurrentPoint');
	if (strcmp(opt,'plot'))
		t_xtide(pt(1,1),pt(1,2));
	elseif (strcmp(opt,'info'))
		info = t_xtide(pt(1,1),pt(1,2),'format','info');
		str{1} = info.station;
		str{2} = ['Position: Lon = ' num2str(info.longitude) '  Lat = ' num2str(info.latitude)];
		str{3} = ['Timezone: UTC ' num2str(info.timezone)];
		str{4} = ['Datum: ' num2str(info.datum)];
		str{5} = ['Number of constit = ' num2str(length(info.freq))];
		msgbox(str,'Station info')
	% elseif (strcmp(opt,'calendar'))
	%     date = clock;
	%     tim = datenum(date(1),date(2),1):1/24:datenum(date(1),date(2),31);
	%     out = t_xtide(pt(1,1),pt(1,2),tim,'format','times');
	end

% -----------------------------------------------------------------------------------------
function mareg_online(obj,eventdata,h, data, opt)
	i = get(gco,'Userdata');
	handles = guidata(h(1));
	dest_fiche = [handles.path_tmp 'lixo.dat'];
	nome = strrep(data.name{i},'_',' ');		pais = strrep(data.country{i},'_',' ');
	code = data.codeSt{i};

	c = clock;
	ds_num = datevecmx(datenummx(c(1:6))-2);
	date_start = sprintf('%04d-%02d-%02d %02d:%02d:%02d', [ds_num(1:5) fix(ds_num(6))]);	% BUG

	if (nargin == 4)	% Get last 2 days of records
		nDays = '2';
	else
		prompt = {'Start date (mind the format)' 'number of days'};
		date_str_tmp = sprintf('%02d-%02d-%04d %02d:%02d', [c(3:-1:1) c(4:5)]);
		resp   = inputdlg(prompt,'Select date',[1 30; 1 30],{date_str_tmp '2'});	pause(0.01);
		if isempty(resp),	return,		end
		nDays = resp{2};
		t = sscanf(resp{1}, '%02d-%02d-%04d %02d:%02d');	t = t(:)';
		t = datevecmx(datenummx([t(3:-1:1) t(4:5) 0]) + str2double(nDays));
		date_start = sprintf('%04d-%02d-%02d %02d:%02d:00', t(1:5));
	end
	url = ['http://www.ioc-sealevelmonitoring.org/bgraph.php?code=' code '&output=asc&period=' nDays '&endtime=' date_start];
	%url = ['http://www.ioc-sealevelmonitoring.org/bgraph.php?output=asc&time=' date_start '&period=' nDays '&par=' code];
	if (ispc),		dos(['wget "' url '" -q --tries=2 --connect-timeout=5 --no-check-certificate -O ' dest_fiche]);
	else,			unix(['wget ''' url ''' -q --tries=2 --connect-timeout=5 --no-check-certificate -O ' dest_fiche]);
	end

	fid = fopen(dest_fiche,'r');
	todos = fread(fid,'*char');
	if (numel(todos) < 100)
		warndlg(sprintf('This station has no data or a file tranfer error occured.\nYou may inquire this web site to see why\n%s',url),'Warning')
		fclose(fid);		return
	end
	ind = strfind(todos(1:128)',sprintf('\n'));
	if (numel(ind) == 1)
		warndlg('Sorry, but the downloaded file is not organized in the standard way. Quiting.','Warning')
		return
	end
	what = strread(todos(ind(1)+1:ind(2)-1),'%s', 'delimiter', '\t');
	indVar = find(strcmp(what, 'prs(m)'));
	if (isempty(indVar))
		indVar = find(strcmp(what, 'rad(m)'));
		if (isempty(indVar))
			fclose(fid);
			errordlg('Sorry, but no ''prs'' or ''rad'' variables in this file. Quiting.','Error')
			return
		end
	end
	todos = todos(ind(2)+1:end);		% Jump the 2 header lines
	yymmdd = strread(todos,'%s', 'delimiter', '\t');
	fclose(fid);    clear todos
	indDate = 1:numel(what):numel(yymmdd);		% Index of the Date lines
	sl = str2double(yymmdd(indDate+indVar-1));
	yymmdd = yymmdd(indDate);
	m = numel(yymmdd);		y = zeros(m,6);
	for (k = 1:m)
		t = sscanf(yymmdd{k}, '%04d-%02d-%02d %02d:%02d');		t = t(:)';
		y(k,1:5) = t;
	end
	serial_date = datenummx(y);
	hf = ecran(handles, serial_date, sl, [nome ' (' pais ')']);
	hAxes = findobj(hf,'Tag','axes1');
	setappdata(hAxes, 'LabelFormatType', 'Date');		% Tell pixval_stsbar to display XX coords as Date-Time
	h = findobj(hf,'Tag','add_uictx');
	cb = get(h, 'Callback');
	feval(cb, h, guidata(hf))			% Call the ecran's add_uictx_CB function

% -----------------------------------------------------------------------------------------
function meteor_info(obj,eventdata,h, name, diameter, dating, exposed, btype)
	i = get(gco,'Userdata');
	nome = strrep(name{i},'_',' ');			dating = strrep(dating{i},'_',' ');
	str = sprintf('Impact name:   %s\nDiameter (km):   %s\nAge (Ma):   %s\nExposed:    %s', nome, diameter{i}, dating, exposed{i});
	if (btype{i}(1) ~= '-')
		btype = strrep(btype{i},'_',' ');
		str = sprintf('%s\nBolid Type:   %s', str, btype);
	end
	msgbox( sprintf(str),'Impact info')

% -----------------------------------------------------------------------------------------
function hydro_info(obj,eventdata, h, desc)
	i = get(gco,'Userdata');
	str = sprintf('Vent name:   %s\nDepth (m):   %s\nActivity:   %s', ...
		desc{i,1}, desc{i,2}, desc{i,3});
	msgbox( sprintf(str),'Hydrothermal info')

% -----------------------------------------------------------------------------------------
function volcano_info(obj,eventdata, h, name, desc, dating)
	i = get(gco,'Userdata');
	msgbox( sprintf(['Volcano name: ' name{i} '\n' 'Volcano type:   ' desc{i} '\n' ...
            'Activity:      ' dating{i}] ),'Volcano info')

% -----------------------------------------------------------------------------------------
function ODP_info(obj,eventdata,h,leg,site,z,penetration)
	i = get(gco,'Userdata');
	tag = get(h,'Tag');     tag = tag{1};
	msgbox( sprintf([[tag ' Leg:    '] '%d\n' [tag ' Site:    '] site{i} '\n' ...
			'Depth:        %d\n' 'Hole penetration: %d'], ...
			double(leg(i)), double(z(i)), double(penetration(i)) ),'ODP info')

% -----------------------------------------------------------------------------------------
function str = Isochrons_Info(obj, evt, data)
% Display the info associated with the current object (a pline). get(OBJ, 'Type') == uimenu
	hLine = gco;
	LineInfo = getappdata(hLine,'LineInfo');	% Try first to get the info directly fom appdata
	if (isempty(LineInfo))						% If it fails, fall back to old mechanism
		i = get(hLine,'Userdata');
		if (isstruct(i))			% This happens when h is ui_edit_polygon(ed)
			i = i.old_ud;
		end
		if (~isempty(i)),	LineInfo = data{i};	end
	end
	if (nargout),	str = LineInfo;	return,		end

	if (~isempty(LineInfo))
		if (isa(LineInfo, 'cell')),		LineInfo = LineInfo{1};		end
		if (strncmp(LineInfo, 'MB ', 3))
			fname = LineInfo(4:end);
			WrapString = textwrap({fname}, 100);
			pato = fileparts(fname);
			if (isempty(pato))
				pato = getappdata(hLine, 'filePath');		fname = [pato filesep fname];
			end
			DefFigPos=get(0,'DefaultFigurePosition');
			h = dialog('Position',[DefFigPos(1:2) 300 130],'Name','Info', 'WindowStyle','normal');
			uicontrol('Parent',h, 'Style','text', 'Position',[10 70 210 50], 'String', WrapString);
			uicontrol('Parent',h, 'Position',[10 40 130 23], 'String','Show as point cloud', 'Callback',{@show_MBinFleder, fname});
			handles = guidata(hLine);
			uicontrol('Parent',h, 'Position',[10 10 130 23], 'String','Open in mbedit', ...
				'Callback',{@show_inMbedit, fname, handles.figure1});
			uicontrol('Parent',h, 'Position',[230 10 60 23], 'String','Cancel', 'Callback','delete(gcf)');
		else
			msgbox(sprintf('%s',LineInfo),'This line info')
		end
	else
		msgbox('Could not find reference to this object','This line info')
	end
	function show_MBinFleder(obj, evt, fname)
		mirone(['-Cshow_mb,guidata(gcf),' fname], '-Xpush_showPC')	% This case is more convoluted because it takes a fname in
	function show_inMbedit(obj, evt, fname, hMirFig)
		fcomm = ['mbedit -I' fname];
		if (isunix)
			resp = unix(fcomm);
			if (resp == 0)
				errordlg('Need to have MB-system installed and a X11 server installed and running.','Error')
			end
		else
			handles = guidata(hMirFig);
			if (handles.IamCompiled),	fcomm = ['start /b ' fcomm];	end
			s = dos(fcomm);
			if (s == 0)
				errordlg('Need to have MB-system installed and a X11 server installed and running.','Error')
			end
		end

% -----------------------------------------------------------------------------------------
function gmtfile_Info(obj,eventdata,h,data)
	agency = [];
	if (ischar(data) && exist(data, 'file') == 2)		% MGD77+ files transmit their names in 'data'
		[data, agency] = aux_funs('mgd77info',data);
	end
	str{1} = ['N_recs = ' num2str(data(1)) ', N_grav = ' num2str(data(2)) ', N_mag = ' num2str(data(3)) ...
		', N_top = ' num2str(data(4))];
	str{2} = ['E: = ' num2str(data(5)) '  W: = ' num2str(data(6))];
	str{3} = ['S: = ' num2str(data(7)) '  N: = ' num2str(data(8))];
	str{4} = ['Start day,month,year: = ' num2str(data(9)) '  ' num2str(data(10)) '  ' num2str(data(11))];
	str{5} = ['End   day,month,year: = ' num2str(data(12)) '  ' num2str(data(13)) '  ' num2str(data(14))];
	if (~isempty(agency)),		str{6} = '';	str{7} = agency;	end
	tag = get(h,'Tag');
	msgbox(str,tag)

% -----------------------------------------------------------------------------------------
function PB_All_Info(obj,eventdata,h,data)
	i = get(gco,'Userdata');
	txt_id = [];    txt_class = [];
	switch char(data(i).pb_id)
		case {'EU-NA','NA-EU'},         txt_id = 'Eurasia-North America';
		case {'AF-NA','NA-AF'},         txt_id = 'Africa-North America';
		case 'EU-AF',                   txt_id = 'Eurasia-Africa';
		case {'AF-AN','AN-AF'},         txt_id = 'Africa-Antartica';
	end

	switch char(data(i).class)      % Make Type-of-Boundary text
		case 'OTF',		txt_class = 'Oceanic Transform Fault';
		case 'OSF',		txt_class = 'Oceanic Spreadin Ridge';
		case 'CRB',		txt_class = 'Continental Rift Boundary';
		case 'CTF',		txt_class = 'Continental Transform Fault';
		case 'CCB',		txt_class = 'Continental Convergent Boundary';
		case 'OCB',		txt_class = 'Oceanic Convergent Boundary';
		case 'SUB',		txt_class = 'Subduction Zone';
	end

	if isempty(txt_id),		txt_id = char(data(i).pb_id);		end		% If id was not decoded, print id
	if isempty(txt_class),	txt_class = char(data(i).class);	end		% Shouldn't happen, but just in case

	msgbox( sprintf(['Plate pairs:           ' txt_id '\n' 'Boundary Type:    ' txt_class '\n' ...
			'Speed (mm/a):       ' sprintf('%g',data(i).vel) '\n' ...
			'Speed Azimuth:      ' sprintf('%g',data(i).azim_vel)] ),'Segment info')

% -----------------------------------------------------------------------------------------
function deleteObj(hTesoura)
% hTesoura is the handle to the 'Tesoura' uitoggletool
% Build the scisors pointer (this was done with the help of an image file)
	pointer = ones(16)*NaN;
	pointer(2,7) = 1;		pointer(2,11) = 1;		pointer(3,7) = 1;		pointer(3,11) = 1;
	pointer(4,7) = 1;		pointer(4,11) = 1;		pointer(5,7) = 1;		pointer(5,8) = 1;
	pointer(5,10) = 1;		pointer(5,11) = 1;		pointer(6,8) = 1;		pointer(6,10) = 1;
	pointer(7,8) = 1;		pointer(7,9) = 1;		pointer(7,10) = 1;		pointer(8,9) = 1;
	pointer(9,8) = 1;		pointer(9,9) = 1;		pointer(9,10) = 1;		pointer(10,8) = 1;
	pointer(10,10) = 1;		pointer(10,11) = 1;		pointer(10,12) = 1;		pointer(11,6) = 1;
	pointer(11,7) = 1;		pointer(11,8) = 1;		pointer(11,10) = 1;		pointer(11,13) = 1;
	pointer(12,5) = 1;		pointer(12,8) = 1;		pointer(12,10) = 1;		pointer(12,13) = 1;
	pointer(13,5) = 1;		pointer(13,8) = 1;		pointer(13,10) = 1;		pointer(13,13) = 1;
	pointer(14,5) = 1;		pointer(14,8) = 1;		pointer(14,11) = 1;		pointer(14,12) = 1;
	pointer(15,6:7) = 1;

	hFig = get(get(hTesoura,'Parent'),'Parent');
	state = uisuspend_j(hFig);
	set(hFig,'Pointer','custom','PointerShapeCData',pointer,'PointerShapeHotSpot',[1 8],...
		'WindowButtonDownFcn',{@wbd_delObj,hFig,hTesoura,state})

function wbd_delObj(obj,event,hFig,hTesoura,state)
	stype = get(hFig,'selectiontype');
	if (stype(1) == 'a')                    % A right click ('alt'), end killing
		uirestore_j(state, 'nochildren');	% Restore the figure's initial state
		set(hTesoura,'State','off')         % Set the Toggle button state to depressed
		return
	end
	h = gco;
	obj_type = get(h,'Type');
	if (strcmp(obj_type,'line') || strcmp(obj_type,'text') || strcmp(obj_type,'patch'))
		del_line([],[],h);
	end
	if (strcmp(obj_type,'text'))
		refresh;    % because of the text elements bug
	end

% -----------------------------------------------------------------------------------------
function del_line(obj,eventdata,h)
% Delete a line (or patch or text obj) but before check if it's in edit mode
	h = gco;    % I have to do this otherwise copied objects will have their del fun applyed to parent handle
	if (~isempty(getappdata(h,'polygon_data')))
		s = getappdata(h,'polygon_data');
		if strcmpi(s.controls,'on')     % Object is in edit mode, so this
			ui_edit_polygon(h)          % call will force out of edit mode
		end
	end
	delete(h);

% -----------------------------------------------------------------------------------------
function del_insideRect(obj, evt, h)
% Delete all lines/patches/text objects that have at least one vertex inside the rectangle

	[hLP, hText] = fish_inside_rect(h);
	delete(hLP);	delete(hText);
	if (~isempty(hText)),	refresh;    end     % Bloody text bug
	return

% -----------------------------------------------------------------------------------------
function trim_withPolygon(obj, evt, h, side)
% Trim lines with the polygon whose handle is H. Can trim inside our outside.
% H is the handle to polyg
% SIDE -> == 1, kill inside, else kill outside

	hL = fish_inside_rect(h, 'line');
	if (isempty(hL)),	return,		end

	hAxes = get(h,'Parent');		hFig = get(hAxes,'Parent');
	x0 = get(h,'XData');			y0 = get(h,'YData');	x0 = x0(:);		y0 = y0(:);

	inside = (side ~= 1);
	for (k = 1:numel(hL))
		x = get(hL(k),'XData');		y = get(hL(k),'YData');
		if (numel(x) < 2),	continue,	end		% Points are not for this game
		[xc, yc, iout] = intersections(x, y, x0, y0, 1);
		if (isempty(xc)),	continue,	end		% It didn't cross the polygon
		[iout, I] = sort(iout);			% Shit but they aren't always sorted
		if (any(diff(I) ~= 1))			% When shit happens we need to reorder the crossings too
			xc = xc(I);		yc = yc(I);
		end
		x_(numel(x) + numel(xc)) = single(0);	y_(numel(x_)) = single(0);
		pos_int = [0; ceil(iout)];			% Position where to add the intersection points
		for (n = 1:numel(iout)-1)
			pos_int(n+2) = pos_int(n+2)+n;
		end
		s = 1;
		for (n = 1:numel(xc))
			ind = pos_int(n)+1:pos_int(n+1)-1;
			x_(ind) = x(s:s+numel(ind)-1);		x_(pos_int(n+1)) = xc(n);
			y_(ind) = y(s:s+numel(ind)-1);		y_(pos_int(n+1)) = yc(n);
			s = s + numel(ind);
		end
		if (pos_int(end) < numel(x_))			% We still have the rest of the line after last intercept
			x_(pos_int(end)+1:end) = x(s:end);
			y_(pos_int(end)+1:end) = y(s:end);
		end

		in_out = inpolygon([x_(1) x_(end)],[y_(1) y_(end)], x0,y0);
		start_outside = ~in_out(1);				% True if line starts outside polyg
		end_outside = ~in_out(2);				% True if line ends outside polyg
		pos_int(1) = [];						% Remove the 0 indice. Not needed anymore
		if (start_outside)
			get_in = pos_int(1:2:end);			% Entry points in polyg.
			get_out = pos_int(2:2:end);			% Exit points out of polyg
		else
			get_out = pos_int(1:2:end);			% Exit points out of polyg.
			get_in  = pos_int(2:2:end);			% Entry points in polyg
		end
		if (inside)					% Construct inside segments
			if (start_outside)
				segs = cell(numel(get_in),2);
				if (isempty(get_out)),	get_out = numel(x_);	end
				for (n = 1:numel(get_in))
					segs{n,1} = x_(get_in(n):get_out(n));
					segs{n,2} = y_(get_in(n):get_out(n));
				end
			else
				segs = cell(1+numel(get_in),2);
				segs{1,1} = x_(1:get_out(1));	segs{1,2} = y_(1:get_out(1));
				if (~end_outside),	get_out(end+1) = numel(x_);	end		% For the algho
				for (n = 1:numel(get_in))
					segs{n+1,1} = x_(get_in(n):get_out(n+1));					
					segs{n+1,2} = y_(get_in(n):get_out(n+1));					
				end				
			end
		else						% Construct outside segments
			if (end_outside),	get_in(end+1) = numel(x_);	end		% For the algho
			if (start_outside)
				segs = cell(numel(get_in),2);
				segs{1,1} = x_(1:get_in(1));	segs{1,2} = y_(1:get_in(1));
				for (n = 1:numel(get_out))
					segs{n+1,1} = x_(get_out(n):get_in(n+1));
					segs{n+1,2} = y_(get_out(n):get_in(n+1));
				end
			else
				segs = cell(numel(get_out),2);
				for (n = 1:numel(get_out))
					segs{n,1} = x_(get_out(n):get_in(n));
					segs{n,2} = y_(get_out(n):get_in(n));
				end				
			end			
		end
		clear x_ y_				% They are reinialized if more lines
		break_trimeds(hL(k), segs, hFig, hAxes)
	end

	function break_trimeds(h, segs, hFig, hAxes)
	% Create the segment lines after the cell matix SEGS containing their coordinates
	set(h, 'XData',segs{1,1}, 'YData',segs{1,2})
	if (size(segs,1) == 1),	return,		end		% We are done, only one segment
	hCopy = copy_line_object([], h, hFig, hAxes);
	set(hCopy, 'XData',[], 'YData',[])			% So that we will not copy potentially big lines for nothing
	for (k = 2:size(segs,1))					% From second to before last segments
		hNext = copy_line_object([], hCopy, hFig, hAxes);
		set(hNext, 'XData',segs{k,1}, 'YData',segs{k,2})
	end
	delete(hCopy)

% -------------------------------------------------------------------------------------------------------
function [hLP, hText] = fish_inside_rect(h, opt)
% Get the handles of those objects that have at least one vertex inside the rectangle whose handle is H
% OPT is used to restrict what object types to search for. It can be
%	'all'   -> Search for 'lines', 'patches' and 'texts'
%	'line'  -> Search only for line type objects
%	'patch' -> Search only for patch type objects
%	'text'  -> Search only for text type objects

	if (nargin == 1),	opt = 'all';	end
	hLines = [];	hPatch = [];	hText = [];

	s = getappdata(h,'polygon_data');
	if (~isempty(s))            % If the rectangle is in edit mode, force it out of edit
		if strcmpi(s.controls,'on'),    ui_edit_polygon(h);     end
	end
	set(h, 'HandleVis','off')           % Make the rectangle handle invisible
	hAxes = get(h,'Parent');

	if (strcmpi(opt, 'all') || strcmpi(opt, 'line'))
		hLines = findobj(hAxes,'Type','line');     % Fish all objects of type line in Mirone figure
	end
	if (strcmpi(opt, 'all') || strcmpi(opt, 'patch'))
		hPatch = findobj(hAxes,'Type','patch');
	end
	if (strcmpi(opt, 'all') || strcmpi(opt, 'text'))
		hText = findobj(hAxes,'Type','text');
	end
	hLP = [hLines(:); hPatch(:)];
	rx = get(h,'XData');        ry = get(h,'YData');
	rx = [min(rx) max(rx)];     ry = [min(ry) max(ry)];
	found = false;
	for (i = 1:numel(hLP))		% Loop over objects to find if any is on edit mode
		s = getappdata(hLP(i),'polygon_data');
		if (~isempty(s))
			if strcmpi(s.controls,'on')     % Object is in edit mode, so this
				ui_edit_polygon(hLP(i))     % call will force out of edit mode
				found = true;
			end
		end
	end
	if (found)		% We have to do it again because some line handles have meanwhile desapeared
		hLines = findobj(hAxes,'Type','line');
		hPatch = findobj(hAxes,'Type','patch');
		hLP = [hLines(:); hPatch(:)];
	end
	set(h, 'HandleVis','on')			% Make the rectangle handle findable again

	f_LP = false(numel(hLP),1);
	for (i = 1:numel(hLP))				% Loop over objects to find out which cross the rectangle
		x = get(hLP(i),'XData');        y = get(hLP(i),'YData');
		if (any( (x >= rx(1) & x <= rx(2)) & (y >= ry(1) & y <= ry(2)) ))
			f_LP(i) = true;
		end
	end
	hTryAgain = hLP(~f_LP);				% Copy these for a second round
	hLP = hLP(f_LP);					% Retain only those that fulfill the condition

	if (~isempty(hTryAgain))			% Check for those beasts that cross but have no vertices inside rect
		f_LP = false(numel(hTryAgain),1);
		for (i = 1:numel(hTryAgain))	% Loop over potential slippers that don't have a vertex inside the rect
			x = get(hTryAgain(i),'XData');		y = get(hTryAgain(i),'YData');
			if (isempty(x)),	f_LP(i) = true;	continue,	end		% A stray empty line. Kill it too
			x0 = min(x);	x1 = max(x);		y0 = min(y);	y1 = max(y);
			rect = aux_funs('rectangle_and', [rx ry], [x0 x1 y0 y1]);
			if (~isempty(rect)),	f_LP(i) = true;		end
		end
		if (any(f_LP))					% OK, these potentially cross but only a more detailed check can tell
			hTryAgain = hTryAgain(f_LP);
			rxx = get(h,'XData');		ryy = get(h,'YData');
			f_LP = false(numel(hTryAgain),1);
			for (i = 1:numel(hTryAgain))% No escap than to compute intersections and see if they are empty or not
				x = get(hTryAgain(i),'XData');	y = get(hTryAgain(i),'YData');
				if (isempty(x)),	f_LP(i) = true;	continue,	end		% A stray empty line. Kill it too
				x0 = intersections(rxx, ryy, x, y, 0);
				if (~isempty(x0)),	f_LP(i) = true;		end
			end
			hTryAgain = hTryAgain(f_LP);
			hLP = [hLP; hTryAgain(:)];
		end
	end

	f_Txt = false(numel(hText),1);
	for (i = 1:numel(hText))	% Text objs are a bit different, so treat them separately
		pos = get(hText(i),'Position');
		if ( (pos(1) >= rx(1) && pos(1) <= rx(2)) && (pos(2) >= ry(1) && pos(2) <= ry(2)) )
			f_Txt(i) = true;
		end
	end
	hText = hText(f_Txt);

% -------------------------------------------------------------------------------------------------------
function changeAxesLabels(opt)
% This function formats the axes labels strings using a geographical notation

	hFig = get(0,'CurrentFigure');      hAxes = get(hFig,'CurrentAxes');
	x_tick = getappdata(hAxes,'XTickOrig');
	if (isa(x_tick, 'cell')),	x_tick = char(x_tick);	end		% Damn TMW never stops breaking compatibility (now in R2016a)
	y_tick = getappdata(hAxes,'YTickOrig');
	if (isa(y_tick, 'cell')),	y_tick = char(y_tick);	end
	n_xtick = size(x_tick,1);	n_ytick = size(y_tick,1);
	sep = ':';
	switch opt
		case 'ToDegDec'
			% This is easy because original Labels where saved in appdata
			str_x = x_tick;			str_y = y_tick;
			setappdata(hAxes,'LabelFormatType','DegDec')       % Save it so zoom can know the label type
		case 'ToDegMin'
			x_str = degree2dms(str2num(ddewhite(x_tick)),'DDMM',0,'str');     % x_str is a structure with string fields
			y_str = degree2dms(str2num(ddewhite(y_tick)),'DDMM',0,'str');
			str_x = [x_str.dd repmat(sep,n_xtick,1) x_str.mm];
			str_y = [y_str.dd repmat(sep,n_ytick,1) y_str.mm];
			setappdata(hAxes,'LabelFormatType','DegMin')        % Save it so zoom can know the label type
		case 'ToDegMinDec'
			x_str = degree2dms(str2num(ddewhite(x_tick)),'DDMM.x',2,'str');    % x_str is a structure with string fields
			y_str = degree2dms(str2num(ddewhite(y_tick)),'DDMM.x',2,'str');
			str_x = [x_str.dd repmat(sep,n_xtick,1) x_str.mm];
			str_y = [y_str.dd repmat(sep,n_ytick,1) y_str.mm];
			setappdata(hAxes,'LabelFormatType','DegMinDec')     % Save it so zoom can know the label type
		case 'ToDegMinSec'
			x_str = degree2dms(str2num(ddewhite(x_tick)),'DDMMSS',0,'str');    % x_str is a structure with string fields
			y_str = degree2dms(str2num(ddewhite(y_tick)),'DDMMSS',0,'str');
			str_x = [x_str.dd repmat(sep,n_xtick,1) x_str.mm repmat(sep,n_xtick,1) x_str.ss];
			str_y = [y_str.dd repmat(sep,n_ytick,1) y_str.mm repmat(sep,n_ytick,1) y_str.ss];
			setappdata(hAxes,'LabelFormatType','DegMinSec')      % Save it so zoom can know the label type
		case 'ToDegMinSecDec'
			x_str = degree2dms(str2num(ddewhite(x_tick)),'DDMMSS.x',1,'str');   % x_str is a structure with string fields
			y_str = degree2dms(str2num(ddewhite(y_tick)),'DDMMSS.x',1,'str');
			str_x = [x_str.dd repmat(sep,n_xtick,1) x_str.mm repmat(sep,n_xtick,1) x_str.ss];
			str_y = [y_str.dd repmat(sep,n_ytick,1) y_str.mm repmat(sep,n_ytick,1) y_str.ss];
			setappdata(hAxes,'LabelFormatType','DegMinSecDec')   % Save it so zoom can know the label type
	end
	set(hAxes,'XTickLabel',str_x, 'XTick', getappdata(hAxes,'XTickOrigNum'));
	set(hAxes,'YTickLabel',str_y, 'YTick', getappdata(hAxes,'YTickOrigNum'))

% -----------------------------------------------------------------------------------------
function sout = ddewhite(s)
%DDEWHITE Double dewhite. Strip both leading and trailing whitespace.
%
%   DDEWHITE(S) removes leading and trailing white space and any null characters
%   from the string S.  A null character is one that has an absolute value of 0.

%   Author:      Peter J. Acklam
%   Time-stamp:  2002-03-03 13:45:06 +0100
	
	if (nargin == 0),	error('DDWHITE ERROR: No input transmitted.');     end
	if ~ischar(s),		error('DDWHITE ERROR: Input must be a string (char array).');     end
	if isempty(s),		sout = s;		return,		end
	
	[r, c] = find(~isspace(s));
	if (size(s, 1) == 1),	sout = s(min(c):max(c));
	else,					sout = s(:,min(c):max(c));
	end

% --------------------------------------------------------------------
function set_transparency(obj,eventdata, h_patch)
% Sets the transparency of a patch object

	if (nargin == 2),	h_patch = gco;	end
	p_fc = get(h_patch,'FaceColor');
	if (strcmpi(p_fc,'none'))
		msg{1} = 'Transparency assumes that the element has a color';
		msg{2} = 'However, as you will agree, that is not the case.';
		msg{3} = 'See "Fill color" in this element properties.';
		warndlg(msg,'Warning')
		return
	end

	handles = guidata(get(0,'CurrentFigure'));
	r_mode = get(handles.figure1,'RendererMode');
	if (~strcmp(r_mode,'auto'))
		set(handles.figure1,'RendererMode','auto')
	end

	set(handles.figure1,'doublebuffer','on')        % I may be wrong, but I think patches are full of bugs

	% Define height of the bar code
	width = 7.0;            % Figure width in cm
	height = 1.5;           % Figure height in cm

	% Create figure for transparency display
	F = figure('Units','centimeters',...
		'Position',[1 10 [width height]],...
		'Toolbar','none', 'Menubar','none',...
		'Numbertitle','off',...
		'Name','Transparency',...
		'RendererMode','auto',...
		'Visible','of',...
		'Color',[.75 .75 .75]);

	T = uicontrol('style','text','string','Transparency  ',...
		'fontweight','bold','horizontalalignment','left',...
		'units','normalized','pos',[0.05  0.2  0.7  0.25],...
		'backgroundcolor',[.75 .75 .75]);

	transp = get(h_patch,'FaceAlpha');     % Get the previous transparency value

	pos=[0.02 0.5 .97 .25];
	S = {@apply_transparency,T,h_patch,handles};
	uicontrol('style','slider','units','normalized','position',pos, 'Callback',S,'min',0,'max',1,'Value',transp);

	set(F,'Visible','on')

% -----------------------------------------------------------------------------------------
function apply_transparency(obj,eventdata,T,h_patch,handles)
	val = get(obj,'Value');
	set(h_patch,'FaceAlpha',val)
	if (val > 0.99)
		h_all = findobj(handles.figure1,'Type','patch');
		set_painters = 1;
		for i = 1:length(h_all)
			if (get(h_all(i),'FaceAlpha') < 0.99)
				set_painters = 0;
			end
		end
		if (set_painters)
			set(handles.figure1,'Renderer','painters', 'RendererMode','auto')
		end
	end
	set(T,'String', sprintf('Opacity = %.2f',val))

% -----------------------------------------------------------------------------------------
function set_telhas_uicontext(h)
% h is a handle to a telhas patch object

	handles = guidata(h(1));	cmenuHand = uicontextmenu('Parent',handles.figure1);
	set(h, 'UIContextMenu', cmenuHand);
	cb_LineWidth = uictx_LineWidth(h);      % there are 5 cb_LineWidth outputs
	cb_solid = 'set(gco, ''LineStyle'', ''-''); refresh';
	cb_dashed = 'set(gco, ''LineStyle'', ''--''); refresh';
	cb_dotted = 'set(gco, ''LineStyle'', '':''); refresh';
	cb_dashdot = 'set(gco, ''LineStyle'', ''-.''); refresh';

	uimenu(cmenuHand, 'Label', 'Delete', 'Callback', 'delete(gco)');
	item_lw = uimenu(cmenuHand, 'Label', 'Line Width', 'Sep','on');
	setLineWidth(item_lw,cb_LineWidth)
	item_ls = uimenu(cmenuHand, 'Label', 'Line Style');
	setLineStyle(item_ls,{cb_solid cb_dashed cb_dotted cb_dashdot})
	item7 = uimenu(cmenuHand, 'Label', 'Line Color');
	cb_color = uictx_color(h,'EdgeColor');      % there are 9 cb_color outputs
	setLineColor(item7,cb_color)

	set_stack_order(cmenuHand)      % Change order in the stackpot

	uimenu(item7, 'Label', 'None', 'Sep','on', 'Callback', 'set(gco, ''EdgeColor'', ''none'');refresh');
	item8 = uimenu(cmenuHand, 'Label','Fill Color', 'Sep','on');
	cb_color = uictx_color(h,'facecolor');      % there are 9 cb_color outputs
	setLineColor(item8,cb_color)
	uimenu(item8, 'Label', 'None', 'Sep','on', 'Callback', 'set(gco, ''FaceColor'', ''none'');refresh');
	uimenu(cmenuHand, 'Label', 'Transparency', 'Callback', @set_transparency);

% -----------------------------------------------------------------------------------------
function set_stack_order(cmenuHand)
% Change order in the stackpot. cmenuHand is what it says. 
	item_order = uimenu(cmenuHand, 'Label', 'Order');
	uimenu(item_order, 'Label', 'Bring to Top', 'Callback','uistack_j(gco,''top'')');
	uimenu(item_order, 'Label', 'Send to Bottom', 'Callback','uistack_j(gco,''bottom'')');
	uimenu(item_order, 'Label', 'Move up', 'Callback','uistack_j(gco,''up'')');
	uimenu(item_order, 'Label', 'Move down', 'Callback','uistack_j(gco,''down'')');

%--------------------------------------------------------------------------------
function fig = getParentFigure(fig)
% Get the Fig handle of graphical object whose handle is FIG
	while ~isempty(fig) && ~strcmp('figure', get(fig,'Type'))
		fig = get(fig,'Parent');
	end

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

% $Id: draw_funs.m 7942 2016-09-05 22:45:44Z j $

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
		cmenuHand = uicontextmenu( 'parent',get(get(hand,'parent'),'parent') );		% We know 'hand' is a text handle
		set(hand, 'UIContextMenu', cmenuHand);
		cb_color = uictx_color(hand);	% there are 9 cb_color outputs
		uimenu(cmenuHand, 'Label', 'Change Font', 'Call', @text_FontSize);
		item_fc = uimenu(cmenuHand, 'Label', 'Font Color');
		setLineColor(item_fc,cb_color)
		uimenu(cmenuHand, 'Label', 'Edit   text', 'Call', 'set(gco, ''Editing'', ''on''); refresh', 'Sep','on');
		uimenu(cmenuHand, 'Label', 'Copy   text', 'Call', @copy_text_object);
		uimenu(cmenuHand, 'Label', 'Delete text', 'Call', 'delete(gco); refresh');
		uimenu(cmenuHand, 'Label', 'Move   text', 'Call', @move_text);
		uimenu(cmenuHand, 'Label', 'Rotate text', 'Call', @rotate_text);
		uimenu(cmenuHand, 'Label', 'Export text', 'Call', @export_text);
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
	else					hAx = get(get(0,'CurrentFigure'), 'CurrentAxes');
	end
	x = get(hLine, 'xdata');	y = get(hLine, 'ydata');
	if (strcmp(get(hLine,'type'), 'line'))
		h = line('xdata',x, 'ydata',y, 'Parent', hAx, 'LineWidth', get(hLine,'LineWidth'), ...
			'LineStyle',get(hLine,'LineStyle'), 'Color',get(hLine,'Color'), 'Tag',get(hLine,'Tag') );
		marker = get(hLine, 'Marker');				markSize = get(hLine,'MarkerSize');
		markFC = get(hLine,'MarkerFaceColor');		markEC = get(hLine,'MarkerEdgeColor');
		% Set the markers only if they are different from the red square markers used for edition
		% However, this is a risky test because if one of them is ever changed the test will fail.
		if ( ~( strcmp(marker, 'square') && strcmp(markFC, 'none') && markSize == 5 && isequal(markEC,[1 0 0]) ) )
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
% 		uimenu(item, 'Label', labels{k}, 'Call', cbs{k});
% 	end

% -----------------------------------------------------------------------------------------
function setLineStyle(item,cbs)
% Set the line Style uicontexts of graphic elements
	uimenu(item, 'Label', 'solid', 'Call', cbs{1});
	uimenu(item, 'Label', 'dashed', 'Call', cbs{2});
	uimenu(item, 'Label', 'dotted', 'Call', cbs{3});
	uimenu(item, 'Label', 'dash-dotted', 'Call', cbs{4});

% -----------------------------------------------------------------------------------------
function setLineColor(item,cbs)
% Set the line color uicontexts of graphic elements
	uimenu(item, 'Label', 'Black', 'Call', cbs{1});
	uimenu(item, 'Label', 'White', 'Call', cbs{2});
	uimenu(item, 'Label', 'Red', 'Call', cbs{3});
	uimenu(item, 'Label', 'Green', 'Call', cbs{4});
	uimenu(item, 'Label', 'Blue', 'Call', cbs{5});
	uimenu(item, 'Label', 'Yellow', 'Call', cbs{6});
	uimenu(item, 'Label', 'Cyan', 'Call', cbs{7});
	uimenu(item, 'Label', 'Magenta', 'Call', cbs{8});
	uimenu(item, 'Label', 'Other...', 'Call', cbs{9});

% -----------------------------------------------------------------------------------------
function setLineWidth(item,cbs)
% Set the line color uicontexts of graphic elements
	uimenu(item, 'Label', '1       pt', 'Call', cbs{1});
	uimenu(item, 'Label', '2       pt', 'Call', cbs{2});
	uimenu(item, 'Label', '3       pt', 'Call', cbs{3});
	uimenu(item, 'Label', '4       pt', 'Call', cbs{4});
	uimenu(item, 'Label', 'Other...', 'Call', cbs{5});

% -----------------------------------------------------------------------------------------
function setSHPuictx(h,opt)
% h is a handle to a shape line object

	handles = guidata(h(1));
	for (i = 1:numel(h))
		cmenuHand = uicontextmenu('Parent',handles.figure1);      set(h(i), 'UIContextMenu', cmenuHand);
		uimenu(cmenuHand, 'Label', 'Save line', 'Call', {@save_formated,h});
		uimenu(cmenuHand, 'Label', 'Delete this line', 'Call', {@del_line,h(i)});
		uimenu(cmenuHand, 'Label', 'Delete class', 'Call', 'delete(findobj(''Tag'',''SHPpolyline''))');

		cb_solid   = 'set(gco, ''LineStyle'', ''-''); refresh';
		cb_dashed  = 'set(gco, ''LineStyle'', ''--''); refresh';
		cb_dotted  = 'set(gco, ''LineStyle'', '':''); refresh';
		cb_dashdot = 'set(gco, ''LineStyle'', ''-.''); refresh';

		item = uimenu(cmenuHand, 'Label', 'Line Width', 'Sep','on');
		uimenu(item, 'Label', 'Other...', 'Call', {@other_LineWidth,h(i)});

		item = uimenu(cmenuHand, 'Label', 'Line Style');
		setLineStyle(item,{cb_solid cb_dashed cb_dotted cb_dashdot})

		item = uimenu(cmenuHand, 'Label', 'Line Color');
		uimenu(item, 'Label', 'Other...', 'Call', {@other_color,h(i)});
		ui_edit_polygon(h(i))
		isPt = getappdata(h(1), 'isPoint');
		if (~isempty(isPt) && ~isPt)	% For points it makes no sense a 'Join lines'
			uimenu(cmenuHand, 'Label', 'Join lines', 'Call', {@join_lines,handles.figure1});
		elseif (~isempty(get(h(i), 'UserData')))		% If we have z info
			uimenu(cmenuHand, 'Label', 'Quick grid', 'Call', {@shp_quick_grd,h(i)}, 'Sep', 'on');			
		end
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
	if (strcmp(get(h,'Type'),'patch')),
		IS_PATCH = true;
		if (IS_PATCH && ~LINE_ISCLOSED),	LINE_ISCLOSED = true;	end
	elseif (strcmp(get(h,'Tag'),'SeismicLine'))
		IS_SEISMICLINE = true;
	end

	handles = guidata(get(h,'Parent'));		% Get Mirone handles

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
		uimenu(cmenuHand, 'Label', 'Delete me', 'Call', {@del_line,h});
		uimenu(cmenuHand, 'Label', 'Delete inside rect', 'Call', {@del_insideRect,h});
		item_tools = uimenu(cmenuHand, 'Label', 'Trim with rect');
		uimenu(item_tools, 'Label', 'inside', 'Call', {@trim_withPolygon,h,1});
		uimenu(item_tools, 'Label', 'outside', 'Call', {@trim_withPolygon,h,0});
		ui_edit_polygon(h)
	elseif (IS_LINE)
		uimenu(cmenuHand, 'Label', 'Delete', 'Call', {@del_line,h});
		ui_edit_polygon(h)			% Set edition functions
	elseif (IS_MBTRACK)				% Multibeam tracks, when deleted, have to delete also the bars
		uimenu(cmenuHand, 'Label', 'Delete track (left-click on it)', 'Call', 'save_track_mb(1);');
		% Old style edit function. New edit is provided by ui_edit_polygon which doesn't work with mbtracks 
		uimenu(cmenuHand, 'Label', 'Edit track (left-click on it)', 'Call', 'edit_track_mb');
	end
	uimenu(cmenuHand, 'Label', label_save, 'Call', {@save_formated,h});
	if (~IS_SEISPOLYG && ~IS_MBTRACK && ~strcmp(get(h,'Tag'),'FaultTrace'))	% Those are not to allowed to copy
		if (~LINE_ISCLOSED && ~IS_ARROW)
			uimenu(cmenuHand, 'Label', 'Join lines', 'Call', {@join_lines,handles.figure1});
		end
		uimenu(cmenuHand, 'Label', 'Copy', 'Call', {@copy_line_object,handles.figure1,handles.axes1});
	end
	if (~IS_SEISPOLYG && ~IS_ARROW && ~IS_RECTANGLE)
		if (numel(x) > 2)
			uimenu(cmenuHand, 'Label', 'Spline Smooth', 'Call', {@smooth_line,h})
		end
		uimenu(cmenuHand, 'Label', label_length, 'Call', @show_LineLength)
		uimenu(cmenuHand, 'Label', label_azim,   'Call', @show_lineAzims)
	end
	if (IS_MBTRACK),	uimenu(cmenuHand, 'Label', 'All MB-tracks length', 'Call', @show_AllTrackLength);	end

	if (LINE_ISCLOSED)
		uimenu(cmenuHand, 'Label', 'Area under polygon', 'Call', @show_Area);
		if (IS_PATCH && ~IS_SEISPOLYG)
			item8 = uimenu(cmenuHand, 'Label','Fill Color');
			setLineColor( item8, uictx_color(h, 'facecolor') )		% there are 9 cb_color outputs
			uimenu(item8, 'Label', 'None', 'Sep','on', 'Call', 'set(gco, ''FaceColor'', ''none'');refresh');
			uimenu(cmenuHand, 'Label', 'Transparency', 'Call', @set_transparency);
		end
		uimenu(cmenuHand, 'Label', 'Create Mask', 'Call', 'poly2mask_fig(guidata(gcbo),gco)');
	end
	if (IS_RECTANGLE && handles.validGrid)
		uimenu(cmenuHand, 'Label', 'Make Chess board', 'Call', {@chessify, h});
	end

	if (handles.image_type ~= 20 && ~LINE_ISCLOSED && strcmp(opt,'line'))
		if ((ndims(get(handles.hImg,'CData')) == 2) || (handles.validGrid))	% Because Track of RGB doesn't know how to save
			if (handles.nLayers > 1)
				cbTrack = 'setappdata(gcf,''TrackThisLine'',gco); mirone(''ExtractProfile_CB'',guidata(gcbo),''3D'')';
				uimenu(cmenuHand, 'Label', '3D interpolation', 'Call', cbTrack);
			end
			cbTrack = 'setappdata(gcf,''TrackThisLine'',gco); mirone(''ExtractProfile_CB'',guidata(gcbo),''point'')';
			uimenu(cmenuHand, 'Label', 'Point interpolation', 'Call', cbTrack);
		end
		cbTrack = 'setappdata(gcf,''TrackThisLine'',gco); mirone(''ExtractProfile_CB'',guidata(gcbo))';
		uimenu(cmenuHand, 'Label', 'Extract profile', 'Call', cbTrack);
	end

	if strcmp(opt,'MBtrack'),	uimenu(cmenuHand, 'Label', 'Show track''s Swath Ratio', 'Call', {@show_swhatRatio,h});	end

	if (IS_RECTANGLE)
		uimenu(cmenuHand, 'Label', 'Rectangle limits (edit)', 'Sep','on', 'Call', @rectangle_limits);
		uimenu(cmenuHand, 'Label', 'Register Image', 'Call', @rectangle_register_img);
		if (~handles.validGrid)
			uimenu(cmenuHand, 'Label', 'Transplant Image here','Call', 'transplants(gco,''image'')');
		end
		if (handles.geog)
			uimenu(cmenuHand, 'Label', 'Grid/Image mosaicer', 'Call', 'mosaicer(gco)');
			uimenu(cmenuHand, 'Label', 'Get image from Web Map Server', 'Call', 'wms_tool(gco)');
			uimenu(cmenuHand, 'Label', 'CMT Catalog (Web download)', 'Call', 'globalcmt(gcf,gco)');
		end

		item_tools = uimenu(cmenuHand, 'Label', 'ROI Crop Tools','Sep','on');
		if (handles.validGrid)    % Option only available to recognized grids
			uimenu(item_tools, 'Label', 'Crop Grid', 'Call', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''CropaGrid_pure'')');
			uimenu(item_tools, 'Label', 'Crop Image', 'Call', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco)');
			uimenu(item_tools, 'Label', 'Crop Image (with coords)', 'Call', ...
				'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''CropaWithCoords'')');
			uimenu(item_tools, 'Label', 'Set to constant', 'Call', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''SetConst'')','Sep','on');
			item_tools3 = uimenu(item_tools, 'Label', 'Stats');
			uimenu(item_tools3,'Label', 'Mean',          'Call', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''GetMean'')');
			uimenu(item_tools3,'Label', 'Median',        'Call', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''GetMedian'')');
			uimenu(item_tools3,'Label', 'STD',           'Call', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''GetSTD'')');
			uimenu(item_tools, 'Label', 'Clip grid',     'Call', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''Clip'')');
			uimenu(item_tools, 'Label', 'Median filter', 'Call', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''MedianFilter'')');
			uimenu(item_tools, 'Label', 'Spline smooth', 'Call', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''SplineSmooth'')');
			uimenu(item_tools, 'Label', 'Histogram (grid)', 'Call', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''CropaGrid_histo'')');
			uimenu(item_tools, 'Label', 'Histogram (image)', 'Call', 'image_histo(guidata(gcbo),gco)');
%			uimenu(item_tools, 'Label', 'Detect Fronts', 'Call', 'cayula_cornillon(guidata(gcbo),gco)');
			uimenu(item_tools, 'Label', 'Power',           'Call', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''CropaGrid_power'')');
			uimenu(item_tools, 'Label', 'Autocorrelation', 'Call', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''CropaGrid_autocorr'')');
			uimenu(item_tools, 'Label', 'FFT tool', 'Call', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''CropaGrid_fftTools'')');
			item_fill = uimenu(item_tools, 'Label', 'Fill gaps');
			uimenu(item_fill, 'Label', 'Fill gaps (surface)','Call', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''FillGaps'',''surface'')');
			uimenu(item_fill, 'Label', 'Fill gaps (cubic)',  'Call', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''FillGaps'',''cubic'');');
			uimenu(item_fill, 'Label', 'Fill gaps (linear)', 'Call', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''FillGaps'',''linear'');');
		else			% We have an Image
			uimenu(item_tools, 'Label', 'Crop Image', 'Call', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco)');
			if (handles.image_type == 3)
					uimenu(item_tools, 'Label', 'Crop Image (with coords)', 'Call', ...
						'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''CropaWithCoords'')');
			end
			uimenu(item_tools, 'Label', 'Histogram', 'Call', 'image_histo(guidata(gcbo),gco)');
		end
		uimenu(item_tools, 'Label', 'Image process', 'Sep','on', 'Call', 'mirone(''DrawClosedPolygon_CB'',guidata(gcbo),gco)');
		deal_opts({'MGG' 'MICROLEV' 'GMT_DB_IDS' 'GMT_SYMBOL'}, cmenuHand);
	end

	if (~IS_SEISPOLYG && LINE_ISCLOSED && ~IS_RECTANGLE)
		item_tools2 = uimenu(cmenuHand, 'Label', 'ROI Crop Tools','Sep','on');
		if (handles.validGrid)    % Option only available to recognized grids
			uimenu(item_tools2, 'Label', 'Crop Grid', 'Call', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''CropaGrid_pure'')');
			uimenu(item_tools2, 'Label', 'Crop Image', 'Call', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco)');
			uimenu(item_tools2, 'Label', 'Crop Image (with coords)', 'Call', ...
				'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''CropaWithCoords'')');
			uimenu(item_tools2, 'Label', 'Set to constant', 'Call', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''ROI_SetConst'')','Sep','on');
			item_tools3 = uimenu(item_tools2, 'Label', 'Stats');
			uimenu(item_tools3, 'Label', 'Mean', 'Call', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''ROI_Mean'')');
			uimenu(item_tools3, 'Label', 'Median', 'Call', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''ROI_Median'')');
			uimenu(item_tools3, 'Label', 'STD', 'Call', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''ROI_STD'')');
			uimenu(item_tools2, 'Label', 'Clip grid', 'Call', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''ROI_Clip'')');
			uimenu(item_tools2, 'Label', 'Median filter', 'Call', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''ROI_MedianFilter'')');
			uimenu(item_tools2, 'Label', 'Spline smooth', 'Call', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''ROI_SplineSmooth'')');
			uimenu(item_tools2, 'Label', 'Histogram (grid)', 'Call', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''CropaGrid_histo'')');
			uimenu(item_tools2, 'Label', 'Histogram (image)','Call', 'image_histo(guidata(gcbo),gco)');
			hP = getappdata(handles.figure1, 'ParentFig');
			if ( ~isempty(hP) && ishandle(hP) && ~isempty(strfind(get(handles.figure1,'Name'), 'spectrum')) )
				uimenu(item_tools2, 'Label', 'Low Pass FFT filter', 'Call', 'mirone(''GridToolsSectrum_CB'',guidata(gcbo), ''lpass'', gco)');
				uimenu(item_tools2, 'Label', 'High Pass FFT filter','Call', 'mirone(''GridToolsSectrum_CB'',guidata(gcbo), ''hpass'', gco)');
			end
		else			% We have an Image
			uimenu(item_tools2, 'Label', 'Crop Image', 'Call', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco)');
			if (handles.image_type == 3)
					uimenu(item_tools2, 'Label', 'Crop Image (with coords)', 'Call', ...
						'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''CropaWithCoords'')');
			end
			uimenu(item_tools2, 'Label', 'Histogram', 'Call', 'image_histo(guidata(gcbo),gco)');
		end
		uimenu(item_tools2, 'Label', 'Image process', 'Sep','on', 'Call', 'mirone(''DrawClosedPolygon_CB'',guidata(gcbo),gco)');

		if (strcmp(get(h,'Tag'),'EulerTrapezium'))
			uimenu(cmenuHand, 'Label', 'Compute Euler Pole', 'Sep','on', 'Call',...
				'calc_bonin_euler_pole(get(gco,''XData''), get(gco,''YData''));' );
		end
		%uimenu(cmenuHand, 'Label', 'Testa patches', 'Sep','on', 'Call', 'patch_options(gco)');
	end

	setLineWidth(uimenu(cmenuHand, 'Label', 'Line Width', 'Sep','on'), cb_LineWidth)
	setLineStyle(uimenu(cmenuHand, 'Label', 'Line Style'), {cb_solid cb_dashed cb_dotted cb_dashdot})
	if (IS_PATCH),		cb_color = uictx_color(h,'EdgeColor');	end      % there are 9 cb_color outputs
	setLineColor(uimenu(cmenuHand, 'Label', 'Line Color'), cb_color)
	set_stack_order(cmenuHand)      % Set options to change order in the stackpot

	if (strcmp(get(h,'Tag'),'FaultTrace'))      % For Okada modeling
		uimenu(cmenuHand, 'Label', 'Okada', 'Sep','on', 'Call', {@okada_model,h,'okada'});    
		uimenu(cmenuHand, 'Label', 'Mansinha', 'Call', {@okada_model,h,'mansinha'});    
	end

	if (IS_SEISPOLYG)                         % Seismicity options
		% gco gives the same handle as h 
		uimenu(cmenuHand, 'Label', 'Save events', 'Call', 'save_seismicity(gcf,[],gco)', 'Sep','on');
		uimenu(cmenuHand, 'Label', 'Find clusters', 'Call', 'find_clusters(gcf,gco)');
		itemHist = uimenu(cmenuHand, 'Label','Histograms');
		uimenu(itemHist, 'Label', 'Guttenberg & Richter', 'Call', 'histos_seis(gco,''GR'')');
		uimenu(itemHist, 'Label', 'Cumulative number', 'Call', 'histos_seis(gco,''CH'')');
		uimenu(itemHist, 'Label', 'Cumulative moment', 'Call', 'histos_seis(gco,''CM'')');
		uimenu(itemHist, 'Label', 'Magnitude', 'Call', 'histos_seis(gco,''MH'')');
		uimenu(itemHist, 'Label', 'Time', 'Call', 'histos_seis(gco,''TH'')');
		uimenu(itemHist, 'Label', 'Display in Table', 'Call', 'histos_seis(gcf,''HT'')','Sep','on');
		%uimenu(itemHist, 'Label', 'Hour of day', 'Call', 'histos_seis(gco,''HH'')');
		itemTime = uimenu(cmenuHand, 'Label','Time series');
		uimenu(itemTime, 'Label', 'Time magnitude', 'Call', 'histos_seis(gco,''TM'')');
		uimenu(itemTime, 'Label', 'Time depth', 'Call', 'histos_seis(gco,''TD'')');
		uimenu(cmenuHand, 'Label', 'Mc and b estimate', 'Call', 'histos_seis(gco,''BV'')');
		uimenu(cmenuHand, 'Label', 'Fit Omori law', 'Call', 'histos_seis(gco,''OL'')');
		%uimenu(cmenuHand, 'Label', 'Skell', 'Call', 'esqueleto_tmp(gco)','Sep','on');
	elseif (IS_SEISMICLINE)
		uimenu(cmenuHand, 'Label', 'Set buffer zone', 'Call', {@seismic_line,h,'buf'}, 'Sep','on');
		uimenu(cmenuHand, 'Label', 'Project seismicity', 'Call', {@seismic_line,h,'proj'}, 'Enable','off');
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
			uimenu(cmenuHand, 'Label', 'Rectangle limits (edit)', 'Call', @rectangle_limits);
			ud = 1;		set(h(k),'UserData',ud)
		end
		if (handles.validGrid)			% Option only available to recognized grids
			uimenu(cmenuHand, 'Label', 'Crop Grid', 'Call', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''CropaGrid_pure'')');
		end
% 		uimenu(cmenuHand, 'Label', 'Adjust to nesting dimensions', 'Sep', 'on', 'Call', 'nesting_sizes(gcbo)');
		if (ud > 1 || handles.validGrid)
			uimenu(cmenuHand, 'Label', 'Show nesting info', 'Call', 'nesting_sizes(gcbo,''Info'')');
			uimenu(cmenuHand, 'Label', 'Create blank grid', 'Call', 'nesting_sizes(gcbo,''Grid'')');
		end
		if (k == numel(h))
			uimenu(cmenuHand, 'Label', 'New nested grid', 'Call', 'nesting_sizes(gcbo,''New'')');
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
	else			cb_color = uictx_color(h);
	end
	setLineWidth(uimenu(cmenuHand, 'Label', 'Line Width', 'Sep','on'), cb_LineWidth)
	setLineStyle(uimenu(cmenuHand, 'Label', 'Line Style'), {cb_solid cb_dashed cb_dotted cb_dashdot})
	setLineColor(uimenu(cmenuHand, 'Label', 'Line Color'), cb_color)

% -----------------------------------------------------------------------------------------
function seismic_line(obj,evt,hL,opt)
% Plot seismicity projected along a polyline and enclosed inside a buffer zone of that pline
	handles = guidata(hL);
	hP = get(hL, 'UserData');	% handle to the buffer zone (if not exist, will be created later)
	if (opt(1) == 'b')		% Create or expand a buffer zone
		resp = inputdlg('Width of buffer zone (deg)','Buffer width',[1 30],{'0.5'});
		if isempty(resp),	return,		end
		[y, x] = buffer_j(get(hL, 'ydata'), get(hL, 'xdata'), str2double(resp{1}), 'out', 13, 1);
		if (isempty(x)),	return,		end
		if (isempty(hP) || ~ishandle(hP))
			hP = patch('XData',x, 'YData',y, 'Parent',handles.axes1, 'EdgeColor',handles.DefLineColor, ...
				'FaceColor','none', 'LineWidth',handles.DefLineThick+1, 'Tag','SeismicBuffer');
			uistack_j(hP,'bottom'),		draw_funs(hP,'line_uicontext')
			set(hL, 'UserData', hP)			% Save it there to ease its later retrival
		else
			set(hP,'XData',x, 'YData',y)	% Buffer already existed, just resize it
		end
		h = findobj(get(obj,'Parent'),'Label','Project seismicity');
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
		ptFrac = ptFrac(:,5) + 1;								% +1 because mapproject is 0 based
		intSeg = fix(ptFrac);
		rd = lineSegDistAcum(intSeg) + (ptFrac - intSeg) .* lineSegDist(intSeg);
		[rd, ind] = sort(rd);

		evt_time = evt_time(ind);	evt_dep = evt_dep(ind);
		figure;
		h = subplot(2,2,1);		plot(rd, evt_time, '.'),	set(h,'ylim',[min(evt_time) max(evt_time)],'xlim',[rd(1) rd(end)]);
		h = subplot(2,2,2);		plot(rd, evt_dep, '.'),		set(h,'ylim',[min(evt_dep) max(evt_dep)+0.01],'xlim',[rd(1) rd(end)]);
		subplot(2,2,3);			histo_m('hist', rd, 0:25:rd(end), [0 rd(end)]);
		h = subplot(2,2,4);		plot(rd, evt_mag, '.'),		set(h,'ylim',[min(evt_mag) max(evt_mag)],'xlim',[rd(1) rd(end)]);
	end

% -----------------------------------------------------------------------------------------
function hCopy = copy_line_object(obj, evt, hFig, hAxes)
% Make a copy of an line object and update some of its properties.
% If HCOPY output arg is specified, return without setting WindowButtonMotionFcn(). In this case one very
% likely need to transmit the original object in EVT, since this function was not called after a mouse event.

	if (~isempty(evt)),		oldH = evt;		% When calling this function by another that just wants a copy.
	else					oldH = gco(hFig);
	end

	newH = copyobj(oldH,hAxes);
	h = findobj(get(newH,'uicontextmenu'),'label','Save line');
	if (~isempty(h))        % Replace the old line handle in the 'Save line' Callback by the just created one
		hFun = get(h,'Call');
		hFun{2} = newH;
		set(h,'Call',hFun)
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
function wbd_MovePolygon(obj,eventdata,h,state)
	uirestore_j(state, 'nochildren');	% Restore the figure's initial state
	ui_edit_polygon(h)					% Reset the edition functions with the correct handle
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
		[x, y, was_closed] = join2lines([hCurrLine hLines(k)]);
		if (~was_closed),	delete(hLines(k)),	end		% Closed polygons are ignored 
		set(hCurrLine, 'XData',x, 'YData',y)
	end

% ---------
function [x, y, was_closed] = join2lines(hLines, TOL)
% Joint the two lines which have handles "hLines" by their closest connection points
% TOL is the max distance that the two lines can be apart and still be joint.
%	If not provided, defaults to Inf, but if it is and the two lines are too further
%	apart than X,Y & WAS_CLOSED are all set to empty. It's callers responsability to check it.

	x1 = get(hLines(1),'XData');		y1 = get(hLines(1),'YData');
	x2 = get(hLines(2),'XData');		y2 = get(hLines(2),'YData');
	if (nargin == 1),	TOL = inf;		end
	
	was_closed = false;
	if ( (x2(1) == x2(end)) && (y2(1) == y2(end)) )		% Ignore closed polygons
		x = x1;		y = y1;
		was_closed = true;
		return
	end

	% Find how segments should be glued. That is find the closest extremities
	dif_x = [(x1(1) - x2(1)); (x1(1) - x2(end)); (x1(end) - x2(1)); (x1(end) - x2(end))];
	dif_y = [(y1(1) - y2(1)); (y1(1) - y2(end)); (y1(end) - y2(1)); (y1(end) - y2(end))];
	dist = sum([dif_x dif_y] .^2 ,2);	% Square of distances between the 4 extremities
	[mimi, I] = min(dist);				% We only care about the min location
	if (mimi > TOL)
		x = [];		y = [];		was_closed = [];
		return
	end
	if (I == 1)				% Line 1 starts near the begining of line 2
		last = 1;
		if (x1(1) == x2(1) && y1(1) == y2(1)),		last = 2;		end
		x = [x2(end:-1:last) x1];	y = [y2(end:-1:last) y1];
	elseif (I == 2)			% Line 2 ends near the begining of line 1 
		x = [x2 x1];				y = [y2 y1];
	elseif (I == 3)			% Line 1 ends near the begining of line 2
		x = [x1 x2];				y = [y1 y2];
	else					% Line 1 ends near the end of line 2
		last = 1;
		if (x1(end) == x2(end) && y1(end) == y2(end)),	last = 2;	end		% If points are equal, do not repeat it
		x = [x1 x2(end:-1:last)];	y = [y1 y2(end:-1:last)];
	end
% -----------------------------------------------------------------------------------------
	
% --------------------------------------------------------------------
function shp_quick_grd(hObj, evt, h)
% Automatically calculate a grid from a PointZ shapefile data
	x = get(h, 'XData');	y = get(h, 'YData');
	z = double(get(h, 'UserData'));
	x_min = min(x);			x_max = max(x);
	y_min = min(y);			y_max = max(y);
	x_inc = abs(median(diff(x(1:100))));
	y_inc = abs(median(diff(y(1:100))));
	if (y_inc < x_inc/20),		y_inc = x_inc;		% Quite likely if data was previously gridded abd dumped.
	elseif (x_inc < x_inc/20),	x_inc = y_inc;
	end
	opt_I = sprintf('-I%.8g/%.8g', x_inc, y_inc);
	if (isa(x, 'single'))
		x = double(x);		y = double(y);
	end
	nx = round((x_max - x_min) / x_inc) + 1;
	ny = round((y_max - y_min) / y_inc) + 1;
	opt_R = sprintf('-R%.12g/%.12g/%.12g/%.12g', x_min, x_min + (nx-1)*x_inc, y_min, y_min + (ny-1)*y_inc);
	[Z, head] = gmtmbgrid_m(x, y, z, opt_I, opt_R, '-Mz', '-C2');
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
		else							lc = 'k';
		end
		if (isfield(struc, 'thick')),	lThick = struc.thick;
		else							lThick = 1;
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
		if min(size(x))==1, n=sqrt(numel(x)); m=n; else [m,n]=size(x); end
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

	if ( isempty(hQuiver) || ~ishandle(hQuiver(1)) )		% No arrows yet.
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
		uimenu(cmenuHand, 'Label', 'Save line', 'Call', @save_line);
		uimenu(cmenuHand, 'Label', 'Delete', 'Call', 'delete(gco)');
		if (handles.validGrid)
			item_ct = uimenu(cmenuHand, 'Label', 'ROI Crop Tools','Sep','on');
			uimenu(item_ct, 'Label', 'Crop Grid', 'Call', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''CropaGrid_pure'')');
			uimenu(item_ct, 'Label', 'Set to constant', 'Call', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''ROI_SetConst'')');
			uimenu(item_ct, 'Label', 'Mean', 'Call', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''ROI_Mean'')');
			uimenu(item_ct, 'Label', 'Clip grid', 'Call', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''ROI_Clip'')');
			uimenu(item_ct, 'Label', 'Median filter', 'Call', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''ROI_MedianFilter'')');
			uimenu(item_ct, 'Label', 'Spline smooth', 'Call', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''ROI_SplineSmooth'')');
			uimenu(item_ct, 'Label', 'Histogram', 'Call', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''CropaGrid_histo'')');
		end
		cb_LineWidth = uictx_LineWidth(h(i));			% there are 5 cb_LineWidth outputs
		item_lw = uimenu(cmenuHand, 'Label', 'Line Width', 'Sep','on');
		uimenu(item_lw, 'Label', '1     pt', 'Call', cb_LineWidth{1});
		uimenu(item_lw, 'Label', 'Other...', 'Call', cb_LineWidth{5});
		item8 = uimenu(cmenuHand, 'Label','Fill Color', 'Sep','on');
		cb_color = uictx_color(h(i),'facecolor');		% there are 9 cb_color outputs
		uimenu(item8, 'Label', 'Other...', 'Call', cb_color{9});
		uimenu(item8, 'Label', 'None', 'Call', 'set(gco, ''FaceColor'', ''none'');refresh');
		uimenu(cmenuHand, 'Label', 'Transparency', 'Call', @set_transparency);
		uimenu(cmenuHand, 'Label', 'Create Mask', 'Call', 'poly2mask_fig(guidata(gcbo),gco)');
		if (handles.validGrid)
		end
		if (handles.image_type ~= 20)
			uimenu(cmenuHand, 'Label', 'Region-Of-Interest', 'Sep','on', 'Call', ...
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
		else                                        min_len = 5000; % Assume meters
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
	uimenu(cmenuHand, 'Label', 'Delete', 'Call', 'delete(gco)');
	cb_Fill_surface = 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''FillGaps'',''surface'');delete(gco)';
	cb_Fill_cubic = 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''FillGaps'',''cubic'');delete(gco)';
	cb_Fill_linear = 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''FillGaps'',''linear'');delete(gco)';
	cb_Fill_sea   = 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''FillGaps'',''sea'');delete(gco)';
	uimenu(cmenuHand, 'Label', 'Fill gaps (surface)', 'Call', cb_Fill_surface);
	uimenu(cmenuHand, 'Label', 'Fill gaps (cubic)', 'Call', cb_Fill_cubic);
	uimenu(cmenuHand, 'Label', 'Fill gaps (linear)', 'Call', cb_Fill_linear);
	uimenu(cmenuHand, 'Label', 'Fill gaps (sea)', 'Call', cb_Fill_sea);

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

	uimenu(cmenuHand, 'Label', 'Delete contour', 'Call',{@remove_singleContour,h});
	uimenu(cmenuHand, 'Label', 'Delete all contours', 'Call', cb_rac);
	% item1 = uimenu(cmenuHand, 'Label', 'Edit contour (left-click on it)', 'Call', cb1);
	uimenu(cmenuHand, 'Label', 'Join contours', 'Call', {@join_lines,handles.figure1});
	uimenu(cmenuHand, 'Label', 'Save contour', 'Call', {@save_formated,h});
	uimenu(cmenuHand, 'Label', 'Contour length', 'Call', {@show_LineLength,[]});
	uimenu(cmenuHand, 'Label', 'Area under contour', 'Call', @show_Area);
	item_lw = uimenu(cmenuHand, 'Label', 'Contour Line Width', 'Sep','on');
	setLineWidth(item_lw,cb_LineWidth)
	item_ls = uimenu(cmenuHand, 'Label', 'Contour Line Style');
	setLineStyle(item_ls,{cb18 cb19 cb20 cb21})
	item_lc = uimenu(cmenuHand, 'Label', 'Contour Line Color');
	setLineColor(item_lc,cb_color)
	cb_CLineWidth = uictx_Class_LineWidth(h);           % there are 5 cb_CLineWidth outputs
	item8 = uimenu(cmenuHand, 'Label', 'All Contours Line Width', 'Sep','on');
	uimenu(item8, 'Label', '1       pt', 'Call', cb_CLineWidth{1});
	uimenu(item8, 'Label', '2       pt', 'Call', cb_CLineWidth{2});
	uimenu(item8, 'Label', '3       pt', 'Call', cb_CLineWidth{3});
	uimenu(item8, 'Label', 'Other...', 'Call', cb_CLineWidth{5});
	cb_CLineStyle = uictx_Class_LineStyle(h);        % there are 4 cb_CLineStyle outputs
	item9 = uimenu(cmenuHand, 'Label', 'All Contours Line Style');
	uimenu(item9, 'Label', 'solid', 'Call', cb_CLineStyle{1});
	uimenu(item9, 'Label', 'dashed', 'Call', cb_CLineStyle{2});
	uimenu(item9, 'Label', 'dotted', 'Call', cb_CLineStyle{3});
	uimenu(item9, 'Label', 'dash-dotted', 'Call', cb_CLineStyle{4});
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
	
	uimenu(cmenuHand, 'Label', label, 'Call', 'delete(gco)');
	uimenu(cmenuHand, 'Label', 'Save coastline', 'Call', {@save_formated,h});

	if (handles.validGrid && strcmp(tag,'CoastLineNetCDF'))		% Options to apply ocean/land masking
		item = uimenu(cmenuHand, 'Label', 'Mask', 'Sep','on');
		uimenu(item, 'Label', 'Land',   'Call', {@apply_grdlandMask,h, 'L'})
		uimenu(item, 'Label', 'Oecean', 'Call', {@apply_grdlandMask,h, 'O'})
	end

	setLineWidth(uimenu(cmenuHand, 'Label', 'Line Width', 'Sep','on'), cb_LineWidth)
	item_ls = uimenu(cmenuHand, 'Label', 'Line Style');
	setLineStyle(item_ls,{cb13 cb14 cb15 cb16})
	setLineColor(uimenu(cmenuHand, 'Label', 'Line Color'), cb_color)
	ui_edit_polygon(h)

% -----------------------------------------------------------------------------------------
function apply_grdlandMask(hObj, evt, h, opt)
% Blank current grid with the default's grdlandmask setting for this coastline resolution
	handles = guidata(h);
	opt_R = sprintf('-R%.10g/%.10g/%.10g/%.10g',handles.head(1), handles.head(2), handles.head(3), handles.head(4));
	opt_I = sprintf('-I%0.10g/%0.10g',handles.head(8), handles.head(9));
	opt_F = ' ';	opt_e = ' ';
	if (handles.head(7)),		opt_F = '-F';	end
	if (handles.IamCompiled),	opt_e = '-e';	end
	opt_D = sprintf('-D%s', getappdata(h,'resolution'));	% Get the resolution as stored in line's appdata
	opt_N = '0/1/0/0/0';		% Continent masking
	mask = c_grdlandmask(opt_R, opt_I, opt_D, opt_N, opt_F, opt_e, '-A0/0/1', '-V');
	if (opt == 'O'),	mask = ~mask;	end			% Ocean masking. Setting '1/0/0/0/0' is not working. A bug.

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

	[X,Y,Z] = load_grd(handles);	Z(mask) = NaN;
	setappdata(handles.figure1,'dem_z',Z);
	handles.have_nans = grdutils(Z,'-N');
	zz = grdutils(Z,'-L');			handles.head(5:6) = [zz(1) zz(2)];
	handles.firstIllum = true;
	set(handles.haveNaNs,'Vis','on')
	guidata(handles.figure1, handles)
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
	uimenu(cmenuHand, 'Label', 'Segment info', 'Call', {@PB_All_Info,h_cur,data_cur});
	uimenu(cmenuHand, 'Label', 'Delete class', 'Call', 'delete(findobj(''Tag'',''PB_All''))', 'Sep','on');
	uimenu(cmenuHand, 'Label', 'Segment length', 'Call', {@show_LineLength,[]});
	item3 = uimenu(cmenuHand, 'Label', 'Line Width', 'Sep','on');
	uimenu(item3, 'Label', '2       pt', 'Call', cb_LineWidth{2});
	uimenu(item3, 'Label', '3       pt', 'Call', cb_LineWidth{3});
	uimenu(item3, 'Label', '4       pt', 'Call', cb_LineWidth{4});
	uimenu(item3, 'Label', 'Other...', 'Call', cb_LineWidth{5});
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
	cmenuHand = uicontextmenu('Parent',handles.figure1);
	set(h, 'UIContextMenu', cmenuHand);
	cb_LineWidth = uictx_LineWidth(h);		% there are 5 cb_LineWidth outputs
	cb_color = uictx_color(h);				% there are 9 cb_color outputs
	cbls1 = 'set(gco, ''LineStyle'', ''-''); refresh';   cbls2 = 'set(gco, ''LineStyle'', ''--''); refresh';
	cbls3 = 'set(gco, ''LineStyle'', '':''); refresh';   cbls4 = 'set(gco, ''LineStyle'', ''-.''); refresh';
	if (~isempty(data) && ~all(isempty(cat(2,data{:}))) )
		uimenu(cmenuHand, 'Label', [tag ' info'], 'Call', {@Isochrons_Info,data});
		uimenu(cmenuHand, 'Label', ['Delete this ' tag ' line'], 'Call', {@del_line,h}, 'Sep','on');
	else
		uimenu(cmenuHand, 'Label', ['Delete this ' tag ' line'], 'Call', {@del_line,h});
	end
	uimenu(cmenuHand, 'Label', ['Delete all ' tag ' lines'], 'Call', {@remove_symbolClass,h});
	uimenu(cmenuHand, 'Label', ['Save this '  tag ' line'],  'Call', @save_line);
	uimenu(cmenuHand, 'Label', ['Save all '   tag ' lines'], 'Call', {@save_line,h});
	uimenu(cmenuHand, 'Label', 'Join lines', 'Call', {@join_lines,handles.figure1});
	uimenu(cmenuHand, 'Label', 'Line azimuths', 'Call', @show_lineAzims);
	uimenu(cmenuHand, 'Label', 'Line length', 'Call', {@show_LineLength,[],'nikles'});
	LINE_ISCLOSED = 0;
	for (i = 1:numel(h))
		x = get(h(i),'XData');      y = get(h(i),'YData');
		if (numel(x) > 2 && (x(1) == x(end)) && (y(1) == y(end)))		% See if we have at least one closed line
			LINE_ISCLOSED = 1;		break
		end
	end
	% If at least one is closed, activate the Area option
	if (LINE_ISCLOSED)
		uimenu(cmenuHand, 'Label', 'Area under polygon', 'Call', @show_Area);
		uimenu(cmenuHand, 'Label', 'Create Mask', 'Call', 'poly2mask_fig(guidata(gcbo),gco)');
	end
	if (handles.image_type ~= 20 && (ndims(get(handles.hImg,'CData')) == 2 || handles.validGrid))
		if (handles.nLayers > 1)
			cbTrack = 'setappdata(gcf,''TrackThisLine'',gco); mirone(''ExtractProfile_CB'',guidata(gcbo),''3D'')';
			uimenu(cmenuHand, 'Label', '3D interpolation', 'Call', cbTrack);
		end
		cbTrack = 'setappdata(gcf,''TrackThisLine'',gco); mirone(''ExtractProfile_CB'',guidata(gcbo),''point'')';
		uimenu(cmenuHand, 'Label', 'Point interpolation', 'Call', cbTrack);
		cbTrack = 'setappdata(gcf,''TrackThisLine'',gco); mirone(''ExtractProfile_CB'',guidata(gcbo))';
		uimenu(cmenuHand, 'Label', 'Extract profile', 'Call', cbTrack);
	end
	setLineWidth(uimenu(cmenuHand, 'Label', 'Line Width', 'Sep','on'), cb_LineWidth)
	setLineStyle(uimenu(cmenuHand, 'Label', 'Line Style'), {cbls1 cbls2 cbls3 cbls4})
	setLineColor(uimenu(cmenuHand, 'Label', 'Color'), cb_color)
	set_stack_order(cmenuHand)		% Set options to change order in the stackpot	
	% --------- Now set the class properties
	cb_ClassColor = uictx_Class_LineColor(h);        % there are 9 cb_color outputs
	setLineColor(uimenu(cmenuHand, 'Label', ['All ' tag ' Color'], 'Sep','on'), cb_ClassColor)
	cb_ClassLineWidth = uictx_Class_LineWidth(h);    % there are 5 cb_ClassLineWidth outputs
	item_Class_lw = uimenu(cmenuHand, 'Label', ['All ' tag ' Line Width']);
	uimenu(item_Class_lw, 'Label', '1       pt', 'Call', cb_ClassLineWidth{1});
	uimenu(item_Class_lw, 'Label', '2       pt', 'Call', cb_ClassLineWidth{2});
	uimenu(item_Class_lw, 'Label', '3       pt', 'Call', cb_ClassLineWidth{3});
	uimenu(item_Class_lw, 'Label', '4       pt', 'Call', cb_ClassLineWidth{4});
	uimenu(item_Class_lw, 'Label', 'Other...', 'Call', cb_ClassLineWidth{5});
	cb_ClassLineStyle = uictx_Class_LineStyle(h);    % there are 4 cb_ClassLineStyle outputs
	item_Class_lt = uimenu(cmenuHand, 'Label', ['All ' tag ' Line Style']);
	setLineStyle(item_Class_lt,{cb_ClassLineStyle{1} cb_ClassLineStyle{2} cb_ClassLineStyle{3} cb_ClassLineStyle{4}})
	item = uimenu(cmenuHand, 'Label', 'Compute pole to neighbor', 'Sep','on');
	%uimenu(item, 'Label', 'Bonin (all in this plate)', 'Sep','on', 'Call', {@pole2neighbor, [], 'Bonin'});
	uimenu(item, 'Label', 'Best-Fit (all in this plate)', 'Call', {@pole2neighbor, [], 'anglefit'});
	uimenu(item, 'Label', 'Best-Fit (only me)', 'Call', {@pole2neighbor, [], 'anglefit', 1});
	uimenu(item, 'Label', 'Best-Fit (only me -> iterate)', 'Call', {@pole2neighbor, [], 'anglefit', 10});
	uimenu(item, 'Label', 'Show Results', 'Call', {@pole2neighbor, [], 'showresults'});
	uimenu(cmenuHand, 'Label', 'Age Grid', 'Call', {@pole2neighbor, [], 'agegrid'});
	uimenu(cmenuHand, 'Label', 'Make Age-script', 'Call', @make_age_script);
	uimenu(cmenuHand, 'Label', 'Euler rotation', 'Sep','on', 'Call', 'euler_stuff(gcf,gco)');
	for (i=1:length(h)),		ui_edit_polygon(h(i)),		end		% Set edition functions

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
	uimenu(cmenuHand, 'Label', [tag ' info'], 'Call', {@gmtfile_Info,h,data});
	uimenu(cmenuHand, 'Label', ['Delete this ' tag ' line'], 'Call', 'delete(gco)', 'Sep','on');
	uimenu(cmenuHand, 'Label', ['Save this ' tag ' line'], 'Call', @save_line);
	uimenu(cmenuHand, 'Label', 'Open with gmtedit', 'Call', {@call_gmtedit,h});
 	uimenu(cmenuHand, 'Label', 'Extract Mag Chunk', 'Call', {@call_gmtedit,h,'nikles'});
	uimenu(cmenuHand, 'Label', 'Create Mask', 'Call', 'poly2mask_fig(guidata(gcbo),gco)');
	deal_opts('mgg_coe', cmenuHand);
	%uimenu(cmenuHand, 'Label', 'Try to relocate', 'Call', {@tryRelocate,h});
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
	if (nargin == 4)		% Call helper window to extract a chunk of the mag anom profile
		hFig = get(get(h,'Parent'),'Parent');
		[xp,yp] = getline_j(hFig);
		if (numel(xp) < 2),		return,		end
		hGL = line('XData', xp, 'YData', yp,'Color','y','Parent',get(h,'Parent'),'LineWidth',3,'Tag','polyline');
		guidelineAzim = azimuth_geo(yp(1), xp(1), yp(end), xp(end));
		mag_synthetic(hFig, h, hGL, guidelineAzim)
		return
	end
	pt = get(get(h,'Parent'), 'CurrentPoint');
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
		if (length(c) > 1),         % That is, if a color was selected
			if isempty(prop)
				set(h,'color',c);   refresh;
			else
				set(h,prop,c);   refresh;
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
	
	uimenu(cmenuHand, 'Label', 'Delete', 'Call', 'delete(gco)');
	uimenu(cmenuHand, 'Label', 'Save line', 'Call', {@save_formated,h});
	uimenu(cmenuHand, 'Label', 'Line length', 'Call', {@show_LineLength,[],'total'});
	uimenu(cmenuHand, 'Label', 'Line azimuth','Call', @show_lineAzims)
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
	
	uimenu(cmenuHand, 'Label', 'Delete', 'Call', 'delete(gco)')
	uimenu(cmenuHand, 'Label', 'Save circle', 'Call', {@save_formated,h})
	uimenu(cmenuHand, 'Label', 'Line length', 'Call', {@show_LineLength,[]})
	if ~strcmp(tag,'CircleEuler')       % "Just" a regular geographical circle
		uimenu(cmenuHand, 'Label', 'Make me Donut', 'Sep','on', 'Call', {@donutify, h(1), 'geog'})
		uimenu(cmenuHand, 'Label', 'Region-Of-Interest', 'Sep','on', 'Call', cb_roi)
	else
		uimenu(cmenuHand, 'Label', 'Compute velocity', 'Sep','on', 'Call', {@report_EulerVel,h})
		uimenu(cmenuHand, 'Label', 'Plot tangent arrow', 'Call', {@report_EulerVel, h, 'tg'})
		uimenu(cmenuHand, 'Label', 'Plot normal arrow',  'Call', {@report_EulerVel, h, 'no'})
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
	uimenu(cmenuHand, 'Label', 'Delete', 'Call', 'delete(gco)')
	uimenu(cmenuHand, 'Label', 'Radial average', 'Call', 'grid_profiler(gcf, gco)', 'Sep', 'on')

% -----------------------------------------------------------------------------------------
function set_circleCart_uicontext(h)
% h is a handle to a circle (in cartesian coords) object
	handles = guidata(h(1));	cmenuHand = uicontextmenu('Parent',handles.figure1);
	set(h, 'UIContextMenu', cmenuHand);
	cb_solid  = 'set(gco, ''LineStyle'', ''-''); refresh';   cb_dashed      = 'set(gco, ''LineStyle'', ''--''); refresh';
	cb_dotted = 'set(gco, ''LineStyle'', '':''); refresh';   cb_dash_dotted = 'set(gco, ''LineStyle'', ''-.''); refresh';
	cb_roi = 'mirone(''DrawClosedPolygon_CB'',guidata(gcbo),gco)';
	
	uimenu(cmenuHand, 'Label', 'Delete', 'Call', 'delete(gco)');
	uimenu(cmenuHand, 'Label', 'Save circle', 'Call', {@save_formated,h});
	uimenu(cmenuHand, 'Label', 'Circle perimeter', 'Call', {@show_LineLength,[]});
	uimenu(cmenuHand, 'Label', 'Move (interactive)', 'Call', {@move_circle,h});
	uimenu(uimenu(cmenuHand, 'Label', 'Change'), 'Label', 'By coordinates', 'Call', {@change_CircCenter1,h});
	uimenu(cmenuHand, 'Label', 'Make me Donut', 'Sep','on', 'Call', {@donutify, h(1), 'cart'})
	uimenu(cmenuHand, 'Label', 'Region-Of-Interest', 'Sep','on', 'Call', cb_roi);
	hp = getappdata(handles.figure1, 'ParentFig');
	if ( ~isempty(hp) && ishandle(hp) && ~isempty(strfind(get(handles.figure1,'Name'), 'spectrum')) )
		uimenu(cmenuHand, 'Label', 'Low Pass FFT filter', 'Call', 'mirone(''GridToolsSectrum_CB'',guidata(gcbo), ''lpass'', gco)');
		uimenu(cmenuHand, 'Label', 'High Pass FFT filter','Call', 'mirone(''GridToolsSectrum_CB'',guidata(gcbo), ''hpass'', gco)');
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

	D2R = pi/180;
	earth_rad = 6371e3;    % Earth radius in km
	plat = plat*D2R;      plon = plon*D2R;      omega = omega*D2R;
	alat = alat*D2R;      alon = alon*D2R;

	x = cos(plat).*sin(plon).*sin(alat) - cos(alat).*sin(alon).*sin(plat);    % East vel
	y = cos(alat).*cos(alon).*sin(plat) - cos(plat).*cos(plon).*sin(alat);    % North vel
	z = cos(plat).*cos(alat).*sin(alon-plon);
	vlon = -sin(alon).*x + cos(alon).*y;
	vlat = -sin(alat).*cos(alon).*x-sin(alat).*sin(alon).*y + cos(alat).*z;
	azim = 90 - atan2(vlat,vlon) / D2R;

	if (azim < 0),		azim = azim + 360;		end		% Give allways the result in the 0-360 range

	x = sin(alat).*sin(plat) + cos(alat).*cos(plat).*cos(plon-alon);
	delta = acos(x);
	vel = omega/1e+4 * earth_rad .* sin(delta);			% to give velocity in cm/Ma
	
	if (nargin == 6)
		azim = azim * D2R;			% Put it back in radians
		vel = vel .* cos(azim);		% X
		azim = vel .* sin(azim);	% Y
	end

% -----------------------------------------------------------------------------------------
function set_vector_uicontext(h)
% h is a handle to a vector object
	h = h(1);
	handles = guidata(h);	cmenuHand = uicontextmenu('Parent',handles.figure1);
	set(h, 'UIContextMenu', cmenuHand);
	uimenu(cmenuHand, 'Label', 'Delete', 'Call', 'delete(gco)');
	uimenu(cmenuHand, 'Label', 'Save line', 'Call', {@save_formated, h});
	uimenu(cmenuHand, 'Label', 'Copy', 'Call', {@copy_line_object, handles.figure1, handles.axes1});
 	uimenu(cmenuHand, 'Label', 'Edit (redraw)', 'Call', {@arrowRedraw, h});
	uimenu(cmenuHand, 'Label', 'Reshape (head)', 'Call', {@arrow_shape,h});
	if (handles.geog)		% No solution yet for cartesian arrows
		uimenu(cmenuHand, 'Label', 'Copy (mirror)', 'Call', {@mirror_arrow, h});
	end
	setLineColor(uimenu(cmenuHand, 'Label', 'Line Color','Sep','on'), uictx_color(h,'EdgeColor'))
	item2 = uimenu(cmenuHand, 'Label','Fill Color');
	setLineColor( item2, uictx_color(h, 'facecolor') )
	uimenu(item2, 'Label', 'None', 'Sep','on', 'Call', 'set(gco, ''FaceColor'', ''none'');refresh');
	uimenu(cmenuHand, 'Label', 'Transparency', 'Call', @set_transparency);
	set(h, 'ButtonDownFcn', 'move_obj(1)')

% -----------------------------------------------------------------------------------------
function arrowRedraw(obj, evt, hVec)
% Edit a previously existing arrow
	hAxes = get(hVec,'parent');
	hFig = get(hAxes,'parent');
	state = uisuspend_j(hFig);		% Remember initial figure state
	set(hFig,'Pointer', 'crosshair');
	hVec(2) = line('XData', [], 'YData', [],'LineWidth',0.5,'Tag','Arrow');
	ud = get(hVec(1),'UserData');
	vectorFirstButtonDown(hFig, hAxes, hVec, state, ud.anchors, ...
		ud.headLength, ud.vFac, ud.aspectRatio)

% -----------------------------------------------------------------------------------------
function mirror_arrow(obj, evt, h)
% Make a copy of the arrow whose handle is H, but mirrored. That is pointing in the oposite sense
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
function chessify(obj, evt, h)
...
	prompt = {'Enter DX patch size:','Enter DX patch size:','One Z value','Other Z value'};
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
	uimenu(item1, 'Label', '1       pt', 'Call', cb_LineWidth{1});
	uimenu(item1, 'Label', '2       pt', 'Call', cb_LineWidth{2});
	uimenu(item1, 'Label', '3       pt', 'Call', cb_LineWidth{3});
	uimenu(item1, 'Label', '4       pt', 'Call', cb_LineWidth{4});
	uimenu(item1, 'Label', '5       pt', 'Call', cb7);
	uimenu(item1, 'Label', 'Other...', 'Call', cb_LineWidth{5});
	item_ls = uimenu(cmenuHand, 'Label', 'Line Style');
	setLineStyle(item_ls,{cb10 cb11 cb12 cb13})
	item_lc = uimenu(cmenuHand, 'Label', 'Line Color');
	setLineColor(item_lc,cb_color)

% -----------------------------------------------------------------------------------------
function cb = uictx_color(h,opt)
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

	function other_color(obj,eventdata,h,opt)
	if (nargin == 3),   opt = [];   end
	c = uisetcolor;
	if (length(c) > 1),			% That is, if a color was selected
		if ~isempty(opt) && ischar(opt)
			set(h,opt,c);   refresh;
		else
			set(h,'color',c);   refresh;
		end
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

	if (nargin >= 5),		pt = anchor;
	else					pt = get(hAxes, 'CurrentPoint');
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
	x = get(h(2), 'XData');		y = get(h(2), 'YData');
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

	region = bg_region('empty');
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
	set(ax,'xlim',new_xlim,'ylim',new_ylim,'YDir','normal')
	handles.head(1:4) = [new_xlim new_ylim];
	resizetrue(handles, [], 'xy');
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
function copy_text_object(obj,eventdata)
    copyobj(gco,gca);
    move_text([],[])

% -----------------------------------------------------------------------------------------
function move_text(obj,eventdata)
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
	if (x<x_lim(1)) || (x>x_lim(2)) || (y<y_lim(1)) || (y>y_lim(2));   return; end
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
		uimenu(cmenuHand, 'Label', 'Hotspot info', 'Call', {@hotspot_info,h,data.name,data.age,[]});
		uimenu(cmenuHand, 'Label', 'Plot name', 'Call', {@hotspot_info,h,data.name,data.age,'text'});
		separator = true;
	elseif strcmp(tag,'volcano')    % DATA must be a structure containing name, description & dating fields
		uimenu(cmenuHand, 'Label', 'Volcano info', 'Call', {@volcano_info,h,data.name,data.desc,data.dating});
		separator = true;
	elseif strcmp(tag,'meteor')		% DATA must be a structure containing name, diameter & dating fields
		uimenu(cmenuHand, 'Label', 'Impact info', 'Call', {@meteor_info,h,data.name,data.diameter,data.dating,data.exposed,data.btype});
		separator = true;
	elseif strcmp(tag,'mar_online')	% DATA must be a structure containing name, code station & country fields
		uimenu(cmenuHand, 'Label', 'Download Mareg (2 days)', 'Call', {@mareg_online,h,data});
		uimenu(cmenuHand, 'Label', 'Download Mareg (Calendar)', 'Call', {@mareg_online,h,data,'cal'});
		separator = true;
	elseif strcmp(tag,'hydro')		% DATA must be a cell array with 5 cols contining description of each Vent
		uimenu(cmenuHand, 'Label', 'Hydrotermal info', 'Call', {@hydro_info,h,data});
		separator = true;	
	elseif ( strcmp(tag,'DSDP') || strcmp(tag,'ODP') || strcmp(tag,'IODP') )	% DATA is a struct with leg, site, z & penetration fields
		uimenu(cmenuHand, 'Label', [tag ' info'], 'Call', {@ODP_info,h,data.leg,data.site,data.z,data.penetration});
		separator = true;
	elseif strcmp(tag,'City_major') || strcmp(tag,'City_other')
		this_not = true;
	elseif strcmp(tag,'Earthquakes')	% DATA is empty because I didn't store any info (they are too many)
		seismicity_options = isappdata(h,'SeismicityTime');
	elseif strcmp(tag,'Pointpolyline')	% DATA is empty because it doesn't have any associated info
		s_info = getappdata(h, 'LineInfo');		% But there might be one in appdata
		if (~isempty(s_info))
			uimenu(cmenuHand, 'Label', 'Symbol info', 'Call', {@Isochrons_Info,s_info});
			sep = 'on';
		end
		this_not = false;
	elseif strcmp(tag,'TTT')			% DATA is empty
		this_not = false;
	elseif strcmp(tag,'TideStation')	% DATA is empty
		tide_options = true;
		separator = false;
	else
		this_not = false;
	end

	if (~this_not)						% non class symbols can be moved
		ui_edit_polygon(h)				% Set edition functions
		uimenu(cmenuHand, 'Label', 'Move (precise)', 'Call', {@change_SymbPos,h}, 'Sep', sep);
	end

	if separator
		if (~more_than_one)				% Single symbol
			uimenu(cmenuHand, 'Label', 'Remove', 'Call', 'delete(gco)', 'Sep','on');
		else							% Multiple symbols
			uimenu(cmenuHand, 'Label', 'Remove this', 'Call', {@remove_one_from_many,h}, 'Sep','on');
		end
	else
		if (~more_than_one)				% Single symbol
			uimenu(cmenuHand, 'Label', 'Remove', 'Call', 'delete(gco)');
		else							% Multiple symbols
			uimenu(cmenuHand, 'Label', 'Remove this', 'Call', {@remove_one_from_many,h});
			uimenu(cmenuHand, 'Label', 'Delete all', 'Call', {@del_line,h});
		end
	end

	if (this_not)						% individual symbols don't belong to a class
		uimenu(cmenuHand, 'Label', 'Remove class', 'Call', {@remove_symbolClass,h});
	end

	if strcmp(tag,'LinkedSymb')			% Symbols plotted by the "Linked" option in either 'ecran' or 'gmtedit'
		uimenu(cmenuHand, 'Label', 'Remove class', 'Call', {@remove_symbolClass,h});
		uimenu(cmenuHand, 'Label', 'Save all linked pts', 'Call', {@save_line,h});
		this_not = true;				% Set this to avoid going inside next IF-test
	end

	if (~this_not)						% class symbols don't export
		uimenu(cmenuHand, 'Label', 'Save', 'Call', {@export_symbol,h});
		if (handles.image_type ~= 20)
			if (strcmp(tag,'Pointpolyline'))	% Allow pure grdtrack interpolation
				cbTrack = 'setappdata(gcf,''TrackThisLine'',gco); mirone(''ExtractProfile_CB'',guidata(gcbo),''point'')';
				uimenu(cmenuHand, 'Label', 'Point interpolation', 'Call', cbTrack, 'Sep','on');
			end
			if (handles.nLayers > 1)
				cbTrack = 'setappdata(gcf,''TrackThisLine'',gco); mirone(''ExtractProfile_CB'',guidata(gcbo),''3D'')';
				uimenu(cmenuHand, 'Label', '3D interpolation', 'Call', cbTrack, 'Sep','on');
			end
		end
	end

	if (seismicity_options)
		uimenu(cmenuHand, 'Label', 'Save events', 'Call', 'save_seismicity(gcf,gco)', 'Sep','on');
		uimenu(cmenuHand, 'Label', 'Seismicity movie', 'Call', 'animate_seismicity(gcf,gco)');
		uimenu(cmenuHand, 'Label', 'Draw polygon', 'Call', ...
			'mirone(''DrawClosedPolygon_CB'',guidata(gcbo),''SeismicPolyg'')');
		uimenu(cmenuHand, 'Label', 'Draw seismic line', 'Call', ...
			'mirone(''DrawLine_CB'',guidata(gcbo),''SeismicLine'')');
		itemHist = uimenu(cmenuHand, 'Label','Histograms');
		uimenu(itemHist, 'Label', 'Guttenberg & Richter', 'Call', 'histos_seis(gco,''GR'')');
		uimenu(itemHist, 'Label', 'Cumulative number', 'Call', 'histos_seis(gco,''CH'')');
		uimenu(itemHist, 'Label', 'Cumulative moment', 'Call', 'histos_seis(gco,''CM'')');
		uimenu(itemHist, 'Label', 'Magnitude', 'Call', 'histos_seis(gco,''MH'')');
		uimenu(itemHist, 'Label', 'Time', 'Call', 'histos_seis(gco,''TH'')');
		uimenu(itemHist, 'Label', 'Display in Table', 'Call', 'histos_seis(gcf,''HT'')','Sep','on');
		%uimenu(itemHist, 'Label', 'Hour of day', 'Call', 'histos_seis(gco,''HM'')');
		itemTime = uimenu(cmenuHand, 'Label','Time series');
		uimenu(itemTime, 'Label', 'Time magnitude', 'Call', 'histos_seis(gco,''TM'')');
		uimenu(itemTime, 'Label', 'Time depth', 'Call', 'histos_seis(gco,''TD'')');
		uimenu(cmenuHand, 'Label', 'Mc and b estimate', 'Call', 'histos_seis(gcf,''BV'')');
		uimenu(cmenuHand, 'Label', 'Fit Omori law', 'Call', 'histos_seis(gcf,''OL'')');
	end

	if (tide_options)
		uimenu(cmenuHand, 'Label', 'Plot tides', 'Call', {@tidesStuff,h,'plot'}, 'Sep','on');
		uimenu(cmenuHand, 'Label', 'Station Info', 'Call', {@tidesStuff,h,'info'});
		%uimenu(cmenuHand, 'Label', 'Tide Calendar', 'Call', {@tidesStuff,h,'calendar'});
	end

	itemSymb = uimenu(cmenuHand, 'Label', 'Symbol', 'Sep','on');
	cb_mark = uictx_setMarker(h,'Marker');              % there are 13 uictx_setMarker outputs
	uimenu(itemSymb, 'Label', 'plus sign', 'Call', cb_mark{1});
	uimenu(itemSymb, 'Label', 'circle', 'Call', cb_mark{2});
	uimenu(itemSymb, 'Label', 'asterisk', 'Call', cb_mark{3});
	uimenu(itemSymb, 'Label', 'point', 'Call', cb_mark{4});
	uimenu(itemSymb, 'Label', 'cross', 'Call', cb_mark{5});
	uimenu(itemSymb, 'Label', 'square', 'Call', cb_mark{6});
	uimenu(itemSymb, 'Label', 'diamond', 'Call', cb_mark{7});
	uimenu(itemSymb, 'Label', 'upward triangle', 'Call', cb_mark{8});
	uimenu(itemSymb, 'Label', 'downward triangle', 'Call', cb_mark{9});
	uimenu(itemSymb, 'Label', 'right triangle', 'Call', cb_mark{10});
	uimenu(itemSymb, 'Label', 'left triangle', 'Call', cb_mark{11});
	uimenu(itemSymb, 'Label', 'five-point star', 'Call', cb_mark{12});
	uimenu(itemSymb, 'Label', 'six-point star', 'Call', cb_mark{13});
	itemSize = uimenu(cmenuHand, 'Label', 'Size');
	cb_markSize = uictx_setMarker(h,'MarkerSize');              % there are 8 uictx_setMarker outputs
	uimenu(itemSize, 'Label', '7       pt', 'Call', cb_markSize{1});
	uimenu(itemSize, 'Label', '8       pt', 'Call', cb_markSize{2});
	uimenu(itemSize, 'Label', '9       pt', 'Call', cb_markSize{3});
	uimenu(itemSize, 'Label', '10     pt', 'Call', cb_markSize{4});
	uimenu(itemSize, 'Label', '12     pt', 'Call', cb_markSize{5});
	uimenu(itemSize, 'Label', '14     pt', 'Call', cb_markSize{6});
	uimenu(itemSize, 'Label', '16     pt', 'Call', cb_markSize{7});
	uimenu(itemSize, 'Label', 'other...', 'Call', cb_markSize{8});
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
function remove_one_from_many(obj,eventdata,h)
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
	if (strcmp(get(h, 'Tag'),'isochron'))	% Special case
		name = create_isoc_name(h);
		if (isempty(name)),		name = '.dat';	end		% Something did not work well
	end
	str1 = {'*.dat;*.DAT', 'Line file (*.dat,*.DAT)'; '*.*', 'All Files (*.*)'};
	[FileName,PathName] = put_or_get_file(handles,str1,'Select Line File name','put', name);
	if isequal(FileName,0),			return,		end
	x = get(h,'XData');				y = get(h,'YData');

	fname = [PathName FileName];
	[PATH,FNAME,EXT] = fileparts([PathName FileName]);
	if isempty(EXT),    fname = [PathName FNAME '.dat'];    end

	if (strcmp(get(h, 'Tag'),'GMT_DBpolyline'))
		save_GMT_DB_asc(h, fname)
		return
	end

	fid = fopen(fname, 'w');
	if (fid < 0),	errordlg(['Can''t open file:  ' fname],'Error'),	return,		end

	if (strcmp(get(h, 'Tag'),'LinkedSymb'))		% Though those are rare cases, we need to check it first
		hAll = findobj(handles.axes1,'-depth',1,'type','line','Tag','LinkedSymb');
		for (k = 1:numel(hAll))
			x = get(hAll(k),'XData');		y = get(hAll(k),'YData');	str = getappdata(hAll(k),'box');
			if (isempty(str))
				fprintf(fid,'%.6f\t%.6f\n',[x(:)'; y(:)']);
			else
				fprintf(fid,'%.6f\t%.6f\t%s\n',[x(:)'; y(:)'], str);
			end
		end
	elseif (~iscell(x))
		LineInfo = getappdata(h,'LineInfo');
		if (~isempty(LineInfo)),	fprintf(fid,'> %s\n',LineInfo);	end
		fprintf(fid,'%.6f\t%.6f\n',[x(:)'; y(:)']);
	else
		for (i = 1:numel(h))
			LineInfo = getappdata(h(i),'LineInfo');
			if (~isempty(LineInfo))
				str = ['> ' LineInfo];
			else
				str = '>';
			end
			fprintf(fid,'%s\n',str);
			fprintf(fid,'%.6f\t%.6f\n',[x{i}(:)'; y{i}(:)']);
		end
	end
	fclose(fid);

function name = create_isoc_name(h)
% Generate the likely name that I use to use for the isochrons
	try
		LineInfo = getappdata(h,'LineInfo');
		[isoc, r] = strtok(LineInfo);
		[plates, r] = strtok(r);		r = ddewhite(r);
		ind = strfind(plates, '/');
		if (strcmp(plates, 'NORTH') || strcmp(plates, 'SOUTH'))		% First word of first plate is this
			ind = strfind(r, '/');
			name = sprintf('c%s_%s%s_%s.dat',isoc, plates(1), r(1), r(ind+1:ind+2));
		elseif (strcmp(plates(ind+1:end), 'NORTH') || strcmp(plates(ind+1:end), 'SOUTH'))	% Second plate is this
			name = sprintf('c%s_%s_%s%s.dat',isoc, plates(1:2), plates(ind+1), r(1));
		else
			name = sprintf('c%s_%s_%s.dat',isoc, plates(1:2), plates(ind+1:ind+2));
		end
	catch
		name = [];
	end

function save_GMT_DB_asc(h, fname)
% Go through all GMT_DB polygons present in figure and save only those who have been edited
	fid = fopen(fname, 'w');
	if (fid < 0),	errordlg(['Can''t open file:  ' fname],'Error'),	return,		end
	h = findobj('Tag','GMT_DBpolyline');
	for (k = 1:numel(h))
		if ( isempty(getappdata(h(k), 'edited')) ),		continue,	end		% Skip because it was not modified
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
	else					doSave_formated(xx, yy, zz)
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
		else					h = h{1};		% Really use this handle.
		end
		xx = get(h,'XData');    yy = get(h,'YData');
		doSave_formated(xx, yy)
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
			fprintf(fid, ['%.2f\t' repmat('%f\t',[1,size(opt.z,2)]) '\n'], [t(:) double(opt.z)]');
			fclose(fid);
		end
	end

% -----------------------------------------------------------------------------------------
function doSave_formated(xx, yy, opt_z)
% Save x,y[,z] vars into a file but taking into account the 'LabelFormatType'
% OPT_Z is what the name says, optional
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
			xy = [xx(:) yy(:)];
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
	end
	
	if (nargin == 3)      
		xy = [xy opt_z(:)];    fmt = [fmt '\t%f'];
	end
	double2ascii(f_name,xy,fmt,'maybeMultis');

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
	if (ispc),		dos(['wget "' url '" -q --tries=2 --connect-timeout=5 -O ' dest_fiche]);
	else			unix(['wget ''' url ''' -q --tries=2 --connect-timeout=5 -O ' dest_fiche]);
	end

	fid = fopen(dest_fiche,'r');
	todos = fread(fid,'*char');
	if (numel(todos) < 100)
		warndlg('This station has no data or a file tranfer error occured.','Warning')
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
	h = findobj(hf,'Tag','add_uictx');
	cb = get(h, 'Call');
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
		msgbox(sprintf('%s',LineInfo),'This line info')
	else
		msgbox('Could not find reference to this object','This line info')
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

% 	s = getappdata(h,'polygon_data');
% 	if (~isempty(s))            % If the rectangle is in edit mode, force it out of edit
% 		if strcmpi(s.controls,'on'),    ui_edit_polygon(h);     end
% 	end
% 	set(h, 'HandleVis','off')           % Make the rectangle handle invisible
% 	hAxes = get(h,'Parent');
% 
% 	hLines = findobj(hAxes,'Type','line');     % Fish all objects of type line in Mirone figure
% 	hPatch = findobj(hAxes,'Type','patch');
% 	hText = findobj(hAxes,'Type','text');
% 	hLP = [hLines(:); hPatch(:)];
% 	rx = get(h,'XData');        ry = get(h,'YData');
% 	rx = [min(rx) max(rx)];     ry = [min(ry) max(ry)];
% 	found = false;
% 	for (i = 1:numel(hLP))		% Loop over objects to find if any is on edit mode
% 		s = getappdata(hLP(i),'polygon_data');
% 		if (~isempty(s))
% 			if strcmpi(s.controls,'on')     % Object is in edit mode, so this
% 				ui_edit_polygon(hLP(i))     % call will force out of edit mode
% 				found = true;
% 			end
% 		end
% 	end
% 	if (found)		% We have to do it again because some line handles have meanwhile desapeared
% 		hLines = findobj(hAxes,'Type','line');
% 		hPatch = findobj(hAxes,'Type','patch');
% 		hLP = [hLines(:); hPatch(:)];
% 	end
% 	set(h, 'HandleVis','on')    % Make the rectangle handle findable again
% 
% 	de_LP = false(numel(hLP),1);
% 	for (i = 1:numel(hLP))		% Loop over objects to find out which cross the rectangle
% 		x = get(hLP(i),'XData');        y = get(hLP(i),'YData');
% 		if ( any( (x >= rx(1) & x <= rx(2)) & (y >= ry(1) & y <= ry(2)) ) )
% 			de_LP(i) = true;
% 			delete(hLP(i))
% 		end
% 	end
% 
% 	found = false;
% 	for (i = 1:numel(hText))	% Text objs are a bit different, so treat them separately
% 		pos = get(hText(i),'Position');
% 		if ( (pos(1) >= rx(1) && pos(1) <= rx(2)) && (pos(2) >= ry(1) && pos(2) <= ry(2)) )
% 			delete(hText(i))
% 			found = true;
% 		end
% 	end
% 	if (found),     refresh;    end     % Bloody text bug

% -----------------------------------------------------------------------------------------
function trim_withPolygon(obj, evt, h, side)
% Intention is to ...
% H is the handle to polyg
% SIDE -> == 1, kill inside, else kill outside

	hL = fish_inside_rect(h, 'line');
	if (isempty(hL)),	return,		end

	hAxes = get(h,'Parent');		hFig = get(hAxes,'Parent');

	x0 = get(h,'XData');			y0 = get(h,'YData');	x0 = x0(:);		y0 = y0(:);

	for (k = 1:numel(hL))
		x = get(hL(k),'XData');		y = get(hL(k),'YData');
		IN = inpolygon(x,y, x0,y0);					% Find which ones are inside the polygon
		if (~side),		IN = ~IN;		end
		if (~any(IN)),		continue,	end			% It didn't cross the polygon
		ind = find(IN);
		if (numel(ind) == numel(x))					% Object entirely inside
			delete(hL(k))							% Not wanted, delete it
			continue
		end
		[xc, yc, iout] = intersections(x, y, x0, y0, 1);
		[iout, I] = sort(iout);			% Shit but they aren't always sorted
		if (any(diff(I) ~= 1))			% When shit happens we need to reorder the crossings too
			xc = xc(I);		yc = yc(I);
		end

% 		if (0 && numel(ind) == 1 && numel(xc) == 1)	% Easier case, we have a line with only one vertex inside.
% 			if (ind == 1)						% A line with only first vertex inside, which we want to kill
% 				x(1) = xc(1);	y(1) = yc(1);
% 			else								% A line with last vertex inside, which we want to kill
% 				x(end) = xc(1);	y(end) = yc(1);
% 			end
% 			set(hL(k),'XData',x, 'YData',y)
% 		else
% 			if (numel(xc) == 1)					% Better, only one crossing
% 				if (IN(1))						% Line starts inside the polygon and ends outside
% 					x(ind(end)) = xc(1);	y(ind(end)) = yc(1);
% 					x(1:end-1) = [];		y(1:end-1) = [];
% 				else							% The other way around. Line ends inside polygon
% 					x(ind(1)) = xc(1);		y(ind(1)) = yc(1);
% 					x(ind(1)+1:end) = [];	y(ind(1)+1:end) = [];
% 				end
% 				set(hL(k),'XData',x, 'YData',y)
% 			elseif (numel(xc) == 2)				% Not so bad yet.
% 				if (IN(1))						% Line starts inside polygon (easier, only need to shrink ends)
% 					n = fix(iout(2)) + 1;		% Indice of first point to remove (after re-entrance)
% 					x(n) = xc(2);			y(n) = yc(2);
% 					x(n+1:end) = [];		y(n+1:end) = [];
% 					n = fix(iout(1));			% Indice of last point to remove (before getting out)
% 					x(n) = xc(1);			y(n) = yc(1);
% 					x(1:n-1) = [];			y(1:n-1) = [];
% 				else							% Line starts outside. Now need to break it into two lines
% 					n = fix(iout(1)) + 1;		% Indice of last point to retain (after entrance)
% 					x(n) = xc(1);			y(n) = yc(1);
% 					x(n+1:end) = [];		y(n+1:end) = [];
% 					hCopy = copy_line_object([], hL(k), hFig, hAxes);
% 					xx = get(hCopy,'XData');	yy = get(hCopy,'YData');
% 					n = fix(iout(2));			% Indice of first point to retain before second exit
% 					xx(n) = xc(2);			yy(n) = yc(2);
% 					xx(1:n-1) = [];			yy(1:n-1) = [];
% 					set(hCopy,'XData',xx, 'YData',yy)
% 				end
% 				set(hL(k),'XData',x, 'YData',y)
% 			else								% WORST, several crossings in and out. We'll need to create new lines
				%ind2 = find(~diff(diff(ind) == 1)) + 1;	% Convoluted but it gives the indices of the to be removed pts
				difa = diff(ind);
				if (any(difa == 1))
					ind2 = find(difa == 1);
					ind2 = ind2(diff(ind2) == 1) + 1;	% Find inner sequence of consecutive numbers
					ind(ind2) = [];						% This guy now has the segment's boundaries
					if (IN(end) && ind(end) == ind(end-1)+1)
						ind(end) = [];					% Last pt is inside but not needed for reuse. Must kill it now
					end
				end
				if (IN(1) && ind(1) ~= fix(iout(1)))	% First point survived but it was not the nearest neighbor to crossing
					ind(1) = [];
				end

				% Search for vertices that are 'singletons' inside. Those need to be duplicated to account for the 2x crossing
				c = false(1, numel(ind));
				for (m = 1:numel(ind))
					if (ind(m) == 1),	continue,	end	% First vertex does not need duplication
					difa = ind(m) - fix(iout);
					ind2 = find((ind(m) - fix(iout)) == 0);
					if (~isempty(ind2) && ind2 > 1 && difa(ind2-1) == 1)
						c(m) = true;
					end
				end
				if (any(c))
					ind2 = ind(c);						% These pt guys need to be duplicated
					ind3 = [1:ind2(1) ind2(1)];
					for (m = 2:numel(ind2))				% Awful, we need to grow the array expensively
						ind3 = [ind3 ind2(m-1)+1:ind2(m) ind2(m)];
					end
					ind3 = [ind3 ind2(end)+1:numel(x)];
					x = x(ind3);	y = y(ind3);
					% Now also increase the 'ind'
					ind = sort([ind ind2]);
					ind2 = find(diff(ind) == 0) + 1;
					for (m = 1:numel(ind2))
						ind(ind2(m):end) = ind(ind2(m):end) + 1;	% Increase the second of the repetition by 1
					end
				end

				if (~IN(1) && ind(1) ~= 1)				% Line starts outside, need to add first index
					ind = [1 ind];
				end

				% Deal with last segment case
				n1 = numel(x);		n2 = numel(ind);	% To shut up the annoying MLint
				if (numel(xc) == numel(ind))
					x(ind) = xc;		y(ind) = yc;
					if (rem(numel(xc), 2) == 1)			% Odd number, line satrts inside but ends outside
						ind(n2+1) = n1;
					end
				elseif (numel(xc) == numel(ind) + 1)	% More convoluted case with only one last point inside
					x(ind) = xc(1:numel(ind));	y(ind) = yc(1:numel(ind));
					if (ind(end) == fix(iout(end)))		% We need to insert an extra pt
						xx = x(end);		yy = y(end);
						x(end) = xc(end);	y(end) = yc(end);
						x(n1+1) = xx;		y(n1+1) = yy;
						ind(n2+1) = ind(end) + 1;
						ind(n2+2) = numel(x);
					else
						x(fix(iout(end))) = xc(end);
						y(fix(iout(end))) = yc(end);
						ind(n2+1) = fix(iout(end));
						ind(n2+2) = numel(x);
					end
				elseif (numel(xc) == numel(ind) - 1)	% Also convoluted case. Line starts outside so we added 1 to ind
					x(ind(2:end)) = xc;		y(ind(2:end)) = yc;
					if (rem(numel(ind), 2))				% Hammering thinking, might very screw up
						ind(n2+1) = n1;
					end
				else
					error('Deu merda1')
				end
				break_trimeds(hL(k), ind, hFig, hAxes, x, y);
%  			end
% 		end

	end

	function break_trimeds(h, ind, hFig, hAxes, x, y)
	% ...
	set(h, 'XData',x(ind(1):ind(2)), 'YData',y(ind(1):ind(2)))
	if (numel(ind) == 2),	return,		end		% We rae done, only one segment
	hCopy = copy_line_object([], h, hFig, hAxes);
	set(hCopy, 'XData',[], 'YData',[])			% So that we will not copy potentially big lines for nothing
	for (k = 3:2:numel(ind)-2)					% From second to before last segments
		hNext = copy_line_object([], hCopy, hFig, hAxes);
		set(hNext, 'XData',x(ind(k):ind(k+1)), 'YData',y(ind(k):ind(k+1)))
	end

	hNext = copy_line_object([], hCopy, hFig, hAxes);
	set(hNext, 'XData',x(ind(end-1):ind(end)), 'YData',y(ind(end-1):ind(end)))	% Last segment
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
	else					sout = s(:,min(c):max(c));
	end

% --------------------------------------------------------------------
function set_transparency(obj,eventdata, h_patch)
% Sets the transparency of a patch object

	if (nargin == 2),	h_patch = gco;	end
	p_fc = get(h_patch,'FaceColor');
	if ( strcmpi(p_fc,'none') )
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
	uicontrol('style','slider','units','normalized','position',pos, 'Call',S,'min',0,'max',1,'Value',transp);

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

	uimenu(cmenuHand, 'Label', 'Delete', 'Call', 'delete(gco)');
	item_lw = uimenu(cmenuHand, 'Label', 'Line Width', 'Sep','on');
	setLineWidth(item_lw,cb_LineWidth)
	item_ls = uimenu(cmenuHand, 'Label', 'Line Style');
	setLineStyle(item_ls,{cb_solid cb_dashed cb_dotted cb_dashdot})
	item7 = uimenu(cmenuHand, 'Label', 'Line Color');
	cb_color = uictx_color(h,'EdgeColor');      % there are 9 cb_color outputs
	setLineColor(item7,cb_color)

	set_stack_order(cmenuHand)      % Change order in the stackpot

	uimenu(item7, 'Label', 'None', 'Sep','on', 'Call', 'set(gco, ''EdgeColor'', ''none'');refresh');
	item8 = uimenu(cmenuHand, 'Label','Fill Color', 'Sep','on');
	cb_color = uictx_color(h,'facecolor');      % there are 9 cb_color outputs
	setLineColor(item8,cb_color)
	uimenu(item8, 'Label', 'None', 'Sep','on', 'Call', 'set(gco, ''FaceColor'', ''none'');refresh');
	uimenu(cmenuHand, 'Label', 'Transparency', 'Call', @set_transparency);

% -----------------------------------------------------------------------------------------
function set_stack_order(cmenuHand)
% Change order in the stackpot. cmenuHand is what it says. 
	item_order = uimenu(cmenuHand, 'Label', 'Order');
	uimenu(item_order, 'Label', 'Bring to Top', 'Call','uistack_j(gco,''top'')');
	uimenu(item_order, 'Label', 'Send to Bottom', 'Call','uistack_j(gco,''bottom'')');
	uimenu(item_order, 'Label', 'Move up', 'Call','uistack_j(gco,''up'')');
	uimenu(item_order, 'Label', 'Move down', 'Call','uistack_j(gco,''down'')');

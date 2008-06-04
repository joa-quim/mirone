function OUT = draw_funs(hand,opt,data)
%   This contains several functions necessary to the "Draw" menu of mirone
%   There are no error checking.
%   HAND    contains the handle to the graphical object
%   OPT     is a string for choosing what action to perform
%   DATA    contains data currently used in the volcanoes, fogspots and some other option
%   OUT     Is currently used only as an option to 'ImportLine'. Data is returned rather than ploted

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

switch opt
    case 'line_uicontext',          set_line_uicontext(hand,'line')
    case 'SHPuictx',                set_SHPline_uicontext(hand)
    case 'ContourLines',            set_ContourLines_uicontext(hand,data)
    case 'MBtrackUictx',            set_line_uicontext(hand,'MBtrack')
    case 'MBbarUictx',              set_bar_uicontext(hand)
    case 'CoastLineUictx',          setCoastLineUictx(hand)
    case 'DeleteObj',               delete_obj(hand);
    case 'DrawGreatCircle'
        h = draw_greateCircle;
        if ~isempty(h)      % when in compiled version h may be empty (why?).
            set_greatCircle_uicontext(h)
        end
    case 'DrawCircleEulerPole'
        h = draw_circleEulerPole(data(1),data(2));
        if ~isempty(h)      % when in compiled version h may be empty (why?).
            set_circleGeo_uicontext(h)
        end
    case 'DrawCartesianCircle'
        h = draw_circleGeo;         % It also draws cartesian circles
        if ~isempty(h)              % when in compiled version h may be empty (why?).
            set_circleCart_uicontext(h)
        end
    case 'SessionRestoreCircle'     % Called by "FileOpenSession" or "DrawGeographicalCircle_CB"
        set_circleGeo_uicontext(hand)
    case 'SessionRestoreCircleCart'     % Called by "FileOpenSession" or "DrawGeographicalCircle_CB"
        set_circleCart_uicontext(hand)
    case 'DrawText'
        cmenuHand = uicontextmenu;
        set(hand, 'UIContextMenu', cmenuHand);
        cb1 = 'set(gco, ''Editing'', ''on''); refresh';
        cb_color = uictx_color(hand);      % there are 9 cb_color outputs
        uimenu(cmenuHand, 'Label', 'Change Font', 'Callback', @text_FontSize);
        item_fc = uimenu(cmenuHand, 'Label', 'Font Color');
        setLineColor(item_fc,cb_color)
        uimenu(cmenuHand, 'Label', 'Edit   text', 'Callback', cb1, 'Separator','on');
        uimenu(cmenuHand, 'Label', 'Copy   text', 'Callback', @copy_text_object);
        uimenu(cmenuHand, 'Label', 'Delete text', 'Callback', 'delete(gco); refresh');
        uimenu(cmenuHand, 'Label', 'Move   text', 'Callback', @move_text);
        uimenu(cmenuHand, 'Label', 'Rotate text', 'Callback', @rotate_text);
        uimenu(cmenuHand, 'Label', 'Export text', 'Callback', @export_text);
    case 'DrawSymbol'
        set_symbol_uicontext(hand)
    case 'ImportLine'                   % read AND plot the line
        fname = hand;
        hFig = get(0,'CurrentFigure');         hAxes = get(hFig,'CurrentAxes');
        [bin,n_column,multi_seg,n_headers] = guess_file(fname);
        % If msgbox exist we have to move it from behind the main window. So get it's handle
        hMsgFig = get(0,'CurrentFigure');
        if (hFig ~= hMsgFig),       figure(hMsgFig);   end   % If msgbox exists, bring it forward
        % If error in reading file
        if isempty(bin) && isempty(n_column) && isempty(multi_seg) && isempty(n_headers)
            errordlg(['Error reading file ' fname],'Error');    return
        end
        if (bin ~= 0)   % NOT ASCII
            errordlg('Sorry, reading binary files is not yet programed','Error');   return
        end
        if (n_column < 2)
            errordlg('File error. Your file doesn''t have at least 2 columns','Error'); return
        end
        if (isempty(n_headers)),    n_headers = NaN;    end
        if (multi_seg)
            [numeric_data,multi_segs_str,headerlines] = text_read(fname,NaN,n_headers,'>');
        else
            [numeric_data,multi_segs_str,headerlines] = text_read(fname,NaN,n_headers);
        end
        
        % Project if we need
        handles = guidata(hFig);
        if (handles.is_projected && handles.defCoordsIn > 0)
            try
                if (iscell(numeric_data))
                    for i=1:numel(numeric_data)
                        numeric_data{i}  = geog2projected_pts(handles,numeric_data{i});
                    end
                else
                    numeric_data = geog2projected_pts(handles,numeric_data);
                end
            catch
                errordlg(lasterr,'ERROR');    return
            end
        end
        
        % If OUT is requested there is nothing left to be done here
        if (nargout)
            OUT = numeric_data;     return
        end
            
        if (hFig ~= hMsgFig);       figure(hFig);    axes(hAxes);   end     % gain access to the drawing figure
        % Get rid of points that are outside the map limits
        tol = 0.5;
        if (iscell(numeric_data))
            n_segments = length(numeric_data);
        else
            n_segments = 1;
        end
        XYlim = getappdata(handles.axes1,'ThisImageLims');
        xx = XYlim(1:2);            yy = XYlim(3:4);
        hold on
        lt = handles.DefLineThick;   lc = handles.DefLineColor;
        for i=1:n_segments
            if (iscell(numeric_data))
                tmpx = numeric_data{i}(:,1);    tmpy = numeric_data{i}(:,2);
            else
                tmpx = numeric_data(:,1);       tmpy = numeric_data(:,2);
            end
            ind = find(tmpx < xx(1)-tol | tmpx > xx(2)+tol);
            tmpx(ind) = [];         tmpy(ind) = [];
            ind = find(tmpy < yy(1)-tol | tmpy > yy(2)+tol);
            tmpx(ind) = [];         tmpy(ind) = [];
            switch data
                case 'AsLine'
                    % The following Tag is very important to tell from MB tracks, which have Tags = MBtrack#
                    lineHand = plot(tmpx,tmpy,'Color',lc,'LineWidth',lt,'Tag','polyline');
                    set_line_uicontext(lineHand,'line')     % Set lines's uicontextmenu
                case 'AsPoint'
                    lineHand = plot(tmpx,tmpy,'ko','MarkerEdgeColor','w','MarkerFaceColor','k', ...
                        'MarkerSize',4,'Tag','Pointpolyline');
                    set_symbol_uicontext(lineHand)          % Set marker's uicontextmenu (tag is very important)
                case 'AsMaregraph'
                    lineHand = plot(tmpx,tmpy,'Marker','o','MarkerFaceColor','y',...
                        'MarkerEdgeColor','k','MarkerSize',10,'Tag','Maregraph');
                    set_symbol_uicontext(lineHand)          % Set marker's uicontextmenu
                case 'FaultTrace'
                    lineHand = plot(tmpx,tmpy,'Color',lc,'LineWidth',lt,'Tag','FaultTrace');
                    set_line_uicontext(lineHand,'line')     % Set lines's uicontextmenu
                    % Create empty patches that will contain the surface projection of the fault plane
                    for (k=1:length(tmpx)-1),   hp(k) = patch('XData', [], 'YData',[]);    end
                    setappdata(lineHand,'PatchHand',hp);
            end
        end
        clear numeric_data;     hold off
    case {'hotspot','volcano','ODP','City_major','City_other','Earthquakes','TideStation'}
        set_symbol_uicontext(hand,data)
    case 'PlateBoundPB',		set_PB_uicontext(hand,data)
    case 'DrawVector'
        h = draw_vector;
        if ~isempty(h)      % when in compiled version h may be empty.
            set_vector_uicontext(h)
        end
    case 'ChngAxLabels',		changeAxesLabels(data)
    case 'MagBarCode',			draw_MagBarCode
    case 'SRTMrect',			set_SRTM_rect_uicontext(hand)
    case 'isochron',			set_isochrons_uicontext(hand,data)
    case 'gmtfile',				set_gmtfile_uicontext(hand,data)
    case 'country_patch',		set_country_uicontext(hand)
    case 'telhas_patch',		set_telhas_uicontext(hand)
    case 'save_xyz',			save_formated([],[],[], data)
    case 'tellAzim',			show_lineAzims([],[], hand);
    case 'tellLLength',			show_LineLength([],[], hand);
    case 'tellArea',			show_Area([],[], hand);
end

% -----------------------------------------------------------------------------------------
function setLineStyle(item,cbs)
	% Set the line Style uicontexts of graphic elements
	uimenu(item, 'Label', 'solid', 'Callback', cbs{1});
	uimenu(item, 'Label', 'dashed', 'Callback', cbs{2});
	uimenu(item, 'Label', 'dotted', 'Callback', cbs{3});
	uimenu(item, 'Label', 'dash-dotted', 'Callback', cbs{4});

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
function set_SHPline_uicontext(h,opt)
% h is a handle to a shape line object

handles = guidata(h(1));
for (i = 1:numel(h))
	cmenuHand = uicontextmenu('Parent',handles.figure1);      set(h(i), 'UIContextMenu', cmenuHand);
    uimenu(cmenuHand, 'Label', 'Save line', 'Callback', {@save_formated,h});
	uimenu(cmenuHand, 'Label', 'Delete this line', 'Callback', {@del_line,h(i)});
	uimenu(cmenuHand, 'Label', 'Delete class', 'Callback', 'delete(findobj(''Tag'',''SHPpolyline''))');
	%ui_edit_polygon(h(i))    % Set edition functions   
	
	cb_solid = 'set(gco, ''LineStyle'', ''-''); refresh';
	cb_dashed = 'set(gco, ''LineStyle'', ''--''); refresh';
	cb_dotted = 'set(gco, ''LineStyle'', '':''); refresh';
	cb_dashdot = 'set(gco, ''LineStyle'', ''-.''); refresh';
	
	item_lw = uimenu(cmenuHand, 'Label', 'Line Width', 'Separator','on');
	uimenu(item_lw, 'Label', 'Other...', 'Callback', {@other_LineWidth,h(i)});
	
	item_ls = uimenu(cmenuHand, 'Label', 'Line Style');
	setLineStyle(item_ls,{cb_solid cb_dashed cb_dotted cb_dashdot})
	
	item7 = uimenu(cmenuHand, 'Label', 'Line Color');
	uimenu(item7, 'Label', 'Other...', 'Callback', {@other_color,h(i)});    
end

% -----------------------------------------------------------------------------------------
function set_line_uicontext(h,opt)
% h is a handle to a line object (that can be closed)

IS_SEISPOLYGON = 0;     % Seismicity polygons have special options
% Check to see if we are dealing with a closed polyline
x = get(h,'XData');   y = get(h,'YData');
if (isempty(x) || isempty(y)),   return;     end     % Line is totally out of the figure
if ( (x(1) == x(end)) && (y(1) == y(end)) )
    LINE_ISCLOSED = 1;
    IS_RECTANGLE = 0;
    if ( length(x) == 5 && (x(1) == x(2)) && (x(3) == x(4)) && (y(1) == y(4)) && (y(2) == y(3)) )
        IS_RECTANGLE = 1;
    end
    if (strcmp(get(h,'Tag'),'SeismicityPolygon')),  IS_SEISPOLYGON = 1;    end
else
    LINE_ISCLOSED = 0;
    IS_RECTANGLE = 0;       % If line is not closed, it cannot be a rectangle
end
if (strcmp(get(h,'Type'),'patch')), IS_PATCH = 1;
else                                IS_PATCH = 0;
end

handles = guidata(get(h,'Parent'));             % Get Mirone handles

% Check to see if we are dealing with a multibeam track
cmenuHand = uicontextmenu('Parent',handles.figure1);
set(h, 'UIContextMenu', cmenuHand);
switch opt
    case 'line'
        label_save = 'Save line';   label_length = 'Line length(s)';   label_azim = 'Line azimuth(s)';
        IS_LINE = 1;    IS_MBTRACK = 0;
    case 'MBtrack'
        label_save = 'Save track';   label_length = 'Track length';   label_azim = 'Track azimuth(s)';
        IS_LINE = 0;    IS_MBTRACK = 1;
end
cb_LineWidth = uictx_LineWidth(h);      % there are 5 cb_LineWidth outputs
cb_solid = 'set(gco, ''LineStyle'', ''-''); refresh';
cb_dashed = 'set(gco, ''LineStyle'', ''--''); refresh';
cb_dotted = 'set(gco, ''LineStyle'', '':''); refresh';
cb_dashdot = 'set(gco, ''LineStyle'', ''-.''); refresh';
cb_color = uictx_color(h);      % there are 9 cb_color outputs

if (IS_RECTANGLE)
	uimenu(cmenuHand, 'Label', 'Delete me', 'Callback', {@del_line,h});
	uimenu(cmenuHand, 'Label', 'Delete inside rect', 'Callback', {@del_insideRect,h});
	ui_edit_polygon(h)
elseif (IS_LINE)
	uimenu(cmenuHand, 'Label', 'Delete', 'Callback', {@del_line,h});
	ui_edit_polygon(h)		% Set edition functions
elseif (IS_MBTRACK)			% Multibeam tracks, when deleted, have to delete also the bars
	uimenu(cmenuHand, 'Label', 'Delete track (left-click on it)', 'Callback', 'save_track_mb(1);');
	% Old style edit function. New edit is provided by ui_edit_polygon which doesn't work with mbtracks 
	uimenu(cmenuHand, 'Label', 'Edit track (left-click on it)', 'Callback', 'edit_track_mb');
end
uimenu(cmenuHand, 'Label', label_save, 'Callback', {@save_formated,h});
if (~IS_SEISPOLYGON && ~IS_MBTRACK && ~strcmp(get(h,'Tag'),'FaultTrace'))     % Those are not to allowed to copy
	uimenu(cmenuHand, 'Label', 'Copy', 'Callback', {@copy_line_object,handles.figure1,handles.axes1});
end
if (~IS_SEISPOLYGON),	uimenu(cmenuHand, 'Label', label_length, 'Callback', @show_LineLength);		end
if (IS_MBTRACK),		uimenu(cmenuHand, 'Label', 'All tracks length', 'Callback', @show_AllTrackLength);	end
if (~IS_SEISPOLYGON),	uimenu(cmenuHand, 'Label', label_azim, 'Callback', @show_lineAzims);	end

if (LINE_ISCLOSED)
    uimenu(cmenuHand, 'Label', 'Area under polygon', 'Callback', @show_Area);
	if (~IS_RECTANGLE && ~handles.validGrid)
		uimenu(cmenuHand, 'Label', 'Crop Image', 'Callback', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco)','Sep','on');
		if (handles.image_type == 3)
				uimenu(cmenuHand, 'Label', 'Crop Image (with coords)', 'Callback', ...
				'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''CropaWithCoords'')');
		end
	end
	if (IS_PATCH && ~IS_SEISPOLYGON)
		item8 = uimenu(cmenuHand, 'Label','Fill Color');
		setLineColor( item8, uictx_color(h, 'facecolor') )		% there are 9 cb_color outputs
		uimenu(item8, 'Label', 'None', 'Separator','on', 'Call', 'set(gco, ''FaceColor'', ''none'');refresh');
		uimenu(cmenuHand, 'Label', 'Transparency', 'Call', @set_transparency);
	end
	uimenu(cmenuHand, 'Label', 'Create Mask', 'Call', 'poly2mask_fig(guidata(gcbo),gco)');
end

if ( ~LINE_ISCLOSED && strcmp(opt,'line') && (ndims(get(handles.hImg,'CData')) == 2 || handles.validGrid) )
	cbTrack = 'setappdata(gcf,''TrackThisLine'',gco); mirone(''ExtractProfile_CB'',guidata(gcbo),''point'')';
	uimenu(cmenuHand, 'Label', 'Point interpolation', 'Call', cbTrack);
	cbTrack = 'setappdata(gcf,''TrackThisLine'',gco); mirone(''ExtractProfile_CB'',guidata(gcbo))';
	uimenu(cmenuHand, 'Label', 'Extract profile', 'Call', cbTrack);
end

if strcmp(opt,'MBtrack'),	uimenu(cmenuHand, 'Label', 'Show track''s Swath Ratio', 'Call', {@show_swhatRatio,h});	end

if (IS_RECTANGLE)
	uimenu(cmenuHand, 'Label', 'Rectangle limits', 'Separator','on', 'Call', @rectangle_limits);
	uimenu(cmenuHand, 'Label', 'Crop Image', 'Call', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco)');
	if (handles.image_type == 3 || handles.validGrid)
		uimenu(cmenuHand, 'Label', 'Crop Image (with coords)', 'Call', ...
			'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''CropaWithCoords'')');
	end
	uimenu(cmenuHand, 'Label', 'Register Image', 'Call', @rectangle_register_img);
	uimenu(cmenuHand, 'Label', 'Transplant Image here', 'Call', @Transplant_Image);
	if (handles.validGrid)    % Option only available to recognized grids
		item_tools = uimenu(cmenuHand, 'Label', 'Crop Tools','Separator','on');
		uimenu(item_tools, 'Label', 'Spline smooth', 'Call', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''SplineSmooth'')');
		uimenu(item_tools, 'Label', 'Median filter', 'Call', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''MedianFilter'')');
		uimenu(item_tools, 'Label', 'Crop Grid', 'Call', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''CropaGrid_pure'')');
		uimenu(item_tools, 'Label', 'Histogram', 'Call', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''CropaGrid_histo'')');
		uimenu(item_tools, 'Label', 'Power', 'Call', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''CropaGrid_power'')');
		uimenu(item_tools, 'Label', 'Autocorrelation', 'Call', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''CropaGrid_autocorr'')');
		uimenu(item_tools, 'Label', 'FFT tool', 'Call', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''CropaGrid_fftTools'')');
		item_fill = uimenu(item_tools, 'Label', 'Fill gaps');
		uimenu(item_fill, 'Label', 'Fill gaps (surface)', 'Call', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''FillGaps'',''surface'')');
		uimenu(item_fill, 'Label', 'Fill gaps (cubic)', 'Call', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''FillGaps'',''cubic'');');
		uimenu(item_fill, 'Label', 'Fill gaps (linear)', 'Call', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''FillGaps'',''linear'');');
		uimenu(item_tools, 'Label','Set to constant', 'Call', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''SetConst'')');
	end
end
item_lw = uimenu(cmenuHand, 'Label', 'Line Width', 'Separator','on');
setLineWidth(item_lw,cb_LineWidth)
item_ls = uimenu(cmenuHand, 'Label', 'Line Style');
setLineStyle(item_ls,{cb_solid cb_dashed cb_dotted cb_dashdot})
item7 = uimenu(cmenuHand, 'Label', 'Line Color');
if (IS_PATCH),		cb_color = uictx_color(h,'EdgeColor');	end      % there are 9 cb_color outputs
setLineColor(item7,cb_color)

set_stack_order(cmenuHand)      % Change order in the stackpot

if (LINE_ISCLOSED && ~IS_SEISPOLYGON)
 	if (handles.validGrid && ~IS_RECTANGLE)    % Option only available to recognized grids
		item_tools2 = uimenu(cmenuHand, 'Label', 'ROI Crop Tools','Separator','on');
		uimenu(item_tools2, 'Label', 'Crop Grid', 'Call', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''CropaGrid_pure'')');
		uimenu(item_tools2, 'Label', 'Set to const', 'Call', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''ROI_SetConst'')');
		uimenu(item_tools2, 'Label', 'Histogram', 'Call', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''CropaGrid_histo'')');
		uimenu(item_tools2, 'Label', 'Median filter', 'Call', 'mirone(''ImageCrop_CB'',guidata(gcbo),gco,''ROI_MedianFilter'')');
    end
    if (strcmp(get(h,'Tag'),'EulerTrapezium'))
        uimenu(cmenuHand, 'Label', 'Compute Euler Pole', 'Separator','on', 'Callback',...
            'calcBoninEulerPole(get(gco,''XData''), get(gco,''YData''));' );
    end
    cb_roi = 'mirone(''DrawClosedPolygon_CB'',guidata(gcbo),gco)';
    uimenu(cmenuHand, 'Label', 'Region-Of-Interest', 'Separator','on', 'Callback', cb_roi);
end

if (strcmp(get(h,'Tag'),'FaultTrace'))      % For Okada modeling
	uimenu(cmenuHand, 'Label', 'Okada', 'Separator','on', 'Callback', {@okada_model,h,'okada'});    
	uimenu(cmenuHand, 'Label', 'Mansinha', 'Call', {@okada_model,h,'mansinha'});    
end

if (IS_SEISPOLYGON)                         % Seismicity options
	% gco gives the same handle as h 
	uimenu(cmenuHand, 'Label', 'Save events', 'Call', 'save_seismicity(gcf,[],gco)', 'Separator','on');
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
end

% -----------------------------------------------------------------------------------------
function copy_line_object(obj,eventdata,hFig,hAxes)
    oldH = gco;
	newH = copyobj(oldH,hAxes);
    h = findobj(get(newH,'uicontextmenu'),'label','Save line');
    if (~isempty(h))        % Replace the old line handle in the 'Save line' Callback by the just created one
        hFun = get(h,'Call');
        hFun{2} = newH;
        set(h,'Call',hFun)
    end
	rmappdata(newH,'polygon_data')          % Remove the parent's ui_edit_polygon appdata
	state = uisuspend_fig(hFig);            % Remember initial figure state
	x_lim = get(hAxes,'xlim');        y_lim = get(hAxes,'ylim');
	current_pt = get(hAxes, 'CurrentPoint');
	setappdata(newH,'old_pt',[current_pt(1,1) current_pt(1,2)])
	
	set(hFig,'WindowButtonMotionFcn',{@wbm_MovePolygon,newH,[x_lim y_lim],hAxes},...
        'WindowButtonDownFcn',{@wbd_MovePolygon,newH,state}, 'Pointer','fleur');

% ---------
function wbm_MovePolygon(obj,eventdata,h,lim,hAxes)
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
	uirestore_fig(state);           % Restore the figure's initial state
	ui_edit_polygon(h)              % Reset the edition functions with the correct handle
% -----------------------------------------------------------------------------------------

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
	cb_LineWidth = uictx_LineWidth(h(i));      % there are 5 cb_LineWidth outputs
	item_lw = uimenu(cmenuHand, 'Label', 'Line Width', 'Separator','on');
	uimenu(item_lw, 'Label', '1     pt', 'Callback', cb_LineWidth{1});
	uimenu(item_lw, 'Label', 'Other...', 'Callback', cb_LineWidth{5});
	item8 = uimenu(cmenuHand, 'Label','Fill Color', 'Separator','on');
	cb_color = uictx_color(h(i),'facecolor');      % there are 9 cb_color outputs
	uimenu(item8, 'Label', 'Other...', 'Callback', cb_color{9});
	uimenu(item8, 'Label', 'None', 'Callback', 'set(gco, ''FaceColor'', ''none'');refresh');
	uimenu(cmenuHand, 'Label', 'Transparency', 'Callback', @set_transparency);
	uimenu(cmenuHand, 'Label', 'Create Mask', 'Call', 'poly2mask_fig(guidata(gcbo),gco)');
	if (handles.image_type ~= 20)
		uimenu(cmenuHand, 'Label', 'Region-Of-Interest', 'Separator','on', 'Callback', ...
			'mirone(''DrawClosedPolygon_CB'',guidata(gcbo),gco)');
	end
end

% -----------------------------------------------------------------------------------------
function okada_model(obj,eventdata,h,opt)
if (nargin == 3),   opt = 'okada';     end
hh = findobj('Tag','FaultTrace');       % Check if we have more than one (multi-segment?)faults
if (isempty(hh)),   errordlg('This is just a line, NOT a fault trace. Can''t you see the difference?','Error'); return; end
h_fig = get(0,'CurrentFigure');     handles = guidata(h_fig);
% Guess minimum length segment that could be due to a bad line drawing
if (handles.geog)
    min_len = 0.05;
else
    imgLims = getappdata(handles.axes1,'ThisImageLims');
    if (abs(imgLims(2) - imgLims(1)) < 5000),   min_len = 5;    % Assume that the grid is in km
    else                                        min_len = 5000; % Assume meters
    end
end
if (length(hh) > 1)
    for k=1:length(hh)
        xx = get(hh(k),'XData');    yy = get(hh(k),'YData');
        dx = diff(xx);  dy = diff(yy);              dr = sqrt(dx.*dx + dy.*dy);
        ind = find(dr < min_len);
        if (~isempty(ind) && length(xx) > 2)     % Remove too short segments
            xx(ind) = [];   yy(ind) = [];
            set(hh(k),'XData',xx,'YData',yy)
        end
        azim = show_lineAzims([],[],hh(k));
        az{k} = azim.az;
    end
    h = hh;
else
    xx = get(h,'XData');    yy = get(h,'YData');
    dx = diff(xx);      dy = diff(yy);          dr = sqrt(dx.*dx + dy.*dy);
    ind = find(dr < min_len);
    if (~isempty(ind) && length(xx) > 2)         % Remove too short segments
        xx(ind) = [];   yy(ind) = [];
        set(h,'XData',xx,'YData',yy)
    end
    azim = show_lineAzims([],[],h);
    az = azim.az;
end

if (strcmp(opt,'okada')),           deform_okada(handles,h,az);
elseif (strcmp(opt,'mansinha')),    deform_mansinha(handles,h,az);
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
	% cb1     = 'mirone(''DrawEditLine_Callback'',gcbo,[],guidata(gcbo))';
	ui_edit_polygon(h)            % Set edition functions
	cb_rac = {@remove_symbolClass,h};   % It will also remove the labels because they have the same tag.
	cb_LineWidth = uictx_LineWidth(h);      % there are 5 cb_LineWidth outputs
	cb18 = 'set(gco, ''LineStyle'', ''-''); refresh';   cb19 = 'set(gco, ''LineStyle'', ''--''); refresh';
	cb20 = 'set(gco, ''LineStyle'', '':''); refresh';   cb21 = 'set(gco, ''LineStyle'', ''-.''); refresh';
	cb_color = uictx_color(h);      % there are 9 cb_color outputs

	uimenu(cmenuHand, 'Label', 'Delete contour', 'Callback',{@remove_singleContour,h});
	uimenu(cmenuHand, 'Label', 'Delete all contours', 'Callback', cb_rac);
	% item1 = uimenu(cmenuHand, 'Label', 'Edit contour (left-click on it)', 'Callback', cb1);
	uimenu(cmenuHand, 'Label', 'Save contour', 'Callback', {@save_formated,h});
	uimenu(cmenuHand, 'Label', 'Contour length', 'Callback', {@show_LineLength,[]});
	uimenu(cmenuHand, 'Label', 'Area under contour', 'Callback', @show_Area);
	item_lw = uimenu(cmenuHand, 'Label', 'Contour Line Width', 'Separator','on');
	setLineWidth(item_lw,cb_LineWidth)
	item_ls = uimenu(cmenuHand, 'Label', 'Contour Line Style');
	setLineStyle(item_ls,{cb18 cb19 cb20 cb21})
	item_lc = uimenu(cmenuHand, 'Label', 'Contour Line Color');
	setLineColor(item_lc,cb_color)
	cb_CLineWidth = uictx_Class_LineWidth(h);           % there are 5 cb_CLineWidth outputs
	item8 = uimenu(cmenuHand, 'Label', 'All Contours Line Width', 'Separator','on');
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
	% h is a handle to a line object
	tag = get(h,'Tag');
	if (strcmp(tag,'CoastLineNetCDF')),         label = 'Delete coastlines';
	elseif (strcmp(tag,'PoliticalBoundaries')), label = 'Delete boundaries';
	elseif (strcmp(tag,'Rivers')),              label = 'Delete rivers';
	end
	handles = guidata(h);
	cmenuHand = uicontextmenu('Parent',handles.figure1);
	set(h, 'UIContextMenu', cmenuHand);
	cb_LineWidth = uictx_LineWidth(h);      % there are 5 cb_LineWidth outputs
	cb13 = 'set(gco, ''LineStyle'', ''-''); refresh';   cb14 = 'set(gco, ''LineStyle'', ''--''); refresh';
	cb15 = 'set(gco, ''LineStyle'', '':''); refresh';   cb16 = 'set(gco, ''LineStyle'', ''-.''); refresh';
	cb_color = uictx_color(h);              % there are 9 cb_color outputs
	
	uimenu(cmenuHand, 'Label', label, 'Call', 'delete(gco)');
	uimenu(cmenuHand, 'Label', 'Edit line (left-click on it)', 'Call', 'edit_line');
	uimenu(cmenuHand, 'Label', 'Save coastline', 'Call', {@save_formated,h});
	
	item_lw = uimenu(cmenuHand, 'Label', 'Line Width', 'Sep','on');
	setLineWidth(item_lw,cb_LineWidth)
	item_ls = uimenu(cmenuHand, 'Label', 'Line Style');
	setLineStyle(item_ls,{cb13 cb14 cb15 cb16})
	item_lc = uimenu(cmenuHand, 'Label', 'Line Color');
	setLineColor(item_lc,cb_color)

% -----------------------------------------------------------------------------------------
function set_PB_uicontext(h,data)
% h is a handle to the lines of the PB_All (P. Bird Plate Boundaries) object

for i = 1:7     % Loop over all Plate Boundaries Types
	h_cur = [];
	switch i
		case 1,            h_cur = h.OSR;  data_cur = data.OSR;    % class = 'OSR'
		case 2,            h_cur = h.OTF;  data_cur = data.OTF;    % class = 'OTF'
		case 3,            h_cur = h.CRB;  data_cur = data.CRB;    % class = 'CRB'
		case 4,            h_cur = h.CTF;  data_cur = data.CTF;    % class = 'CTF'
		case 5,            h_cur = h.CCB;  data_cur = data.CCB;    % class = 'CCB'
		case 6,            h_cur = h.OCB;  data_cur = data.OCB;    % class = 'OCB'
		case 7,            h_cur = h.SUB;  data_cur = data.SUB;    % class = 'SUB'
	end
	if (isempty(h_cur)),	continue,	end
	cmenuHand = uicontextmenu;
	set(h_cur, 'UIContextMenu', cmenuHand);
	cb_LineWidth = uictx_Class_LineWidth(h_cur);    % there are 5 cb_PB_LineWidth outputs
	cb_color = uictx_Class_LineColor(h_cur);        % there are 9 cb_PB_color outputs
	uimenu(cmenuHand, 'Label', 'Segment info', 'Callback', {@PB_All_Info,h_cur,data_cur});
	uimenu(cmenuHand, 'Label', 'Delete class', 'Callback', 'delete(findobj(''Tag'',''PB_All''))', 'Separator','on');
	uimenu(cmenuHand, 'Label', 'Segment length', 'Callback', {@show_LineLength,[]});
	item3 = uimenu(cmenuHand, 'Label', 'Line Width', 'Separator','on');
	uimenu(item3, 'Label', '2       pt', 'Callback', cb_LineWidth{2});
	uimenu(item3, 'Label', '3       pt', 'Callback', cb_LineWidth{3});
	uimenu(item3, 'Label', '4       pt', 'Callback', cb_LineWidth{4});
	uimenu(item3, 'Label', 'Other...', 'Callback', cb_LineWidth{5});
	item_lc = uimenu(cmenuHand, 'Label', 'Color');
	setLineColor(item_lc,cb_color)
end

% -----------------------------------------------------------------------------------------
function set_isochrons_uicontext(h,data)
% h are handles to the lines of isochrons (or other lines with a info)
	tag = get(h,'Tag');
	if (iscell(tag)),   tag = tag{1};   end

	handles = guidata(get(h(1),'Parent'));             % Get Mirone handles
	cmenuHand = uicontextmenu('Parent',handles.figure1);
	set(h, 'UIContextMenu', cmenuHand);
	cb_LineWidth = uictx_LineWidth(h);       % there are 5 cb_LineWidth outputs
	cb_color = uictx_color(h);               % there are 9 cb_color outputs
	cbls1 = 'set(gco, ''LineStyle'', ''-''); refresh';   cbls2 = 'set(gco, ''LineStyle'', ''--''); refresh';
	cbls3 = 'set(gco, ''LineStyle'', '':''); refresh';   cbls4 = 'set(gco, ''LineStyle'', ''-.''); refresh';
	if (~all(isempty(cat(2,data{:}))))
		uimenu(cmenuHand, 'Label', [tag ' info'], 'Callback', {@Isochrons_Info,data});
		uimenu(cmenuHand, 'Label', ['Delete this ' tag ' line'], 'Callback', {@del_line,h}, 'Separator','on');
	else
		uimenu(cmenuHand, 'Label', ['Delete this ' tag ' line'], 'Callback', {@del_line,h});
	end
	uimenu(cmenuHand, 'Label', ['Delete all ' tag ' lines'], 'Callback', {@remove_symbolClass,h});
	uimenu(cmenuHand, 'Label', ['Save this ' tag ' line'], 'Callback', @save_line);
	uimenu(cmenuHand, 'Label', ['Save all ' tag ' lines'], 'Callback', {@save_line,h});
	uimenu(cmenuHand, 'Label', 'Line azimuths', 'Callback', @show_lineAzims);
	uimenu(cmenuHand, 'Label', 'Line length', 'Callback', {@show_LineLength,[],'nikles'});
	LINE_ISCLOSED = 0;
	for i=1:length(h)
		x = get(h(i),'XData');      y = get(h(i),'YData');
		if ( (x(1) == x(end)) && (y(1) == y(end)) && length(x) > 1)      % See if we have at least one closed line
			LINE_ISCLOSED = 1;
		end
	end
	if (LINE_ISCLOSED)      % If at least one is closed, activate the Area option
		uimenu(cmenuHand, 'Label', 'Area under polygon', 'Callback', @show_Area);
	end
	item_lw = uimenu(cmenuHand, 'Label', 'Line Width', 'Separator','on');
	setLineWidth(item_lw,cb_LineWidth)
	item_ls = uimenu(cmenuHand, 'Label', 'Line Style');
	setLineStyle(item_ls,{cbls1 cbls2 cbls3 cbls4})
	item_lc = uimenu(cmenuHand, 'Label', 'Color');
	setLineColor(item_lc,cb_color)
	% --------- Now set the class properties
	cb_ClassColor = uictx_Class_LineColor(h);        % there are 9 cb_color outputs
	item_Class_lc = uimenu(cmenuHand, 'Label', ['All ' tag ' Color'], 'Separator','on');
	setLineColor(item_Class_lc,cb_ClassColor)
	cb_ClassLineWidth = uictx_Class_LineWidth(h);    % there are 5 cb_ClassLineWidth outputs
	item_Class_lw = uimenu(cmenuHand, 'Label', ['All ' tag ' Line Width']);
	uimenu(item_Class_lw, 'Label', '1       pt', 'Callback', cb_ClassLineWidth{1});
	uimenu(item_Class_lw, 'Label', '2       pt', 'Callback', cb_ClassLineWidth{2});
	uimenu(item_Class_lw, 'Label', '3       pt', 'Callback', cb_ClassLineWidth{3});
	uimenu(item_Class_lw, 'Label', '4       pt', 'Callback', cb_ClassLineWidth{4});
	uimenu(item_Class_lw, 'Label', 'Other...', 'Callback', cb_ClassLineWidth{5});
	cb_ClassLineStyle = uictx_Class_LineStyle(h);    % there are 4 cb_ClassLineStyle outputs
	item_Class_lt = uimenu(cmenuHand, 'Label', ['All ' tag ' Line Style']);
	setLineStyle(item_Class_lt,{cb_ClassLineStyle{1} cb_ClassLineStyle{2} cb_ClassLineStyle{3} cb_ClassLineStyle{4}})
	uimenu(cmenuHand, 'Label', 'Euler rotation', 'Separator','on', 'Callback', 'euler_stuff(gcf,gco)');
	for i=1:length(h),   ui_edit_polygon(h(i));     end		% Set edition functions

% -----------------------------------------------------------------------------------------
function set_gmtfile_uicontext(h,data)
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
	uimenu(cmenuHand, 'Label', ['Delete this ' tag ' line'], 'Callback', 'delete(gco)', 'Separator','on');
	uimenu(cmenuHand, 'Label', ['Save this ' tag ' line'], 'Callback', @save_line);
	%uimenu(cmenuHand, 'Label', 'Open with gmtedit', 'Callback', ['gmtedit(getappdata(gco,''FullName''))']);
	uimenu(cmenuHand, 'Label', 'Open with gmtedit', 'Callback', {@call_gmtedit,h});
	item_lw = uimenu(cmenuHand, 'Label', 'Line Width', 'Separator','on');
	setLineWidth(item_lw,cb_LineWidth)
	item_ls = uimenu(cmenuHand, 'Label', 'Line Style');
	setLineStyle(item_ls,{cbls1 cbls2 cbls3 cbls4})
	item_lc = uimenu(cmenuHand, 'Label', 'Color');
	setLineColor(item_lc,cb_color)

% -----------------------------------------------------------------------------------------
function call_gmtedit(obj,eventdata,h)
	pt = get(gca, 'CurrentPoint');
	gmtedit(getappdata(h,'FullName'),['-P' sprintf('%.6f',pt(1,1)) '/' sprintf('%.6f',pt(1,2))]);

% -----------------------------------------------------------------------------------------
function cb = uictx_setMarker(h,prop)
	% Set uicontext colors in a PB object class hose handles are contained in h
	% PROP is either "Marker" or "MarkerSize". OPT is either the symbol or it's size
	if (strcmp(prop,'Marker'))
        cb{1} = {@other_Marker,h,prop,'+'};       cb{2} = {@other_Marker,h,prop,'o'};
        cb{3} = {@other_Marker,h,prop,'*'};       cb{4} = {@other_Marker,h,prop,'.'};
        cb{5} = {@other_Marker,h,prop,'x'};       cb{6} = {@other_Marker,h,prop,'s'};
        cb{7} = {@other_Marker,h,prop,'d'};       cb{8} = {@other_Marker,h,prop,'^'};
        cb{9} = {@other_Marker,h,prop,'v'};       cb{10} = {@other_Marker,h,prop,'>'};
        cb{11} = {@other_Marker,h,prop,'<'};       cb{12} = {@other_Marker,h,prop,'p'};
        cb{13} = {@other_Marker,h,prop,'h'};
	elseif (strcmp(prop,'MarkerSize'))
        cb{1} = {@other_Marker,h,prop,7};       cb{2} = {@other_Marker,h,prop,8};
        cb{3} = {@other_Marker,h,prop,9};       cb{4} = {@other_Marker,h,prop,10};
        cb{5} = {@other_Marker,h,prop,12};      cb{6} = {@other_Marker,h,prop,14};
        cb{7} = {@other_Marker,h,prop,16};      cb{8} = {@other_SymbSize,h};
	end
	
	function other_Marker(obj,eventdata,h,prop,opt)
	set(h,prop,opt);    refresh
% -----------------------------------------------------------------------------------------

% -----------------------------------------------------------------------------------------
function cb = uictx_Class_LineWidth(h)
% Set uicontext LineWidths in a PB object class hose handles are contained in h
cb{1} = {@other_Class_LineWidth,h,1};      cb{2} = {@other_Class_LineWidth,h,2};
cb{3} = {@other_Class_LineWidth,h,3};      cb{4} = {@other_Class_LineWidth,h,4};
cb{5} = {@other_Class_LineWidth,h,[]};

function other_Class_LineWidth(obj,eventdata,h,opt)
% If individual Lines were previously removed (by "Remove Line") h has invalid
% handles, so make sure all handles are valid
h=h(ishandle(h));
if ~isempty(opt)
    set(h,'LineWidth',opt);        refresh;
else
    prompt = {'Enter new line width (pt)'};     dlg_title = 'Line width';
    num_lines= [1 30];
    resp  = inputdlg(prompt,dlg_title,num_lines);
    if isempty(resp);    return;     end
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
    if isempty(prop)    % line
        set(h,'color',cor);   refresh;
    else                % marker
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
    else
        return
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
	h=h(ishandle(h));
	set(h,'LineStyle',opt);        refresh;
% -----------------------------------------------------------------------------------------

% -----------------------------------------------------------------------------------------
function set_greatCircle_uicontext(h)
	% h is a handle to a great circle arc (in geog coords) object
	handles = guidata(h(1));	cmenuHand = uicontextmenu('Parent',handles.figure1);
	set(h, 'UIContextMenu', cmenuHand);
	cb_LineWidth = uictx_LineWidth(h);      % there are 5 cb_LineWidth outputs
	cb_solid  = 'set(gco, ''LineStyle'', ''-''); refresh';   cb_dashed      = 'set(gco, ''LineStyle'', ''--''); refresh';
	cb_dotted = 'set(gco, ''LineStyle'', '':''); refresh';   cb_dash_dotted = 'set(gco, ''LineStyle'', ''-.''); refresh';
	cb_color = uictx_color(h);      % there are 9 cb_color outputs
	
	uimenu(cmenuHand, 'Label', 'Delete', 'Callback', 'delete(gco)');
	uimenu(cmenuHand, 'Label', 'Save line', 'Callback', {@save_formated,h});
	uimenu(cmenuHand, 'Label', 'Line length', 'Callback', {@show_LineLength,[],'total'});
	item_lw = uimenu(cmenuHand, 'Label', 'Line Width', 'Separator','on');
	setLineWidth(item_lw,cb_LineWidth)
	item_ls = uimenu(cmenuHand, 'Label', 'Line Style');
	setLineStyle(item_ls,{cb_solid cb_dashed cb_dotted cb_dash_dotted})
	item_lc = uimenu(cmenuHand, 'Label', 'Line Color');
	setLineColor(item_lc,cb_color)

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
	cb_color = uictx_color(h);      % there are 9 cb_color outputs
	% cb_MoveCircle        = {@move_circle,h};
	% cb_ChangeCircCenter1 = {@change_CircCenter1,h};
	cb_roi = 'mirone(''DrawClosedPolygon_CB'',guidata(gcbo),gco)';
	
	uimenu(cmenuHand, 'Label', 'Delete', 'Callback', 'delete(gco)');
	uimenu(cmenuHand, 'Label', 'Save circle', 'Callback', {@save_formated,h});
	uimenu(cmenuHand, 'Label', 'Line length', 'Callback', {@show_LineLength,[]});
	% item_MoveCenter = uimenu(cmenuHand, 'Label', 'Move (interactive)', 'Callback', cb_MoveCircle);
	% item_SetCenter0 = uimenu(cmenuHand, 'Label', 'Change');
	% item_SetCenter1 = uimenu(item_SetCenter0, 'Label', 'By coordinates', 'Callback', cb_ChangeCircCenter1);
	if ~strcmp(tag,'CircleEuler')       % "Just" a regular geographical circle
        uimenu(cmenuHand, 'Label', 'Region-Of-Interest', 'Separator','on', 'Callback', cb_roi);
	else
        uimenu(cmenuHand, 'Label', 'Compute velocity', 'Separator','on', 'Callback', {@compute_EulerVel,h});
	end
	item_lw = uimenu(cmenuHand, 'Label', 'Line Width', 'Separator','on');
	cb_LineWidth = uictx_LineWidth(h);      % there are 5 cb_LineWidth outputs
	setLineWidth(item_lw,cb_LineWidth)
	item_ls = uimenu(cmenuHand, 'Label', 'Line Style');
	setLineStyle(item_ls,{cb_solid cb_dashed cb_dotted cb_dash_dotted})
	item_lc = uimenu(cmenuHand, 'Label', 'Line Color');
	setLineColor(item_lc,cb_color)

% -----------------------------------------------------------------------------------------
function set_circleCart_uicontext(h)
	% h is a handle to a circle (in cartesian coords) object
	handles = guidata(h(1));	cmenuHand = uicontextmenu('Parent',handles.figure1);
	set(h, 'UIContextMenu', cmenuHand);
	cb_solid  = 'set(gco, ''LineStyle'', ''-''); refresh';   cb_dashed      = 'set(gco, ''LineStyle'', ''--''); refresh';
	cb_dotted = 'set(gco, ''LineStyle'', '':''); refresh';   cb_dash_dotted = 'set(gco, ''LineStyle'', ''-.''); refresh';
	cb_color = uictx_color(h);      % there are 9 cb_color outputs
	cb_roi = 'mirone(''DrawClosedPolygon_CB'',guidata(gcbo),gco)';
	
	uimenu(cmenuHand, 'Label', 'Delete', 'Callback', 'delete(gco)');
	uimenu(cmenuHand, 'Label', 'Save circle', 'Callback', {@save_formated,h});
	uimenu(cmenuHand, 'Label', 'Circle perimeter', 'Callback', {@show_LineLength,[]});
	uimenu(cmenuHand, 'Label', 'Move (interactive)', 'Callback', {@move_circle,h});
	item_SetCenter0 = uimenu(cmenuHand, 'Label', 'Change');
	uimenu(item_SetCenter0, 'Label', 'By coordinates', 'Callback', {@change_CircCenter1,h});
	uimenu(cmenuHand, 'Label', 'Region-Of-Interest', 'Separator','on', 'Callback', cb_roi);
	item_lw = uimenu(cmenuHand, 'Label', 'Line Width', 'Separator','on');
	cb_LineWidth = uictx_LineWidth(h);      % there are 5 cb_LineWidth outputs
	setLineWidth(item_lw,cb_LineWidth)
	item2 = uimenu(cmenuHand, 'Label', 'Line Style');
	setLineStyle(item2,{cb_solid cb_dashed cb_dotted cb_dash_dotted})
	item3 = uimenu(cmenuHand, 'Label', 'Line Color');
	setLineColor(item3,cb_color)

% -----------------------------------------------------------------------------------------
function compute_EulerVel(obj,eventdata,h)
% alat & alon are the point coords. plat, plon & omega are the pole parameters
D2R = pi/180;
earth_rad = 6371e3;    % Earth radius in km
s = get(h,'Userdata');
plat = s.clat*D2R;      plon = s.clon*D2R;      omega = s.omega;
alat = s.rlat*D2R;      alon = s.rlon*D2R;

x = cos(plat)*sin(plon)*sin(alat) - cos(alat)*sin(alon)*sin(plat);    % East vel
y = cos(alat)*cos(alon)*sin(plat) - cos(plat)*cos(plon)*sin(alat);    % North vel
z = cos(plat)*cos(alat)*sin(alon-plon);
vlon = -sin(alon)*x + cos(alon)*y;
vlat = -sin(alat)*cos(alon)*x-sin(alat)*sin(alon)*y + cos(alat)*z;
azim = 90 - atan2(vlat,vlon) / D2R;

if (azim < 0)       % Give allways the result in the 0-360 range
    azim = azim + 360;
end

x = sin(alat)*sin(plat) + cos(alat)*cos(plat)*cos(plon-alon);
delta = acos(x);
vel = omega*D2R/1e+4 * earth_rad * sin(delta);      % to give velocity in cm/Ma

msg = ['Pole name:  ', s.plates, sprintf('\n'), ...
        'Pole lon = ', sprintf('%3.3f',s.clon), sprintf('\n'), ...
        'Pole lat = ', sprintf('%2.3f',s.clat), sprintf('\n'), ...
        'Pole rate = ', sprintf('%.3f',omega), sprintf('\n'), ...
        'At point: ',sprintf('\n'), ...
        'Lon = ', sprintf('%3.3f',s.rlon), sprintf('\n'), ...
        'Lat = ', sprintf('%2.3f',s.rlat), sprintf('\n'), ...
        'Speed (cm/yr) =   ', sprintf('%2.2f',vel), sprintf('\n'), ...
        'Azimuth (degrees cw from North) = ', sprintf('%3.1f',azim)];
msgbox(msg,'Euler velocity')

% -----------------------------------------------------------------------------------------
function set_vector_uicontext(h)
	% h is a handle to a vector object
	handles = guidata(h(1));	cmenuHand = uicontextmenu('Parent',handles.figure1);
	set(h, 'UIContextMenu', cmenuHand);
	uimenu(cmenuHand, 'Label', 'Delete', 'Callback', {@delete_vector,h});
	%set(h,'ButtonDownFcn','selectmoveresize')

% -----------------------------------------------------------------------------------------
function delete_vector(obj,eventdata,h)
    % Brute force delete a vector whose handle is h. delete(gco) just doesn't kill the head
    try,	delete(h),	end		% Use a try because a something else may have deleted one part only

% -----------------------------------------------------------------------------------------
function fill_Polygon(obj,eventdata,h)
	% Turn a closed polygon into a patch and fill it in light gray by default
	% EXPERIMENTAL CODE. NOT IN USE.
	x = get(h,'XData');      y = get(h,'YData');
	patch(x,y,0,'FaceColor',[.7 .7 .7], 'EdgeColor','k');

% -----------------------------------------------------------------------------------------
function show_swhatRatio(obj,eventdata,h)
    msgbox(['Swath Ratio for this track is: ' sprintf('%g',getappdata(h,'swathRatio'))],'')

% -----------------------------------------------------------------------------------------
function show_Area(obj,eventdata,h)
	% Compute area under line and insult the user if the line is not closed
	% NOTE that H is optional. Use only when want to make sure that this fun
	% uses that handle (does notwork with copyied objects)

	if (nargin == 3)
		if (size(h,1) >= 2 && size(h,2) == 2)
			x = h(:,1);     y = h(:,2);
		elseif (ishandle(h))
			x = get(h,'XData');    y = get(h,'YData');
		end
	elseif (nargin == 2 || isempty(h) || length(h) > 1)
		h = gco;
		x = get(h,'XData');    y = get(h,'YData');
	end

	% Contour lines for example have NaNs and not at the same x,y positions (???)
	ix = isnan(x);
	x(ix) = [];             y(ix) = [];
	iy = isnan(y);
	x(iy) = [];             y(iy) = [];
	handles = guidata(get(0,'CurrentFigure'));
	if ~( (x(1) == x(end)) && (y(1) == y(end)) )
        msg{1} = 'This is not a closed line. Result is therefore probably VERY idiot';
	else
        msg{1} = '';
	end
	if (handles.geog)
        area = area_geo(y,x);    % Area is reported on the unit sphere
        area = area * 4 * pi * (handles.EarthRad^2);
        msg{2} = ['Area = ' sprintf('%g',area) ' km^2'];
        msgbox(msg,'Area')
	else
        area = polyarea(x,y);   % Area is reported in map user unites
        msg{2} = ['Area = ' sprintf('%g',area) ' map units ^2'];
        msgbox(msg,'Area')
	end
	refresh

% -----------------------------------------------------------------------------------------
function ll = show_LineLength(obj,eventdata,h, opt)
% Line length (perimeter if it is a closed polyline). If output argument, return a structure
% ll.len and ll.type, where "len" is line length and "type" is either 'geog' or 'cart'.
% For polylines ll.len contains only the total length.
% 22-09-04  Added OPT option. If it exists, report only total length (for nargout == 0)
% 22-10-05  H is now only to be used if we whant to specificaly use that handle. Otherwise use []
% to fish it with gco (MUST use this form to work with copied objects)
% 16-08-07  H can contain a Mx2 column vector with the line vertices.

	n_args = nargin;
	if (n_args <= 3),   opt = [];   end
	if (n_args == 3)
        if (size(h,1) >= 2 && size(h,2) == 2)
            x = h(:,1);     y = h(:,2);
        elseif (ishandle(h))
            x = get(h,'XData');    y = get(h,'YData');
        end
	elseif (n_args == 2 || isempty(h) || length(h) > 1)
        h = gco;
        x = get(h,'XData');    y = get(h,'YData');
	end

msg = [];               handles = guidata(get(0,'CurrentFigure'));

% Contour lines for example have NaNs and not at the same x,y positions (???)
ix = isnan(x);      x(ix) = [];     y(ix) = [];
iy = isnan(y);      x(iy) = [];     y(iy) = [];
if (handles.geog)
    lat_i = y(1:length(y)-1);   lat_f = y(2:length(y));     clear y;
    lon_i = x(1:length(x)-1);   lon_f = x(2:length(x));     clear x;
	tmp = vdist(lat_i,lon_i,lat_f,lon_f,handles.DefineEllipsoide([1 3]));
    
    switch handles.DefineMeasureUnit
        case 'n'        % Nautical miles
            scale = 1852;   str_unit = ' NM';
        case 'k'        % Kilometers
            scale = 1000;   str_unit = ' kilometers';
        case {'m','u'}  % Meters or user unites
            scale = 1;   str_unit = ' meters(?)';
    end
    total_len = sum(tmp) / scale;
    if (nargout == 0 && isempty(opt))
        len_i = zeros(1,length(tmp));
        for i = 1:length(tmp)
            len_i(i) = tmp(i) / scale;
            msg = [msg; {['Length' sprintf('%d',i) '  =  ' sprintf('%.5f',len_i(i)) str_unit]}];
        end
        if (length(msg) < 20)
            msg = [msg; {['Total length = ' sprintf('%.5f',total_len) str_unit]}];
        else
            msg = {['Total length = ' sprintf('%.5f',total_len) str_unit]};
        end
        msgbox(msg,'Line(s) length')
    elseif (nargout == 0 && ~isempty(opt))
        msgbox(['Total length = ' sprintf('%.5f',total_len) str_unit],'Line length')
    else        % Should we also out output also the partial lengths?
        ll.len = total_len;   ll.type = 'geog';
    end
else
    dx = diff(x);   dy = diff(y);
    total_len = sum(sqrt(dx.*dx + dy.*dy));
    if (nargout == 0 && isempty(opt))
        len_i = zeros(1,length(dx));
        for i = 1:length(dx)
            len_i(i) = sqrt(dx(i)^2 + dy(i)^2);
            msg = [msg; {['Length' sprintf('%d',i) '  =  ' sprintf('%.5f',len_i(i)) ' map units']}];
        end
        if (length(msg) < 20)
            msg = [msg; {['Total length = ' sprintf('%.5f',total_len) ' map units']}];
        else
            msg = {['Total length = ' sprintf('%.5f',total_len) ' map units']};
        end
        msgbox(msg,'Line(s) length')
    elseif (nargout == 0 && ~isempty(opt))
        msgbox(['Total length = ' sprintf('%.5f',total_len) ' map units'],'Line length')
    else        % The same question as in the geog case
        ll.len = total_len;   ll.type = 'cart';
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
        msgbox(['Total tracks length = ' sprintf('%g',len) ' NM'])
	end

% -----------------------------------------------------------------------------------------
function azim = show_lineAzims(obj,eventdata,h)
	% Works either in geog or cart coordinates. Otherwise the result is a non-sense
	% If output argument, return a structure % azim.az and azim.type, where "len" is line
	% azimuth and "type" is either 'geog' or 'cart'.
	% 22-10-05  H is now only to be used if we whant to specificaly use that handle. Otherwise use
	% either [] or don't pass the H argument to fish it with gco (MUST use this form to work with copied objects)
	% 16-08-07  H can contain a Mx2 column vector with the line vertices.

	if (nargin == 3)
        if (size(h,1) >= 2 && size(h,2) == 2)
            x = h(:,1);     y = h(:,2);
        elseif (ishandle(h))
            x = get(h,'XData');    y = get(h,'YData');
        end
	elseif (nargin == 2 || isempty(h) || length(h) > 1)
        h = gco;
        x = get(h,'XData');    y = get(h,'YData');
	end
	
	handles = guidata(get(0,'CurrentFigure'));
	if (handles.geog)
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
	msg = cell(numel(az),1);
	for (i = 1:numel(az))
        if (nargout == 0)
            msg{i} = ['Azimuth' sprintf('%d',i) '  =  ' sprintf('%3.1f',az(i)) '  degrees'];
        else
            azim.az = az;
        end
	end
	if (nargout == 0)
        if (numel(az) > 1)
            msg{end+1} = '';
            id = (az > 270);    az(id) = az(id) - 360;
            az_mean = mean(az);
            msg{end+1} = ['Mean azimuth = ' sprintf('%.1f',az_mean) '  degrees'];
        end
        if (numel(az) > 15),   msg = msg(end-2:end);   end
        msgbox(msg,'Line(s) Azimuth')
	end
	refresh

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
if (length(c) > 1),         % That is, if a color was selected
    if ~isempty(opt) && ischar(opt)
        set(h,opt,c);   refresh;
    else
        set(h,'color',c);   refresh;
    end
else
    return;
end
% -----------------------------------------------------------------------------------------

% -----------------------------------------------------------------------------------------
function cb = uictx_LineWidth(h)
	% Set uicontext colors in object hose handle is gco (or h for "other color")
	cb{1} = 'set(gco, ''LineWidth'', 1);refresh';   cb{2} = 'set(gco, ''LineWidth'', 2);refresh';
	cb{3} = 'set(gco, ''LineWidth'', 3);refresh';   cb{4} = 'set(gco, ''LineWidth'', 4);refresh';
	cb{5} = {@other_LineWidth,h};

function other_LineWidth(obj,eventdata,h)
	prompt = {'Enter new line width (pt)'};     dlg_title = 'Line width';
	num_lines= [1 30];
	resp  = inputdlg(prompt,dlg_title,num_lines);
	if isempty(resp);    return;     end
	set(h,'LineWidth',str2double(resp));        refresh
% -----------------------------------------------------------------------------------------

% -----------------------------------------------------------------------------------------
function hVec = draw_vector
    hFig = get(0,'CurrentFigure');          handles = guidata(hFig);
	hVec(1) = line('XData', [], 'YData', [],'Color',handles.DefLineColor,'LineWidth',handles.DefLineThick,'Tag','Arrow');
	hVec(2) = patch('XData', [], 'YData', [],'FaceColor',handles.DefLineColor,'EdgeColor',handles.DefLineColor,'Tag','Arrow');
	state = uisuspend_fig(hFig);        % Remember initial figure state
	set(hFig,'Pointer', 'crosshair');
	w = waitforbuttonpress;
	if w == 0       % A mouse click
        vectorFirstButtonDown(hFig,handles.axes1,hVec,state)
	else
        set(hFig,'Pointer', 'arrow');    hVec = [];
	end

function vectorFirstButtonDown(hFig,hAxes,h,state)
    pt = get(hAxes, 'CurrentPoint');
    set(hFig,'WindowButtonMotionFcn',{@wbm_vector,[pt(1,1) pt(1,2)],h,hAxes},'WindowButtonDownFcn',{@wbd_vector,h,state});

function wbm_vector(obj,eventdata,origin,h,hAxes)
	pt = get(hAxes, 'CurrentPoint');
	x  = [origin(1) pt(1,1)];   y = [origin(2) pt(1,2)];
	dx = diff(x);               dy = diff(y);
	set(h(1),'XData',x, 'YData',y)
	ax = getappdata(hAxes,'ThisImageLims');
	lx = diff(ax(1:2))*25e-3;   ly = diff(ax(3:4))*10e-3;
	phi = atan2(dy,dx);
	head = rotate2([-lx 0 -lx; ly 0 -ly], [0; 0], phi);
	set(h(2),'XData',head(1,:)+pt(1,1), 'YData',head(2,:)+pt(1,2));

function wbd_vector(obj,eventdata,h,state)
    uirestore_fig(state);           % Restore the figure's initial state
    ui_edit_polygon(h(2))

function newpoints = rotate2(points,orig,phi)
	%ROTATE2  rotate points PHI radians around the ORIG point
	A=[cos(phi) -sin(phi); sin(phi) cos(phi)];
	newpoints=A*points+orig*ones(1,size(points,2));
% -----------------------------------------------------------------------------------------

% -----------------------------------------------------------------------------------------
function h_gcirc = draw_greateCircle
    hFig = get(0,'CurrentFigure');          handles = guidata(hFig);
    h_gcirc = line('XData', [], 'YData', [],'Color',handles.DefLineColor,'LineWidth',handles.DefLineThick);
	state = uisuspend_fig(hFig);     % Remember initial figure state
	set(hFig,'Pointer', 'crosshair'); % to avoid the compiler BUG
	w = waitforbuttonpress;                             %
	if w == 0       % A mouse click                     %
        gcircFirstButtonDown(hFig,h_gcirc,state)             %
	else                                                %
        set(hFig,'Pointer', 'arrow'); %
        h_gcirc = [];                                   %
	end                                                 %
%---------------
function gcircFirstButtonDown(hFig,h,state)
	hAxes = get(hFig,'CurrentAxes');	pt = get(hAxes, 'CurrentPoint');
	set(hFig,'WindowButtonMotionFcn',{@wbm_gcircle,[pt(1,1) pt(1,2)],h,hFig,hAxes},'WindowButtonDownFcn',{@wbd_gcircle,h,state});
%---------------
function wbm_gcircle(obj,eventdata,first_pt,h,hFig,hAxes)
	pt = get(hAxes, 'CurrentPoint');
	[x,y] = gcirc(first_pt(1),first_pt(2),pt(1,1),pt(1,2));
	% Find the eventual Date line discontinuity and insert a NaN on it
	% ind = find(abs(diff(x)) > 100);   % 100 is good enough
	% if (~isempty(ind))
	%     if (length(ind) == 2)
	%         x = [x(1:ind(1)) NaN x(ind(1)+1:ind(2)) NaN x(ind(2)+1:end)];
	%         y = [y(1:ind(1)) NaN y(ind(1)+1:ind(2)) NaN y(ind(2)+1:end)];
	%     elseif (length(ind) == 1)
	%         x = [x(1:ind) NaN x(ind+1:end)];   y = [y(1:ind) NaN y(ind+1:end)];
	%     end
	% end
	set(h, 'XData', x, 'YData', y,'Userdata',[first_pt [x y]]);
%---------------
function wbd_gcircle(obj,eventdata,h,state)
	lons_lats = get(h,'UserData');    setappdata(h,'LonLatRad',lons_lats)   % save this in appdata
	set(h,'Tag','GreatCircle')
	uirestore_fig(state);           % Restore the figure's initial state
% -----------------------------------------------------------------------------------------

% -----------------------------------------------------------------------------------------
function h_circ = draw_circleGeo
	% THIS IS NOW ONLY USED NOW WITH CARTESIAN CIRCLES
	% Given one more compiler BUG, (WindowButtonDownFcn cannot be redefined)
	% I found the following workaround.
	hFig = get(0,'CurrentFigure');          handles = guidata(hFig);
	h_circ = line('XData', [], 'YData', [],'Color',handles.DefLineColor,'LineWidth',handles.DefLineThick);
	%set(hFig,'WindowButtonDownFcn',{@circFirstButtonDown,h_circ}, 'Pointer', 'crosshair');
	state = uisuspend_fig(hFig);     % Remember initial figure state
	set(hFig,'Pointer', 'crosshair'); % to avoid the compiler BUG
	w = waitforbuttonpress;                             %
	if w == 0       % A mouse click                     %
        circFirstButtonDown(h_circ,state)               %
	else                                                %
        set(get(0,'CurrentFigure'),'Pointer', 'arrow'); %
        h_circ = [];                                    %
	end                                                 %

%---------------
%function circFirstButtonDown(obj,eventdata,h)      % For non compiled version
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
	lon_lat_rad = get(h,'UserData');    setappdata(h,'LonLatRad',lon_lat_rad)   % save this in appdata
	set(h,'Tag','circleCart')
    rmappdata(h,'X');           rmappdata(h,'Y');
	uirestore_fig(state);           % Restore the figure's initial state
% -----------------------------------------------------------------------------------------

% -----------------------------------------------------------------------------------------
function h_circ = draw_circleEulerPole(lon,lat)
	% Draw a circle (or arc of a circle) about the Euler Pole (or any other origin)
	% See notes above for the reason why waitforbuttonpress is used.
	h_circ = line('XData', [], 'YData', []);
	hFig = get(0,'CurrentFigure');         hAxes = get(hFig,'CurrentAxes');
	set(hFig,'Pointer', 'crosshair');
	w = waitforbuttonpress;
	if w == 0       % A mouse click
        set(hFig,'WindowButtonMotionFcn',{@wbm_circle,[lon lat],h_circ,hAxes},'WindowButtonDownFcn',{@wbd_circle,h_circ});
	else
        set(hFig,'Pointer', 'arrow');
        h_circ = [];
	end

% -----------------------------------------------------------------------------------------
function move_circle(obj,eventdata,h)
% ONLY FOR CARTESIAN CIRCLES.
hFig = get(0,'CurrentFigure');  hAxes = get(hFig,'CurrentAxes');
state = uisuspend_fig(hFig);      % Remember initial figure state
np = numel(get(h,'XData'));     x = linspace(-pi,pi,np);
setappdata(h,'X',cos(x));       setappdata(h,'Y',sin(x))    % Save unit circle coords
center = getappdata(h,'LonLatRad');
set(hFig,'WindowButtonMotionFcn',{@wbm_MoveCircle,h,center,hAxes},...
    'WindowButtonDownFcn',{@wbd_MoveCircle,h,state,hAxes},'Pointer', 'crosshair');

function wbm_MoveCircle(obj,eventdata,h,center,hAxes)
	pt = get(hAxes, 'CurrentPoint');
	x = getappdata(h,'X');          y = getappdata(h,'Y');
	x = pt(1,1) + center(3)*x;      y = pt(1,2) + center(3)*y;
	set(h, 'XData', x, 'YData', y,'Userdata',[pt(1,1) pt(1,2) center(3)]);

function wbd_MoveCircle(obj,eventdata,h,state,hAxes)
	% check if x,y is inside of axis
	pt = get(hAxes, 'CurrentPoint');  x = pt(1,1);    y = pt(1,2);
	x_lim = get(hAxes,'xlim');        y_lim = get(hAxes,'ylim');
	if (x<x_lim(1)) || (x>x_lim(2)) || (y<y_lim(1)) || (y>y_lim(2));   return; end
	lon_lat_rad = get(h,'UserData');    setappdata(h,'LonLatRad',lon_lat_rad)   % save this in appdata
    rmappdata(h,'X');           rmappdata(h,'Y');
	uirestore_fig(state);           % Restore the figure's initial state

% -----------------------------------------------------------------------------------------
function change_CircCenter1(obj,eventdata,h)
	% Change the Circle's center by asking it's coordinates
	% ONLY FOR CARTESIAN CIRCLES.
	lon_lat_rad = getappdata(h,'LonLatRad');
	prompt = {'Enter new lon (or x)' ,'Enter new lat (or y)', 'Enter new radius'};     dlg_title = 'Change circle';
	num_lines= [1 30; 1 30; 1 30];
	def = {num2str(lon_lat_rad(1)) num2str(lon_lat_rad(2)) num2str(lon_lat_rad(3))};
	resp  = inputdlg(prompt,dlg_title,num_lines,def);
	if isempty(resp);    return;     end
    np = numel(get(h,'XData'));     x = linspace(-pi,pi,np);
	y = sin(x);                     x = cos(x);    % unit circle coords
	x = str2double(resp{1}) + str2double(resp{3}) * x;
	y = str2double(resp{2}) + str2double(resp{3}) * y;
	set(h, 'XData', x, 'YData', y);
	setappdata(h,'LonLatRad',[str2double(resp{1}) str2double(resp{2}) str2double(resp{3})])
% -----------------------------------------------------------------------------------------

% -----------------------------------------------------------------------------------------
function rectangle_limits(obj,eventdata)
% Change the Rectangle's limits by asking it's corner coordinates
% res = check_IsRectangle(h);
% if ~res
%     errordlg('This no longer a rectangle. The edition deformed it','Error');     return;
% end
%if (isempty(h) | nargin == 2)   h = gco;    end
h = gco;
x = get(h,'XData');     y = get(h,'YData');

region = bg_region('with_limits',[x(1) x(3) y(1) y(3)]);
if isempty(region),    return;  end     % User gave up
x_min = region(1);      x_max = region(2);
y_min = region(3);      y_max = region(4);

set(h, 'XData', [x_min,x_min,x_max,x_max,x_min], 'YData', [y_min,y_max,y_max,y_min,y_min]);
% -----------------------------------------------------------------------------------------

% -----------------------------------------------------------------------------------------
function rectangle_register_img(obj,event)
	% Prompt user for rectangle corner coordinates and use them to register the image
	h = gco;
	handles = guidata(get(h,'Parent'));
	rect_x = get(h,'XData');   rect_y = get(h,'YData');		% Get rectangle limits

	region = bg_region('empty');
	if isempty(region),    return;  end     % User gave up
	x_min = region(1);      x_max = region(2);
	y_min = region(3);      y_max = region(4);
	handles.geog = aux_funs('guessGeog',region(1:4));		% Trast more in the test here
	ax = handles.axes1;

	x(1) = rect_x(1);     x(2) = rect_x(2);     x(3) = rect_x(3);
	y(1) = rect_y(1);     y(2) = rect_y(2);
	img = get(handles.hImg,'CData');
	% Transform the ractangle limits into row-col limits
	limits = getappdata(handles.axes1,'ThisImageLims');
	r_c = cropimg(limits(1:2), limits(3:4), img, [x(1) y(1) (x(3)-x(2)) (y(2)-y(1))], 'out_precise');
	% Find if we are dealing with a image with origin at upper left (i.e. with y positive down)
	if(strcmp(get(ax,'YDir'),'reverse'))
		img = flipdim(img,1);
		% We have to invert the row count to account for the new origin in lower left corner
		tmp = r_c(1);
		r_c(1) = size(img,1) - r_c(2) + 1;
		r_c(2) = size(img,1) - tmp + 1;
	end
	% Compute and apply the affine transformation
	base  = [x_min y_min; x_min y_max; x_max y_max; x_max y_min];
	input = [r_c(3) r_c(1); r_c(3) r_c(2); r_c(4) r_c(2); r_c(4) r_c(1)];
	% tform = cp2tform(input,base,'affine');
	% [new_xlim,new_ylim] = tformfwd(tform,[1 size(img,2)],[1 size(img,1)]);

	trans = AffineTransform(input,base);
	x_pt = [1; size(img,2)];    y_pt = [1; size(img,1)];    % For more X points, change accordingly
	X1 = [x_pt y_pt ones(size(x_pt,1),1)];
	U1 = X1 * trans;
	new_xlim = U1(:,1)';        new_ylim = U1(:,2)';

	% Rebuild the image with the new limits.
	% Here, I found no way to just update the image with the new limits.
	% The command: image(new_xlim,new_ylim,img); just added a new image (CData) to the
	% existing one (that's not what the manual says). So I have to kill the old image
	% and redraw it again. That's f... stupid.
	[m,n,k] = size(img);
	[new_xlim,new_ylim] = aux_funs('adjust_lims',new_xlim,new_ylim,m,n);
	delete(handles.hImg);
	handles.hImg = image(new_xlim,new_ylim,img,'Parent',handles.axes1);
	set(ax,'xlim',new_xlim,'ylim',new_ylim,'YDir','normal')
	handles.head(1:4) = [new_xlim new_ylim];
	resizetrue(handles, [], 'xy');
	setappdata(ax,'ThisImageLims',[get(ax,'XLim') get(ax,'YLim')])
	handles.old_size = get(handles.figure1,'Pos');      % Save fig size to prevent maximizing
	handles.origFig = img;

	% [m,n,k] = size(img);
	% [new_xlim,new_ylim] = aux_funs('adjust_lims',new_xlim,new_ylim,m,n);
	% set(gca,'XLim',new_xlim,'YLim',new_ylim,'YDir','normal')
	% set(h_img,'XData',new_xlim,'YData',new_ylim)
	% resizetrue(handles, [])
	% x = [x_min x_min x_max x_max x_min];        y = [y_min y_max y_max y_min y_min];
	% set(h,'XData',x,'YData',y)

	% Redraw the rectangle that meanwhile has gone to the ether togheter with gca.
	lt = handles.DefLineThick;  lc = handles.DefLineColor;
	x = [x_min x_min x_max x_max x_min];        y = [y_min y_max y_max y_min y_min];
	h = line('XData',x,'YData',y,'Color',lc,'LineWidth',lt,'Parent',handles.axes1);
	if (handles.image_type == 2)                    % Lets pretend that we have a GeoTIFF image
		handles.image_type = 3;
		Hdr.LL_prj_xmin = new_xlim(1);      Hdr.LR_prj_xmax = new_xlim(2);
		Hdr.LL_prj_ymin = new_ylim(1);      Hdr.UR_prj_ymax = new_ylim(2);
		Hdr.projection = 'linear';          Hdr.datum = 'unknown';
		set(handles.figure1,'UserData',Hdr);        % Minimalist Hdr to allow saving as a GeoTIFF image
	end
	x_inc = (new_xlim(2)-new_xlim(1)) / (size(img,2) - 1);
	y_inc = (new_ylim(2)-new_ylim(1)) / (size(img,1) - 1);
	%handles.head = [new_xlim(1) new_xlim(2) new_ylim(1) new_ylim(2) 0 255 0 x_inc y_inc];     % TEMP and ...
	handles.head(8:9) = [x_inc y_inc];   

	handles.fileName = [];			% Not loadable in session
	if (handles.validGrid)
		new_xlim = linspace(new_xlim(1),new_xlim(2),size(img,2));		new_ylim = linspace(new_ylim(1),new_ylim(2),size(img,1));
		setappdata(handles.figure1,'dem_x',new_xlim);  	setappdata(handles.figure1,'dem_y',new_ylim);
	end

	if (handles.geog)
		mirone('SetAxesNumericType',handles,[])          % Set axes uicontextmenus
	end
	guidata(handles.figure1, handles);
	draw_funs(h,'line_uicontext')       % Set lines's uicontextmenu

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
	U = uv;         % just solve for the first two columns of T

	% We know that X * T = U
	if rank(X) >= K
		Tinv = X \ U;
	else
		msg = 'At least %d non-collinear points needed to infer %s transform.';
		errordlg(sprintf(msg,K,'affine'),'Error');
	end

	Tinv(:,3) = [0 0 1]';       % add third column
	trans = inv(Tinv);
	trans(:,3) = [0 0 1]';

% -----------------------------------------------------------------------------------------
function Transplant_Image(obj,eventdata)
% Cirurgy Imagery operation. An external image will be inplanted inside the
% rectangular zone defined by the rectangle whose handle is h.
% Notice that we have to forsee the possibility of transplanting RGB images
% into indexed bg images and vice-versa.

h = gco;
hFig = get(0,'CurrentFigure');  hAxes = get(hFig,'CurrentAxes');
out = implanting_img(findobj(hFig,'Type','image'),h,get(hAxes,'xlim'),get(hAxes,'ylim'));
if isempty(out),   return;      end
h_img = findobj(get(hFig,'Children'),'Type','image');     % Get background image handle
zz = get(h_img,'CData');

% Find if Implanting image needs to be ud fliped
if(strcmp(get(hAxes,'XDir'),'normal') && strcmp(get(hAxes,'YDir'),'reverse'))
        flip = 0;
else    flip = 1;
end

[nl_ip,nc_ip,n_planes_ip] = size(out.ip_img);       % Get dimensions of implanting image
[nl_bg,nc_bg,n_planes_bg] = size(zz);               % Get dimensions of bg image
if (n_planes_ip == 3),  indexed_ip = 0;     else   indexed_ip = 1;     end
if (n_planes_bg == 3),  indexed_bg = 0;     else   indexed_bg = 1;     end

if (out.resizeIP)
    % We have to interpolate the Ip image to fit exactly with the rectangle dimensions.
    %nl_new = linspace(1,nl_ip,(out.r_c(2)-out.r_c(1)+1));
    %nc_new = linspace(1,nc_ip,(out.r_c(4)-out.r_c(3)+1));
    %[X,Y] = meshgrid(nc_new,nl_new);
    head = [1 nc_ip 1 nl_ip 0 255 0 1 1];
    opt_N = ['-N' num2str(out.r_c(4)-out.r_c(3)+1) '/' num2str(out.r_c(2)-out.r_c(1)+1)]; % option for grdsample
    if (~indexed_ip)                                % Implanting image is of RGB type
        for i=1:3
            %ZI(:,:,i) = interp2(double(out.ip_img(:,:,i)),X,Y,'*cubic');
            ZI(:,:,i) = grdsample_m(single(out.ip_img(:,:,i)),head,opt_N);
        end
    else
        if isempty(out.ip_cmap)
            errordlg('Implanting image has no colormap. Don''t know what to do.','Sorry');  return
        end
        %ZI = interp2(double(out.ip_img),X,Y,'*cubic');
        ZI = grdsample_m(single(out.ip_img),head,opt_N);
    end
    if (flip),   ZI = flipdim(ZI,1);    end
elseif (out.resizeIP == 10) % So pra nao funcionar (da erro na penultima linha)
    %nl_new = linspace(1,nl_bg,(out.bg_size_updated(1)));
    %nc_new = linspace(1,nc_bg,(out.bg_size_updated(2)));
    %[X,Y] = meshgrid(nc_new,nl_new);
    %if (~indexed_bg)                            % Background image is of RGB type
        %for (i=1:3)
            %zz(:,:,i) = interp2(double(zz(:,:,i)),X,Y,'*cubic');
            %zz(:,:,i) = grdsample_m(zz(:,:,i),head,opt_N);
        %end
    %else
        %zz = interp2(double(zz),X,Y,'*cubic');
        %zz = grdsample_m(zz,head,opt_N);
    %end
    %if (flip)    out.ip_img = flipdim(out.ip_img,1);    end
end

if (indexed_ip && ~indexed_bg)           % Implanting indexed image on a RGB bg image
    I = ind2rgb8(out.ip_img,out.ip_cmap);    % Transform implanting image to RGB
elseif (indexed_ip && indexed_bg)        % Shit, both ip & bg images are indexed. We have to RGB them
    zz = ind2rgb8(zz,colormap);
    I = ind2rgb8(out.ip_img,out.ip_cmap);
elseif (~indexed_ip && ~indexed_bg)      % Nice, nothing to do
elseif (~indexed_ip && indexed_bg)       % Implanting RGB image on a indexed bg image.
    zz = ind2rgb8(zz,colormap);      % Transform bg image to RGB
end

zz(out.r_c(1):out.r_c(2), out.r_c(3):out.r_c(4), :) = uint8(ZI);
set(h_img,'CData',zz)

% -----------------------------------------------------------------------------------------
function copy_text_object(obj,eventdata)
    copyobj(gco,gca);
    move_text([],[])

% -----------------------------------------------------------------------------------------
function move_text(obj,eventdata)
	h = gco;
    hFig = get(0,'CurrentFigure');  hAxes = get(hFig,'CurrentAxes');
	state = uisuspend_fig(hFig);     % Remember initial figure state
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
	uirestore_fig(state);           % Restore the figure's initial state
% -----------------------------------------------------------------------------------------

% -----------------------------------------------------------------------------------------
function rotate_text(obj,eventdata)
	prompt = {'Enter angle of rotation'};     dlg_title = '';
	num_lines= [1 30];
	resp  = inputdlg(prompt,dlg_title,num_lines);
	if isempty(resp);    return;     end
	h = gco;
	set(h,'Rotation',str2double(resp))
	refresh

% -----------------------------------------------------------------------------------------
function text_FontSize(obj,eventdata)
h = gco;
ft = uisetfont(h,'Change Font');
if (~isstruct(ft) && ft == 0), return;   end
set(h,'FontName',ft.FontName,'FontUnits',ft.FontUnits,'FontSize',ft.FontSize, ...
    'FontWeight',ft.FontWeight,'FontAngle',ft.FontAngle)
refresh

% -----------------------------------------------------------------------------------------
function export_text(obj,eventdata)
h = gco;
pos = get(h,'Position');    font = get(h,'FontName');      size = get(h,'FontSize');
str = get(h,'String');      angle = get(h,'Rotation');

handles = guidata(gcbo);        % I hope I don't get into troubles because of this!
cd(handles.work_dir)
[FileName,PathName] = uiputfile({ ...
    '*.txt;*.TXT', 'Text file (*.txt,*.TXT)'; '*.*', 'All Files (*.*)'}, 'Select Text File name');
cd(handles.home_dir);       % allways come home to avoid troubles
if isequal(FileName,0);     refresh;  return;     end
pause(0.01)
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
tag = get(h,'Tag');
if (length(h) == 1 && length(get(h,'Xdata')) > 1)
    more_than_one = 1;     % Flags that h points to a multi-vertice object
else
    more_than_one = 0;
end

handles = guidata(h(1));
cmenuHand = uicontextmenu('Parent',handles.figure1);	set(h, 'UIContextMenu', cmenuHand);
separator = 0;
this_not = 0;       % for class symbols "this_not = 1";
seismicity_options = 0;
tide_options = 0;

if strcmp(tag,'hotspot')    % Then DATA must be a structure containing name & age for each hotspot
    uimenu(cmenuHand, 'Label', 'Hotspot info', 'Callback', {@hotspot_info,h,data.name,data.age,[]});
    uimenu(cmenuHand, 'Label', 'Plot name', 'Callback', {@hotspot_info,h,data.name,data.age,'text'});
    separator = 1;
    this_not = 1;           % It is used for not seting some options inapropriate to class symbols
elseif strcmp(tag,'volcano')    % Then DATA must be a structure containing name, description & dating for each volcano
    uimenu(cmenuHand, 'Label', 'Volcano info', 'Callback', {@volcano_info,h,data.name,data.desc,data.dating});
    separator = 1;
    this_not = 1;           % It is used for not seting some options inapropriate to class symbols
elseif strcmp(tag,'ODP')    % Then DATA must be a structure with leg, site, z, & penetration for each site
    cb_ODPinfo = {@ODP_info,h,data.leg,data.site,data.z,data.penetration};
    uimenu(cmenuHand, 'Label', 'ODP info', 'Callback', cb_ODPinfo);
    separator = 1;
    this_not = 1;           % It is used for not seting some options inapropriate to class symbols
elseif strcmp(tag,'DSDP')   % Then DATA must be a structure with leg, site, z, & penetration for each site
    cb_ODPinfo = {@ODP_info,h,data.leg,data.site,data.z,data.penetration};
    uimenu(cmenuHand, 'Label', 'DSDP info', 'Callback', cb_ODPinfo);
    separator = 1;
    this_not = 1;           % It is used for not seting some options inapropriate to class symbols
elseif strcmp(tag,'City_major') || strcmp(tag,'City_other')
    this_not = 1;
elseif strcmp(tag,'Earthquakes')    % DATA is empty because I didn't store any info (they are too many)
    seismicity_options = isappdata(h,'SeismicityTime');
    this_not = 1;
elseif strcmp(tag,'Pointpolyline')  % DATA is empty because it doesn't have any associated info
    this_not = 0;
elseif strcmp(tag,'TTT')            % DATA is empty
    this_not = 0;
elseif strcmp(tag,'TideStation')    % DATA is empty
    tide_options = 1;
    separator = 0;
    this_not = 1;           % It is used for not seting some options inapropriate to class symbols
end

if (~this_not)   % non class symbols can be moved
    ui_edit_polygon(h)    % Set edition functions
    uimenu(cmenuHand, 'Label', 'Move (precise)', 'Callback', {@change_SymbPos,h});
end

if separator
    if (~more_than_one)         % Single symbol
        uimenu(cmenuHand, 'Label', 'Remove', 'Callback', 'delete(gco)', 'Separator','on');
    else                        % Multiple symbols
        uimenu(cmenuHand, 'Label', 'Remove this', 'Callback', {@remove_one_from_many,h}, 'Separator','on');
    end
else
    if (~more_than_one)         % Single symbol
        uimenu(cmenuHand, 'Label', 'Remove', 'Callback', 'delete(gco)');
    else                        % Multiple symbols
        uimenu(cmenuHand, 'Label', 'Remove this', 'Callback', {@remove_one_from_many,h});
    end
end
if (this_not)           % individual symbols don't belong to a class
    uimenu(cmenuHand, 'Label', 'Remove class', 'Callback', {@remove_symbolClass,h});
end
if (~this_not)          % class symbols don't export
    uimenu(cmenuHand, 'Label', 'Export', 'Callback', {@export_symbol,h});
    if (strcmp(tag,'Pointpolyline'))    % Allow pure grdtrack interpolation
        cbTrack = 'setappdata(gcf,''TrackThisLine'',gco); mirone(''ExtractProfile_CB'',guidata(gcbo),''point'')';
        uimenu(cmenuHand, 'Label', 'Point interpolation', 'Callback', cbTrack, 'Separator','on');
    end
end
if (seismicity_options)
    uimenu(cmenuHand, 'Label', 'Save events', 'Callback', 'save_seismicity(gcf,gco)', 'Separator','on');
    uimenu(cmenuHand, 'Label', 'Seismicity movie', 'Callback', 'animate_seismicity(gcf,gco)');
    uimenu(cmenuHand, 'Label', 'Draw polygon', 'Call', ...
        'mirone(''DrawClosedPolygon_CB'',guidata(gcbo),''SeismicityPolygon'')');
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
    uimenu(cmenuHand, 'Label', 'Plot tides', 'Callback', {@tidesStuff,h,'plot'}, 'Separator','on');
    uimenu(cmenuHand, 'Label', 'Station Info', 'Callback', {@tidesStuff,h,'info'});
    %uimenu(cmenuHand, 'Label', 'Tide Calendar', 'Callback', {@tidesStuff,h,'calendar'});
end
itemSymb = uimenu(cmenuHand, 'Label', 'Symbol', 'Separator','on');
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
cb_color = uictx_Class_LineColor(h,'MarkerFaceColor');              % there are 9 cb_PB_color outputs
itemFColor = uimenu(cmenuHand, 'Label', 'Fill Color');
setLineColor(itemFColor,cb_color)
cb_color = uictx_Class_LineColor(h,'MarkerEdgeColor');              % there are 9 cb_PB_color outputs
itemEColor = uimenu(cmenuHand, 'Label', 'Edge Color');
setLineColor(itemEColor,cb_color)

% -----------------------------------------------------------------------------------------
function move_symbol(obj,eventdata,h)
hFig = get(0,'CurrentFigure');  hAxes = get(hFig,'CurrentAxes');
state = uisuspend_fig(hFig);     % Remember initial figure state
pt = get(hAxes, 'CurrentPoint');
x = get(h,'XData');
if (length(x) > 1)              % We have a class symbol
    [c,i]=min(abs(x-pt(1,1)));  % Get the index of the current symbol
else                            % We have a individual symbol
    i = 1;
end
set(hFig,'WindowButtonMotionFcn',{@wbm_symb,h,i,hAxes},'WindowButtonDownFcn',{@wbd_symb,h,state,hAxes});

function wbm_symb(obj,eventdata,h,i,hAxes)
	pt = get(hAxes, 'CurrentPoint');
	pos_x = pt(1,1);        pos_y = pt(1,2);
	x = get(h,'XData');     x(i) = pos_x;
	y = get(h,'YData');     y(i) = pos_y;
	set(h,'XData',x,'YData',y);

function wbd_symb(obj,eventdata,h,state,hAxes)
	% check if x,y is inside of axis
	pt = get(hAxes, 'CurrentPoint');  x = pt(1,1);    y = pt(1,2);
	x_lim = get(hAxes,'xlim');      y_lim = get(hAxes,'ylim');
	if (x<x_lim(1)) || (x>x_lim(2)) || (y<y_lim(1)) || (y>y_lim(2));   return; end
	uirestore_fig(state);           % Restore the figure's initial state

% -----------------------------------------------------------------------------------------
function change_SymbPos(obj,eventdata,h)
% Change the Symbol position by asking it's coordinates

tag = get(h,'Tag');
if (strcmp(tag,'Pointpolyline') || strcmp(tag,'Maregraph'))
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
            xx = num2str(xx);    yy = num2str(yy);
    end
else
    xx = num2str(xx);    yy = num2str(yy);
end

prompt = {'Enter new lon (or x)' ,'Enter new lat (or y)'};
%resp  = inputdlg(prompt,'Move symbol',[1 30; 1 30],{num2str(xx) num2str(yy)});
resp  = inputdlg(prompt,'Move symbol',[1 30; 1 30],{xx yy});
if isempty(resp);    return;     end

val_x = test_dms(resp{1});           % See if coords were given in dd:mm or dd:mm:ss format
val_y = test_dms(resp{2});
x = 0;     y = 0;
for (k = 1:length(val_x)),  x = x + sign(str2double(val_x{1}))*abs(str2double(val_x{k})) / (60^(k-1));    end
for (k = 1:length(val_y)),  y = y + sign(str2double(val_y{1}))*abs(str2double(val_y{k})) / (60^(k-1));    end

if (is_single)      % Individual symbol
    xp = x;         yp = y;
else                % Picked symbol from a list
    xp(i) = x;      yp(i) = y;
end
set(h, 'XData', xp, 'YData', yp);

% -----------------------------------------------------------------------------------------
function remove_one_from_many(obj,eventdata,h)
%Delete one symbol that belongs to a class (in fact a vertex of a polyline)
pt = get(gca,'CurrentPoint');
xp = get(h,'XData');    yp = get(h,'YData');
l = length(xp);
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

% -----------------------------------------------------------------------------------------
function other_SymbSize(obj,eventdata,h)
prompt = {'Enter new size (pt)'};     dlg_title = 'Symbol Size';
num_lines= [1 30];
resp  = inputdlg(prompt,dlg_title,num_lines);
if isempty(resp);    return;     end
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
h=h(ishandle(h));

tag = get(h,'Tag');
if iscell(tag)          % When several symbol of class "tag" exists
    h_all = findobj(gca,'Tag',tag{1});
else                    % Only one symbol of class "tag" exists
    h_all = findobj(gca,'Tag',tag);
end
delete(h_all)

% -----------------------------------------------------------------------------------------
function remove_singleContour(obj,eventdata,h)
% Delete an individual contour and its eventual label(s)
labHand = getappdata(h,'LabelHands');
if (~isempty(labHand))
    try     delete(labHand);   end
end
delete(h)

% -----------------------------------------------------------------------------------------
function save_line(obj,eventdata,h)
% Save either individual as well as class lines. The latter uses the ">" symbol to separate segments
if (nargin == 3),   h = h(ishandle(h));     end
handles = guidata(gcbo);        % I hope I don't get into troubles because of this!
cd(handles.work_dir)
[FileName,PathName] = uiputfile({ ...
    '*.dat;*.DAT', 'Symbol file (*.dat,*.DAT)'; '*.*', 'All Files (*.*)'}, 'Select Symbol File name');
cd(handles.home_dir);       % allways go home to avoid troubles
if isequal(FileName,0);   return;     end
pause(0.01)
if (nargin == 2)    % Save only one line
    x = get(gco,'XData');    y = get(gco,'YData');
else                % Save a line class
    x = get(h,'XData');      y = get(h,'YData');
end

fname = [PathName FileName];
[PATH,FNAME,EXT] = fileparts([PathName FileName]);
if isempty(EXT),    fname = [PathName FNAME '.dat'];    end

fid = fopen(fname, 'w');
if (fid < 0),   errordlg(['Can''t open file:  ' fname],'Error');    return;     end
if (~iscell(x))
   	fprintf(fid,'%.5f\t%.5f\n',[x(:)'; y(:)']);
else
    for i=1:length(h)
        fprintf(fid,'%s\n','>');
        fprintf(fid,'%.5f\t%.5f\n',[x{i}(:)'; y{i}(:)']);
    end
end
fclose(fid);

% -----------------------------------------------------------------------------------------
function export_symbol(obj,eventdata,h, opt)
% If OPT is given than it must contain a Mx3 array with the x,y,z data to be saved

if (nargin == 3)
	h = h(ishandle(h));
	tag = get(h,'Tag');
	xx = get(h,'XData');    yy = get(h,'YData');
	if (length(xx) > 1 && ~strcmp(tag,'Pointpolyline') && ~strcmp(tag,'Maregraph'))     % (don't remember why)
        % Points and Maregraphs may be many but don't belong to a class
        msgbox('Only individual symbols may be exported and this one seams to belong to a class of symbols. Exiting','Warning')
        return
	end
    doSave_formated(xx, yy)
else
    errordlg('export_symbol called with a wrong number of arguments.','ERROR')
end

% -----------------------------------------------------------------------------------------
function save_formated(obj,eventdata, h, opt)
	% Save x,y[,z] vars into a file but taking into account the 'LabelFormatType'
	% If OPT is given than it must contain a Mx3 array with the x,y,z data to be saved

	if (nargin == 3)
		h = gco;
		xx = get(h,'XData');    yy = get(h,'YData');
        doSave_formated(xx, yy)
	elseif (nargin == 4)
        if (size(opt,2) ~= 3)
            errordlg('save_formated: variable must contain a Mx3 array.','ERROR')
            return
        end
        doSave_formated(opt(:,1), opt(:,2), opt(:,3))
	else
        errordlg('save_formated: called with a wrong number of arguments.','ERROR')
	end

% -----------------------------------------------------------------------------------------
function doSave_formated(xx, yy, opt_z)
	% Save x,y[,z] vars into a file but taking into account the 'LabelFormatType'
	% OPT_Z is what the name says, optional
	hFig = get(0,'CurrentFigure');
	handles = guidata(hFig);
	cd(handles.work_dir)
	[FileName,PathName] = uiputfile({ ...
        '*.dat;*.DAT', 'Symbol file (*.dat,*.DAT)'; '*.*', 'All Files (*.*)'}, 'Select Symbol File name');
	cd(handles.home_dir);       % allways come home to avoid troubles
	if isequal(FileName,0),   return;     end
	pause(0.01)

	[PATH,FNAME,EXT] = fileparts([PathName FileName]);
	if isempty(EXT),    f_name = [PathName FNAME '.dat'];
	else                f_name = [PathName FNAME EXT];       end

	% Save data with a format determined by axes format
	labelType = getappdata(handles.axes1,'LabelFormatType');             % find the axes label format
	if isempty(labelType),      labelType = ' ';        end     % untempered matlab axes labels
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
    else
        return
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
        msgbox(str,'Satation info')
	% elseif (strcmp(opt,'calendar'))
	%     date = clock;
	%     tim = datenum(date(1),date(2),1):1/24:datenum(date(1),date(2),31);
	%     out = t_xtide(pt(1,1),pt(1,2),tim,'format','times');
	end

% -----------------------------------------------------------------------------------------
function volcano_info(obj,eventdata,h,name,desc,dating)
	i = get(gco,'Userdata');
	msgbox( sprintf(['Volcano name: ' name{i} '\n' 'Volcano type:   ' desc{i} '\n' ...
            'Activity:      ' dating{i}] ),'Volcano info')

% -----------------------------------------------------------------------------------------
function ODP_info(obj,eventdata,h,leg,site,z,penetration)
	i = get(gco,'Userdata');
	tag = get(h,'Tag');     tag = tag{1};
	msgbox( sprintf([[tag ' Leg:    '] leg{i} '\n' [tag ' Site:    '] site{i} '\n' ...
			'Depth:      ' z{i} '\n' 'Hole penetration: ' penetration{i}] ),'ODP info')

% -----------------------------------------------------------------------------------------
function Isochrons_Info(obj,eventdata,data)
	i = get(gco,'Userdata');
	if (isstruct(i))    % This happens when h is ui_edit_polygon(ed)
		i = i.old_ud;
	end
	tag = data{i};
	msgbox( sprintf(tag),'This line info')

% -----------------------------------------------------------------------------------------
function gmtfile_Info(obj,eventdata,h,data)
str{1} = ['N_recs = ' num2str(data(1)) ', N_grav = ' num2str(data(2)) ', N_mag = ' num2str(data(3)) ...
    ', N_top = ' num2str(data(4))];
str{2} = ['E: = ' num2str(data(5)) '  W: = ' num2str(data(6))];
str{3} = ['S: = ' num2str(data(7)) '  N: = ' num2str(data(8))];
str{4} = ['Start day,month,year: = ' num2str(data(9)) '  ' num2str(data(10)) '  ' num2str(data(11))];
str{5} = ['End   day,month,year: = ' num2str(data(12)) '  ' num2str(data(13)) '  ' num2str(data(14))];
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
    case 'OTF',        txt_class = 'Oceanic Transform Fault';
    case 'OSF',        txt_class = 'Oceanic Spreadin Ridge';
    case 'CRB',        txt_class = 'Continental Rift Boundary';
    case 'CTF',        txt_class = 'Continental Transform Fault';
    case 'CCB',        txt_class = 'Continental Convergent Boundary';
    case 'OCB',        txt_class = 'Oceanic Convergent Boundary';
    case 'SUB',        txt_class = 'Subduction Zone';
end

if isempty(txt_id),     txt_id = char(data(i).pb_id);   end      % If id was not decoded, print id
if isempty(txt_class),  txt_class = char(data(i).class);   end   % Shouldn't happen, but just in case

msgbox( sprintf(['Plate pairs:           ' txt_id '\n' 'Boundary Type:    ' txt_class '\n' ...
        'Speed (mm/a):       ' sprintf('%g',data(i).vel) '\n' ...
        'Speed Azimuth:      ' sprintf('%g',data(i).azim_vel)] ),'Segment info')

% -----------------------------------------------------------------------------------------
function delete_obj(hTesoura)
	% hTesoura is the handle to the 'Tesoura' uitoggletool
	% Build the scisors pointer (this was done with the help of an image file)
	pointer = ones(16)*NaN;
	pointer(2,7) = 1;       pointer(2,11) = 1;      pointer(3,7) = 1;       pointer(3,11) = 1;
	pointer(4,7) = 1;       pointer(4,11) = 1;      pointer(5,7) = 1;       pointer(5,8) = 1;
	pointer(5,10) = 1;      pointer(5,11) = 1;      pointer(6,8) = 1;       pointer(6,10) = 1;
	pointer(7,8) = 1;       pointer(7,9) = 1;       pointer(7,10) = 1;      pointer(8,9) = 1;
	pointer(9,8) = 1;       pointer(9,9) = 1;       pointer(9,10) = 1;      pointer(10,8) = 1;
	pointer(10,10) = 1;     pointer(10,11) = 1;     pointer(10,12) = 1;     pointer(11,6) = 1;
	pointer(11,7) = 1;      pointer(11,8) = 1;      pointer(11,10) = 1;     pointer(11,13) = 1;
	pointer(12,5) = 1;      pointer(12,8) = 1;      pointer(12,10) = 1;     pointer(12,13) = 1;
	pointer(13,5) = 1;      pointer(13,8) = 1;      pointer(13,10) = 1;     pointer(13,13) = 1;
	pointer(14,5) = 1;      pointer(14,8) = 1;      pointer(14,11) = 1;     pointer(14,12) = 1;
	pointer(15,6:7) = 1;
    
    hFig = get(get(hTesoura,'Parent'),'Parent');
    state = uisuspend_fig(hFig);
    set(hFig,'Pointer','custom','PointerShapeCData',pointer,'PointerShapeHotSpot',[1 8],...
        'WindowButtonDownFcn',{@wbd_delObj,hFig,hTesoura,state})
    
function wbd_delObj(obj,event,hFig,hTesoura,state)
    stype = get(hFig,'selectiontype');
    if (stype(1) == 'a')                    % A right click ('alt'), end killing
        uirestore_fig(state)
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
function del_insideRect(obj,eventdata,h)
    % Delete all lines/patches/text objects that have at least one vertex inside the rectangle
    
    s = getappdata(h,'polygon_data');
    if (~isempty(s))            % If the rectangle is in edit mode, force it out of edit
        if strcmpi(s.controls,'on'),    ui_edit_polygon(h);     end
    end
    set(h, 'HandleVis','off')           % Make the rectangle handle invisible
    hAxes = get(h,'Parent');
    
    hLines = findobj(hAxes,'Type','line');     % Fish all objects of type line in Mirone figure
    hPatch = findobj(hAxes,'Type','patch');
    hText = findobj(hAxes,'Type','text');
    hLP = [hLines(:); hPatch(:)];
    rx = get(h,'XData');        ry = get(h,'YData');
    rx = [min(rx) max(rx)];     ry = [min(ry) max(ry)];
    found = false;
    for (i=1:numel(hLP))    % Loop over objects to find if any is on edit mode
        s = getappdata(hLP(i),'polygon_data');
        if (~isempty(s))
            if strcmpi(s.controls,'on')     % Object is in edit mode, so this
                ui_edit_polygon(hLP(i))     % call will force out of edit mode
                found = true;
            end
        end
    end
    if (found)      % We have to do it again because some line handles have meanwhile desapeared
        hLines = findobj(hAxes,'Type','line');
        hPatch = findobj(hAxes,'Type','patch');
        hLP = [hLines(:); hPatch(:)];
    end
    for (i=1:numel(hLP))        % Loop over objects to find out which cross the rectangle
        x = get(hLP(i),'XData');        y = get(hLP(i),'YData');
        if ( any( (x >= rx(1) & x <= rx(2)) & (y >= ry(1) & y <= ry(2)) ) )
            delete(hLP(i))
        end
    end

    found = false;
    for (i=1:numel(hText))      % Text objs are a bit different, so treat them separately
        pos = get(hText(i),'Position');
        if ( (pos(1) >= rx(1) && pos(1) <= rx(2)) && (pos(2) >= ry(1) && pos(2) <= ry(2)) )
            delete(hText(i))
            found = true;
        end
    end
    if (found),     refresh;    end     % Bloody text bug
    set(h, 'HandleVis','on')    % Make the rectangle handle findable again

% -------------------------------------------------------------------------------------------------------
function changeAxesLabels(opt)
% This function formats the axes labels strings using a geographical notation
hFig = get(0,'CurrentFigure');      hAxes = get(hFig,'CurrentAxes');
x_tick = getappdata(hAxes,'XTickOrig');
y_tick = getappdata(hAxes,'YTickOrig');
n_xtick = size(x_tick,1);                   n_ytick = size(y_tick,1);
switch opt
    case 'ToDegDec'
        % This is easy because original Labels where saved in appdata
        set(hAxes,'XTickLabel',getappdata(hAxes,'XTickOrig'));
        set(hAxes,'YTickLabel',getappdata(hAxes,'YTickOrig'))
        setappdata(hAxes,'LabelFormatType','DegDec')       % Save it so zoom can know the label type
    case 'ToDegMin'
        x_str = degree2dms(str2num( ddewhite(x_tick) ),'DDMM',0,'str');     % x_str is a structure with string fields
        y_str = degree2dms(str2num( ddewhite(y_tick) ),'DDMM',0,'str');
        str_x = [x_str.dd repmat(' ',n_xtick,1) x_str.mm];
        str_y = [y_str.dd repmat(' ',n_ytick,1) y_str.mm];
        set(hAxes,'XTickLabel',str_x);        set(hAxes,'YTickLabel',str_y)
        setappdata(hAxes,'LabelFormatType','DegMin')        % Save it so zoom can know the label type
    case 'ToDegMinDec'
        x_str = degree2dms(str2num( ddewhite(x_tick) ),'DDMM.x',2,'str');    % x_str is a structure with string fields
        y_str = degree2dms(str2num( ddewhite(y_tick) ),'DDMM.x',2,'str');
        str_x = [x_str.dd repmat(' ',n_xtick,1) x_str.mm];
        str_y = [y_str.dd repmat(' ',n_ytick,1) y_str.mm];
        set(hAxes,'XTickLabel',str_x);        set(hAxes,'YTickLabel',str_y)
        setappdata(hAxes,'LabelFormatType','DegMinDec')     % Save it so zoom can know the label type
    case 'ToDegMinSec'
        x_str = degree2dms(str2num( ddewhite(x_tick) ),'DDMMSS',0,'str');    % x_str is a structure with string fields
        y_str = degree2dms(str2num( ddewhite(y_tick) ),'DDMMSS',0,'str');
        str_x = [x_str.dd repmat(' ',n_xtick,1) x_str.mm repmat(' ',n_xtick,1) x_str.ss];
        str_y = [y_str.dd repmat(' ',n_ytick,1) y_str.mm repmat(' ',n_ytick,1) y_str.ss];
        set(hAxes,'XTickLabel',str_x);        set(hAxes,'YTickLabel',str_y)
        setappdata(hAxes,'LabelFormatType','DegMinSec')      % Save it so zoom can know the label type
    case 'ToDegMinSecDec'
        x_str = degree2dms(str2num( ddewhite(x_tick) ),'DDMMSS.x',1,'str');   % x_str is a structure with string fields
        y_str = degree2dms(str2num( ddewhite(y_tick) ),'DDMMSS.x',1,'str');
        str_x = [x_str.dd repmat(' ',n_xtick,1) x_str.mm repmat(' ',n_xtick,1) x_str.ss];
        str_y = [y_str.dd repmat(' ',n_ytick,1) y_str.mm repmat(' ',n_ytick,1) y_str.ss];
        set(hAxes,'XTickLabel',str_x);        set(hAxes,'YTickLabel',str_y)
        setappdata(hAxes,'LabelFormatType','DegMinSecDec')   % Save it so zoom can know the label type
end

% -----------------------------------------------------------------------------------------
function sout = ddewhite(s)
	%DDEWHITE Double dewhite. Strip both leading and trailing whitespace.
	%
	%   DDEWHITE(S) removes leading and trailing white space and any null characters
	%   from the string S.  A null character is one that has an absolute value of 0.
	
	%   Author:      Peter J. Acklam
	%   Time-stamp:  2002-03-03 13:45:06 +0100
	
	error(nargchk(1, 1, nargin));
	if ~ischar(s),   error('DDWHITE ERROR: Input must be a string (char array).');     end
	if isempty(s),   sout = s;      return;     end
	
	[r, c] = find(~isspace(s));
	if (size(s, 1) == 1),   sout = s(min(c):max(c));
	else                    sout = s(:,min(c):max(c));
	end

% --------------------------------------------------------------------
function draw_MagBarCode
handles = guidata(get(0,'CurrentFigure'));

% Define height of the bar code
tmax = 165;             % Max time (in Ma) represented
width = 8.0;            % Figure width in cm
height = 17;            % Figure height in cm
barHeight_1M = 0.75;    % Scale factor for the height a 1 Ma bar will show in the display
tscal = tmax * barHeight_1M;
dy = 1;             % dy is used to shift down all bar code. It represents 1 Ma in the time scale

% Create figure for bar code display
F = figure('Units','centimeters',...
    'Position',[1 1 [width height]],...
    'Toolbar','none',...
    'Menubar','none',...
    'Numbertitle','off',...
    'Name','Geomagnetic Bar Code',...
    'doublebuffer','on',...
    'Visible','of',...
    'renderer','Zbuffer',...    % Otherwise the patch command below would set it to OpenGL and
    'Color','k');               % despite what TMW says, in R13 it is awfully bugy.

% Create axis for bar code display
axes('Units','normalized',...
    'Position',[0 0 1 1],...
    'XLim',[0 width],...
    'YLim',[0 height],...
    'Color','w',...
    'XTick',[],...
    'YTick',[],...
    'YDir','reverse',...
    'DataAspectRatio',[1 1 1],...
    'DefaultPatchEdgecolor', 'none');

pos=[0.96 0 .04 1];
S=['set(gca,''ylim'',[' num2str(tscal-height) ' ' num2str(tscal) ']-get(gcbo,''value''))'];
uicontrol('style','slider','units','normalized','position',pos,...
    'callback',S,'min',0,'max',tscal-height,'Value',tscal-height);

%fid = fopen([handles.path_data 'BarCode_Cox_direct.dat'],'r');
fid = fopen([handles.path_data 'Cande_Kent_95.dat'],'r');
todos = fread(fid,'*char');     [chron age_start age_end age_txt] = strread(todos,'%s %f %f %s');
fclose(fid);    clear todos

y = barHeight_1M * ([age_start'; age_end'] + dy);
y = y(:);
x = [repmat(0.1,length(y),1); repmat(0.4,length(y),1)]*width; 
y = [y; y(end:-1:1)];
vert = [x y];

n_ages = length(age_start);
n2 = 2 * n_ages;
c1 = (1:n2-1)';     c3 = n2*2 - c1;
c2 = c3 + 1;        c4 = c1 + 1;
faces = [c1 c2 c3 c4];

cor = repmat([0 0 0; 1 1 1],n_ages-1,1);    cor = [cor; [0 0 0]];

patch('Faces',faces,'Vertices',vert,'FaceVertexCData',cor,'FaceColor','flat');

for (i=1:n_ages)
    if (~strcmp(age_txt(i),'a'))   % Plot anomaly names
        text(width*.42, barHeight_1M*(age_start(i)+dy),age_txt(i),'VerticalAlignment','cap')
    end
end

% Draw a vertical line for time ruler
line('XData',width*[.55 .55],'YData',barHeight_1M*([0+dy tmax+dy]),'LineWidth',2)  

for i=0:5:tmax
    line('XData',width*[.55 .58],'YData',barHeight_1M * ([i i]+dy),'LineWidth',.5)   % Time tick marks
    text(width*.6, barHeight_1M * (i + dy), sprintf('%g%s',i,' Ma'))
end
for i=0:tmax
    line('XData',width*[.55 .565],'YData',barHeight_1M * ([i i]+dy),'LineWidth',.5)   % 1 Ma ticks
end

% Draw two vertical lines for known geomagnetic periods
line('XData',width*[.015 .015],'YData',barHeight_1M*([0+dy 5.4+dy]),'LineWidth',1)  
line('XData',width*[.08 .08],'YData',barHeight_1M*([0+dy 5.4+dy]),'LineWidth',1)  
text(width*.04, barHeight_1M*(.35 + dy), 'Bru','rotation',90,'HorizontalAlignment','center')
line('XData',width*[.015 .08],'YData',barHeight_1M*([.73 .73]+dy),'LineWidth',.5)   % period separator
text(width*.04, barHeight_1M*(1.62 + dy), 'Mathu','rotation',90,'HorizontalAlignment','center')
line('XData',width*[.015 .08],'YData',barHeight_1M*([2.5 2.5]+dy),'LineWidth',.5)   % period separator
text(width*.04, barHeight_1M*(2.95 + dy), 'Gau','rotation',90,'HorizontalAlignment','center')
line('XData',width*[.015 .08],'YData',barHeight_1M*([3.4 3.4]+dy),'LineWidth',.5)   % period separator
text(width*.04, barHeight_1M*(4.4 + dy), 'Gilbert','rotation',90,'HorizontalAlignment','center')
line('XData',width*[.015 .08],'YData',barHeight_1M*([5.4 5.4]+dy),'LineWidth',.5)   % period separator

%--------------------------
% Draw two vertical lines for geological periods
line('XData',width*[.79 .79],'YData',barHeight_1M*([0+dy tmax+dy]),'LineWidth',1)  
line('XData',width*[.88 .88],'YData',barHeight_1M*([0+dy tmax+dy]),'LineWidth',1)  

% Plistocene
text(width*.84, barHeight_1M*(1 + dy), ' Plistocene','rotation',90,'HorizontalAlignment','center')
line('XData',width*[.79 .88],'YData',barHeight_1M*([2 2]+dy),'LineWidth',.5)   % period separator

% Pliocene
text(width*.84, barHeight_1M*(3.5 + dy), 'Pliocene','rotation',90,'HorizontalAlignment','center')
line('XData',width*[.79 .88],'YData',barHeight_1M*([5 5]+dy),'LineWidth',.5)   % period separator

% Miocene
text(width*.84, barHeight_1M*(14.75 + dy), 'Miocene','rotation',90,'HorizontalAlignment','center')
line('XData',width*[.79 .88],'YData',barHeight_1M*([24.5 24.5]+dy),'LineWidth',.5)   % period separator

% Oligocene
text(width*.84, barHeight_1M*(31.25 + dy), 'Oligocene','rotation',90,'HorizontalAlignment','center')
line('XData',width*[.79 .88],'YData',barHeight_1M*([38.0 38.0]+dy),'LineWidth',.5)   % period separator

% Eocene
text(width*.84, barHeight_1M*(46.5 + dy), 'Eocene','rotation',90,'HorizontalAlignment','center')
line('XData',width*[.79 .88],'YData',barHeight_1M*([55.0 55.0]+dy),'LineWidth',.5)   % period separator

% Paleocene
text(width*.84, barHeight_1M*(60.0 + dy), 'Eocene','rotation',90,'HorizontalAlignment','center')
line('XData',width*[.79 .88],'YData',barHeight_1M*([65.0 65.0]+dy),'LineWidth',.5)   % period separator

% Maastricthian
text(width*.84, barHeight_1M*(69.0 + dy), 'Maastricthian','rotation',90,'HorizontalAlignment','center')
line('XData',width*[.79 .88],'YData',barHeight_1M*([73.0 73.0]+dy),'LineWidth',.5)   % period separator

% Campanian
text(width*.84, barHeight_1M*(78.0 + dy), 'Campanian','rotation',90,'HorizontalAlignment','center')
line('XData',width*[.79 .88],'YData',barHeight_1M*([83.0 83.0]+dy),'LineWidth',.5)   % period separator

% Santonian
text(width*.84, barHeight_1M*(85.2 + dy), 'Santonian','rotation',90,'HorizontalAlignment','center')
line('XData',width*[.79 .88],'YData',barHeight_1M*([87.4 87.4]+dy),'LineWidth',.5)   % period separator

% Coniacian
text(width*.84, barHeight_1M*(87.95 + dy), 'Con','rotation',90,'HorizontalAlignment','center')
line('XData',width*[.79 .88],'YData',barHeight_1M*([88.5 88.5]+dy),'LineWidth',.5)   % period separator

% Turonian
text(width*.84, barHeight_1M*(89.75 + dy), 'Turonian','rotation',90,'HorizontalAlignment','center')
line('XData',width*[.79 .88],'YData',barHeight_1M*([91.0 91.0]+dy),'LineWidth',.5)   % period separator

% Cenomanian
text(width*.84, barHeight_1M*(94.25 + dy), 'Cenomanian','rotation',90,'HorizontalAlignment','center')
line('XData',width*[.79 .88],'YData',barHeight_1M*([97.5 97.5]+dy),'LineWidth',.5)   % period separator

% Albian
text(width*.84, barHeight_1M*(100.25 + dy), 'Albian','rotation',90,'HorizontalAlignment','center')
line('XData',width*[.79 .88],'YData',barHeight_1M*([103.0 103.0]+dy),'LineWidth',.5)   % period separator

% Aptian
text(width*.84, barHeight_1M*(110.0 + dy), 'Aptian','rotation',90,'HorizontalAlignment','center')
line('XData',width*[.79 .88],'YData',barHeight_1M*([119.0 119.0]+dy),'LineWidth',.5)   % period separator

% Barremian
text(width*.84, barHeight_1M*(122.0 + dy), 'Barremian','rotation',90,'HorizontalAlignment','center')
line('XData',width*[.79 .88],'YData',barHeight_1M*([125.0 125.0]+dy),'LineWidth',.5)   % period separator

% Hauterivian
text(width*.84, barHeight_1M*(128.0 + dy), 'Hauterivian','rotation',90,'HorizontalAlignment','center')
line('XData',width*[.79 .88],'YData',barHeight_1M*([131.0 131.0]+dy),'LineWidth',.5)   % period separator

% Valanginian
text(width*.84, barHeight_1M*(134.5 + dy), 'Valanginian','rotation',90,'HorizontalAlignment','center')
line('XData',width*[.79 .88],'YData',barHeight_1M*([138.0 138.0]+dy),'LineWidth',.5)   % period separator

% Berriasian
text(width*.84, barHeight_1M*(141.0 + dy), 'Berriasian','rotation',90,'HorizontalAlignment','center')
line('XData',width*[.79 .88],'YData',barHeight_1M*([144.0 144.0]+dy),'LineWidth',.5)   % period separator

% Tithonian
text(width*.84, barHeight_1M*(147.0 + dy), 'Tithonian','rotation',90,'HorizontalAlignment','center')
line('XData',width*[.79 .88],'YData',barHeight_1M*([150.0 150.0]+dy),'LineWidth',.5)   % period separator

% Kimmeridgian
text(width*.84, barHeight_1M*(153.0 + dy), 'Kimmeridgian','rotation',90,'HorizontalAlignment','center')
line('XData',width*[.79 .88],'YData',barHeight_1M*([156.0 156.0]+dy),'LineWidth',.5)   % period separator

% Oxfordian
text(width*.84, barHeight_1M*(159.5 + dy), 'Oxfordian','rotation',90,'HorizontalAlignment','center')
line('XData',width*[.79 .88],'YData',barHeight_1M*([163.0 163.0]+dy),'LineWidth',.5)   % period separator

set(F,'Visible','on')

% --------------------------------------------------------------------
function set_transparency(obj,eventdata)
% Sets the transparency of a patch object

h_patch = gco;
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
if ~(strcmp(r_mode,'auto'))
    set(handles.figure1,'RendererMode','auto')
end

set(handles.figure1,'doublebuffer','on')        % I may be wrong, but I think patches are full of bugs

% Define height of the bar code
width = 8.0;            % Figure width in cm
height = 2.0;           % Figure height in cm

% Create figure for transparency display
F = figure('Units','centimeters',...
    'Position',[1 10 [width height]],...
    'Toolbar','none', 'Menubar','none',...
    'Numbertitle','off',...
    'Name','Transparency',...
    'doublebuffer','on',...
    'RendererMode','auto',...
    'Visible','of',...
    'Color',[.75 .75 .75]);

% Create axis for transparency display
axes('Units','normalized',...
    'Position',[0 0 1 1],...
    'XLim',[0 width], 'YLim',[0 height],...
    'Color',[.75 .75 .75],...
    'XTick',[],'YTick',[],...
    'DataAspectRatio',[1 1 1],...
    'DefaultPatchEdgecolor', 'none');

T = uicontrol('style','text','string','Transparency  ',...
	'fontweight','bold','horizontalalignment','left',...
	'units','normalized','pos',[0.05  0.2  0.7  0.18],...
	'backgroundcolor',[.75 .75 .75]);

transp = get(h_patch,'FaceAlpha');     % Get the previous transparency value

pos=[0.05 0.5 .9 .2];
S = {@apply_transparency,T,h_patch,handles};
uicontrol('style','slider','units','normalized','position',pos,...
    'callback',S,'min',0,'max',1,'Value',transp);

set(F,'Visible','on')

% -----------------------------------------------------------------------------------------
function apply_transparency(obj,eventdata,T,h_patch,handles)
val = get(obj,'Value');
if (val > 0.01)
    set(h_patch,'FaceAlpha',val)
% else
%     h_all = findobj(handles.figure1,'Type','patch');
%     set_painters = 1;
%     for i = 1:length(h_all)
%         if (get(h_all(i),'FaceAlpha') > 0.04)
%             set_painters = 0;
%         end
%     end
%     if (set_painters)
%         set(handles.figure1,'Renderer','painters')
%     end
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
	item_lw = uimenu(cmenuHand, 'Label', 'Line Width', 'Separator','on');
	setLineWidth(item_lw,cb_LineWidth)
	item_ls = uimenu(cmenuHand, 'Label', 'Line Style');
	setLineStyle(item_ls,{cb_solid cb_dashed cb_dotted cb_dashdot})
	item7 = uimenu(cmenuHand, 'Label', 'Line Color');
	cb_color = uictx_color(h,'EdgeColor');      % there are 9 cb_color outputs
	setLineColor(item7,cb_color)

	set_stack_order(cmenuHand)      % Change order in the stackpot

	uimenu(item7, 'Label', 'None', 'Separator','on', 'Callback', 'set(gco, ''EdgeColor'', ''none'');refresh');
	item8 = uimenu(cmenuHand, 'Label','Fill Color', 'Separator','on');
	cb_color = uictx_color(h,'facecolor');      % there are 9 cb_color outputs
	setLineColor(item8,cb_color)
	uimenu(item8, 'Label', 'None', 'Separator','on', 'Callback', 'set(gco, ''FaceColor'', ''none'');refresh');
	uimenu(cmenuHand, 'Label', 'Transparency', 'Callback', @set_transparency);

% -----------------------------------------------------------------------------------------
function set_stack_order(cmenuHand)
	% Change order in the stackpot. cmenuHand is what it says. 
	item_order = uimenu(cmenuHand, 'Label', 'Order');
	uimenu(item_order, 'Label', 'Bring to Top', 'Callback','uistack_j(gco,''top'')');
	uimenu(item_order, 'Label', 'Send to Bottom', 'Callback','uistack_j(gco,''bottom'')');
	uimenu(item_order, 'Label', 'Move up', 'Callback','uistack_j(gco,''up'')');
	uimenu(item_order, 'Label', 'Move down', 'Callback','uistack_j(gco,''down'')');

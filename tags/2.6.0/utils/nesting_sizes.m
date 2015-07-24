function nesting_sizes(hand, opt)
% Compute the exact limits needed to create nested grids in COMCOT and NSWING
%
%	Initial data is fetch from the rectangles previously loaded and that are recognized
%	by having a 'NEST' Tag.
%	OPT may either be 'New' to create a new refined level rectangle or 
%	    'whatever' to print the info needed for COMCOT's .ctl and grid limits -R style

%	Copyright (c) 2004-2015 by J. Luis
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

	handles = guidata(hand);			% The Mirone handles
	hRects = findobj('tag','NEST');
	ud = get(hRects,'UserData');
	if (isa(ud, 'cell')),	ud = (cat(1,ud{:}));	end
	[ud,IX] = sort(ud);
	hRects = hRects(IX);			% Reorder them at the order they were plotted (bigger ones first)

	limits = cell(numel(hRects),1);
	for (k = 1:numel(hRects))
		limits{k} = [min(get(hRects(k),'XData')) max(get(hRects(k),'XData')) ...
			min(get(hRects(k),'YData')) max(get(hRects(k),'YData'))];
		ud = getappdata(hRects(k),'LineInfo');
		if (isempty(ud))		% First rectangle and when it was NOT imported via load_xyz
			limits{k}(5:7) = handles.head(7:9);
		else
			reg = '0';						% Default is grid registration
			[x_inc,r] = strtok(ud);
			y_inc = x_inc;					% Start by assuming this
			if (~isempty(r))
				if (strcmp(r,'0') || strcmp(r,'1'))		% Have X_INC and REG info (Y_INC = X_INC)
					reg = r;
					y_inc = x_inc;
				else
					[y_inc,r] = strtok(r);
					if (~isempty(r) && (strcmp(r,'0') || strcmp(r,'1')))
						reg = r;
					end
				end
			end
			limits{k}(5) = str2double(reg);
			limits{k}(6) = str2double(x_inc);
			limits{k}(7) = str2double(y_inc);
		end
	end

	if (nargin == 2)
		if (strcmp(opt,'New'))
			hNew = make_new_nested(handles.axes1, hRects, limits);
			if (isempty(hNew)),		return,		end
			nesting_sizes(hNew)		% Hate recursivity but it makes life easier here
		else
			hLine = hRects(hRects == gco);
			info = getappdata(hLine, 'info');
			opt_R = sprintf('-R%.12g/%.12g/%.12g/%.12g -I%.12g/%.12g\n', info.head(1:4), info.head(6:7));
			if (handles.geog)
				lims_dd = fix(info.head(1:4));
				lims_dec = abs(info.head(1:4) - lims_dd);
				lims_mm  = fix(lims_dec * 60);
				lims_dec = lims_dec * 60 - lims_mm;
				lims_ss  = lims_dec * 60;		% ss.xxx
				opt_Rddmmss = sprintf('-R%d:%02d:%.8g/%d:%02d:%.8g/%d:%02d:%.8g/%d:%02d:%.8g\n', ...
					lims_dd(1), lims_mm(1), lims_ss(1), lims_dd(2), lims_mm(2), lims_ss(2), ...
					lims_dd(3), lims_mm(3), lims_ss(3), lims_dd(4), lims_mm(4), lims_ss(4));
			else
				opt_Rddmmss = [];
			end
			inds = sprintf('x_start = %d\nx_end = %d\ny_start = %d\ny_end = %d', ...
				info.idx_min, info.idx_max, info.idy_min, info.idy_max);
			message_win('create',[opt_R opt_Rddmmss inds],'edit','yes')
		end
		return
	end

	resize2nesting_size(handles, hRects, limits)

% --------------------------------------------------------------------------
function resize2nesting_size(handles, hRects, limits)
% Resize all inner rectangles so they obey to the nesting size rules.
% Also store nesting info in appdata, which is available for "Show nesting info" 

	for (k = 1:numel(hRects))
		if (k == 1 && ~handles.validGrid)		% Parent rectangle in an empty region
			continue
		elseif (k == 1 && handles.validGrid)
			parent_lims = [handles.head(1:4) handles.head(7:9)];
		else
			parent_lims = head;					% "head" was computed in the previous iteration
		end
		nx = ((parent_lims(2) - parent_lims(1))/ parent_lims(6))+1;
		ny = ((parent_lims(4) - parent_lims(3))/ parent_lims(7))+1;

		%get parent grid coord vectors
		X = linspace(parent_lims(1),parent_lims(2),round(nx));
		Y = linspace(parent_lims(3),parent_lims(4),round(ny));

		%look for closest element
		xmin = find_nearest(X, parent_lims(6), limits{k}(6), limits{k}(1));
		xmax = find_nearest(X, parent_lims(6), limits{k}(6), limits{k}(2));
		ymin = find_nearest(Y, parent_lims(7), limits{k}(7), limits{k}(3));
		ymax = find_nearest(Y, parent_lims(7), limits{k}(7), limits{k}(4));

		this_xmin = xmin.val - (parent_lims(6)/2) + (limits{k}(6)/2);
		this_xmax = xmax.val + (parent_lims(6)/2) - (limits{k}(6)/2);
		this_ymin = ymin.val - (parent_lims(7)/2) + (limits{k}(7)/2);
		this_ymax = ymax.val + (parent_lims(7)/2) - (limits{k}(7)/2);
		head = [this_xmin this_xmax this_ymin this_ymax limits{k}(5:7)];

		resp.head = head;
		resp.idx_min = xmin.idx;	resp.idx_max = xmax.idx;
		resp.idy_min = ymin.idx;	resp.idy_max = ymax.idx;

		set(hRects(k),'XData', [this_xmin this_xmin this_xmax this_xmax this_xmin], ...
		'YData', [this_ymin this_ymax this_ymax this_ymin this_ymin])
		setappdata(hRects(k),'info',resp)
	end

% --------------------------------------------------------------------------
function resp = find_nearest(X, dxParent, dxChild, pt)
% Find the element on the Parent grid coord vector (& index) that is closest to PT

	difa = (X - pt);
	ind = find(difa >= 0);		% Find index of closest point on coord vector that is on the RIGHT of 'pt'

	if ( difa(ind(1)) <= abs(difa(ind(1)-1)) )	% Point is closer to X(ind(1))
		resp.val = X(ind(1));
		resp.idx = ind(1);
	else
		resp.val = X(ind(1)-1);
		resp.idx = ind(1)-1;
	end

% --------------------------------------------------------------------------
function hNew = make_new_nested(hAx, hRects, limits)
% Create a new rctangle with half the width of its parent.

	resp = fix(abs(str2double(inputdlg({'Enter refinement factor'},'Refinement factor',[1 30],{'5'}))));
	if (isempty(resp) || isnan(resp) || resp == 0),		hNew = [];	return,		end
	x = limits{end}(1:2);		y = limits{end}(3:4);
	x = x + [diff(x) -diff(x)] / 4;
	y = y + [diff(y) -diff(y)] / 4;
	x_inc = limits{end}(6) / resp;
	y_inc = limits{end}(7) / resp;
	
	cor = get(hRects(end),'Color');
	lThick = get(hRects(end),'Linewidth');
	hNew = line('XData',[x(1) x(1) x(2) x(2) x(1)],'YData',[y(1) y(2) y(2) y(1) y(1)],'Parent',hAx,'Linewidth',lThick,...
				'Color',cor,'Userdata',numel(limits)+1);

	setappdata(hNew,'LineInfo',[sprintf('%.16g %.16g',x_inc, y_inc), limits{end}(5)]);
	draw_funs([],'set_recTsu_uicontext', hNew)

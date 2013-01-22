function nesting_sizes(hand, opt)
% Compute the exact limits needed to create nested grids in COMCOT and TINTOL
%
%	Initial data is fetch from the rectangles previously loaded and that are recognized
%	by having a 'NEST' Tag.
%	OPT may either be 'New' to create a new refined level rectangle or 
%	    'whatever' to print the info needed for COMCOT's .ctl and grid limits -R style

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

% $Id: $

	handles = guidata(hand);			% The Mirone handles
	hRects = findobj('tag','NEST');
	ud = get(hRects,'UserData');
	if (isa(ud, 'cell'))
		ud = (cat(1,ud{:}));
	end
	hRects = hRects(ud);				% Reorder them at the order they were plotted (bigger ones first)

	for (k = 1:numel(hRects))
		limits{k} = [min(get(hRects(k),'XData')) max(get(hRects(k),'XData')) ...
			min(get(hRects(k),'YData')) max(get(hRects(k),'YData'))];
		ud = getappdata(hRects(k),'LineInfo');
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
		limits{k}(5) = str2double(x_inc);
		limits{k}(6) = str2double(y_inc);
		limits{k}(7) = str2double(reg);
	end

	if (nargin == 2)
		if(strcmp(opt,'New'))
			make_new_nested(handles.axes1, hRects, limits)
		else
			ind = find(hRects == gco);
			hLine = hRects(ind);
			info = getappdata(hLine, 'info');
			opt_R = sprintf('-R%.12g/%.12g/%.12g/%.12g\n', info.head(1:4));
			inds = sprintf('x_start = %d\nx_end = %d\ny_start = %d\ny_end = %d', ...
				info.idx_min, info.idx_max, info.idy_min, info.idy_max);
			message_win('create',[opt_R inds],'edit','yes')
		end
		return
	end

	% ..
	for (k = 2:numel(hRects))
	    nx = ((limits{k-1}(2) - limits{k-1}(1))/ limits{k-1}(5))+1;
		ny = ((limits{k-1}(4) - limits{k-1}(3))/ limits{k-1}(6))+1;

		%get parent grid coord vectors
		X = linspace(limits{k-1}(1),limits{k-1}(2),nx);
		Y = linspace(limits{k-1}(3),limits{k-1}(4),ny);

		%look for closest element
		xmin = find_nearest(X, limits{k-1}(5), limits{k}(5), limits{k}(1));
		xmax = find_nearest(X, limits{k-1}(5), limits{k}(5), limits{k}(2));
		ymin = find_nearest(Y, limits{k-1}(6), limits{k}(6), limits{k}(3));
		ymax = find_nearest(Y, limits{k-1}(6), limits{k}(6), limits{k}(4));

		this_xmin = xmin.val - (limits{k-1}(5)/2) + (limits{k}(5)/2);
		this_xmax = xmax.val + (limits{k-1}(5)/2) - (limits{k}(5)/2);
		this_ymin = ymin.val - (limits{k-1}(6)/2) + (limits{k}(6)/2);
		this_ymax = ymax.val + (limits{k-1}(6)/2) - (limits{k}(6)/2);
		h = [this_xmin this_xmax this_ymin this_ymax limits{k}(5:7)];

		resp.head = h;
		resp.idx_min = xmin.idx;
		resp.idx_max = xmax.idx;
		resp.idy_min = ymin.idx;
		resp.idy_max = ymax.idx;

		set(hRects(k),'XData', [this_xmin this_xmin this_xmax this_xmax this_xmin])
		set(hRects(k),'YData', [this_ymin this_ymax this_ymax this_ymin this_ymin])
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
function make_new_nested(hAx, hRects, limits)
% Create a new rctangle with half the width of its parent.
	
	resp = fix(abs(str2double(inputdlg({'Enter refinement factor'},'Refinement factor',[1 30],{'5'}))));
	if (isnan(resp) || resp == 0),		return,		end
	x = limits{end}(1:2);		y = limits{end}(3:4);
	x = x + [diff(x) -diff(x)] / 4;
	y = y + [diff(y) -diff(y)] / 4;
	hLine = copyobj(hRects(end), hAx);
	set(hLine, 'XData',[x(1) x(1) x(2) x(2) x(1)], 'YData',[y(1) y(2) y(2) y(1) y(1)],'Userdata',numel(limits)+1)
	setappdata(hLine,'LineInfo',[sprintf('%d',resp), limits{end}(7)]);
	h = findobj(get(hLine,'uicontextmenu'),'label','Save line');
	if (~isempty(h))        % Replace the old line handle in the 'Save line' Callback by the just created one
		hFun = get(h,'Call');
		hFun{2} = hLine;
		set(h,'Call',hFun)
	end
	if (isappdata(hLine,'polygon_data'))
		rmappdata(hLine,'polygon_data')		% Remove the parent's ui_edit_polygon appdata
	end

	ui_edit_polygon(hLine)
	
function [xx, yy, zz] = grid_profiler(hFig, xp, yp, point_int, do_dynamic)
% Compute profiles on the grid displayed by the Mirone figure.
%
% POINT_INT  true|false		If true interpolate only on [xp, yp] vertex
% DO_DYNAMIC true|false 	If true do dynamic profiling

	handles = guidata(hFig);
	[X,Y,Z,head] = load_grd(handles,'silent');

	n_nodes = numel(xp);
	if (~point_int)         % Profile interp
		dx = X(2) - X(1);   dy = Y(2) - Y(1);
		dx = handles.head(8);		dy = handles.head(9);
		% Here I don't realy know what is the good increment for interpolation, so I'll
		% interpolate at the half grid spacing for each dimension.
		% Construct the vectors with the points where to interpolate the profile
		xx = [];	yy = [];
		for (i = 1:n_nodes-1)
			n_int = round( max( abs(xp(i+1)-xp(i))/(dx/2), abs(yp(i+1)-yp(i))/(dy/2) ) );			% find ...
			xx = [xx linspace(xp(i),xp(i+1),n_int)];     yy = [yy linspace(yp(i),yp(i+1),n_int)];	% at nodes, values are repeated
		end
	else					% Interpolation at line vertex
		xx = xp;    yy = yp;
	end

	% Interpolate
	if ~isempty(getappdata(handles.figure1,'dem_x'))		% Grid is in memory
        if (~getappdata(handles.figure1,'PixelMode'))		% Interpolation mode
			zz = grdtrack_m(Z,head,[xx' yy'],'-Z')';		% It uses less memory (and non doubles)
			%zz = interp2(X,Y,double(Z),xx,yy,'*cubic');
		else												% NEARNEIGBOR mode
			[rows,cols] = size(Z);
			rp = aux_funs('axes2pix',rows, get(handles.hImg,'YData'),yy);
			cp = aux_funs('axes2pix',cols, get(handles.hImg,'XData'),xx);
			r = min(rows, max(1, round(rp)));	c = min(cols, max(1, round(cp)));
			rc = (c - 1) * rows + r;
			zz = double(Z(rc));
        end
	else								% grid was loaded here (big according to preferences), so interp linearly
		zz = bi_linear(X,Y,Z,xx,yy);
	end

	if (do_dynamic)

		hDynProfAx = getappdata(handles.axes1,'dynProfile');		% handle of the dynamic profile Axes

		if (isempty(hDynProfAx))			% First time
			% Create a new axes on the Mirone figure
			unit = get(handles.axes1,'Units');		set(handles.axes1,'Units', 'pixels');
			axSize = get(handles.axes1,'Pos');		set(handles.axes1,'Units', unit);
			hAxes = axes('Parent',hFig,'Units','pixels','Position',[axSize(1) (axSize(2)+axSize(4)/2) axSize(3) min(axSize(4),200)], ...
				'Tag','axDynProf','Visible','off');
			r = dist_along_profile(xx, yy);
			hLine = line('Parent',hAxes,'XData',r,'YData',zz,'Color',handles.DefLineColor,'LineWidth',2);
			set(hAxes,'xlim', [min(r) max(r)])
			setappdata(handles.axes1, 'dynProfile', hAxes)
			setappdata(hAxes,'theLine',hLine)
		else
			r = dist_along_profile(xx, yy);
			hLine = getappdata(hDynProfAx,'theLine');
			set(hLine,'XData',r,'YData',zz,'UserData',[xx(:) yy(:)])
			set(hDynProfAx,'xlim', [min(r) max(r)])
		end
	end

% -------------------------------------------------------------------------------------
function rd = dist_along_profile(x, y)
    xd = diff(x);		yd = diff(y);
    tmp = sqrt(xd.*xd + yd.*yd);
	rd = [0; cumsum(tmp(:))];

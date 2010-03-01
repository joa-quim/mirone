function [xx, yy, zz] = grid_profiler(hFig, xp, yp, point_int, do_dynamic, do_stack)
% Compute profiles on the grid displayed by the Mirone figure.
%
% POINT_INT  true|false		If true interpolate only on [xp, yp] vertex
% DO_DYNAMIC true|false 	If true do dynamic profiling

	handles = guidata(hFig);
	[X,Y,Z,head] = load_grd(handles,'silent');
	if (getappdata(handles.figure1,'PixelMode')),	do_stack = false;	end		% Case no programed
	
	if (nargin == 5),	do_stack = false;	end
	if (~point_int)         % Profile interp
		[xx, yy] = make_track(xp, yp, handles.head(8), handles.head(9), do_stack);
	else					% Interpolation at line vertex
		xx = xp;    yy = yp;
	end

	% Interpolate
	if ~isempty(getappdata(handles.figure1,'dem_x'))		% Grid is in memory
        if (~getappdata(handles.figure1,'PixelMode'))		% Interpolation mode
			if ((size(xx,1) == 1) || (size(xx,2) == 1))		% Typical, one tack interpolation only
				zz = grdtrack_m(Z,head,[xx' yy'],'-Z')';
			else							% Multi-track interpolation and stacking
				zz = grdtrack_m(Z,head,[xx(1,:)' yy(1,:)'],'-Z')';
				for (k = 2:size(xx,1))
					zz = zz + grdtrack_m(Z,head,[xx(k,:)' yy(k,:)'],'-Z')';
				end
				zz = zz / size(xx,1);
			end
			xx = xx(round(size(xx,1)/2),:);		yy = yy(round(size(xx,1)/2),:);		% Mid track coordinates
		else												% NEARNEIGBOR mode
			[rows,cols] = size(Z);
			rp = aux_funs('axes2pix',rows, get(handles.hImg,'YData'),yy);
			cp = aux_funs('axes2pix',cols, get(handles.hImg,'XData'),xx);
			r = min(rows, max(1, round(rp)));	c = min(cols, max(1, round(cp)));
			rc = (c - 1) * rows + r;
			zz = double(Z(rc));
        end
	elseif (isempty(Z) && ndims(get(handles.hImg,'CData')) == 2)	% Gray image inperp (linear)
		Z = get(handles.hImg,'CData');
		img_lims = getappdata(handles.axes1,'ThisImageLims');	% Get limits and correct them for the pix reg problem
		x_inc = (img_lims(2)-img_lims(1)) / size(Z,2);		y_inc = (img_lims(4)-img_lims(3)) / size(Z,1);
		img_lims = img_lims + [x_inc -x_inc y_inc -y_inc]/2;	% Remember that the Image is ALWAYS pix reg
		X = linspace(img_lims(1),img_lims(2),size(Z,2));	Y = linspace(img_lims(3),img_lims(4),size(Z,1));
		zz = bi_linear(X,Y,Z,xx,yy);
		if ((size(zz,1) ~= 1) && (size(zz,2) ~= 1))		% A thickned line
			z = 0;
			for (k = 1:size(zz,1))
				z = z + zz(k,:);
			end
			zz = z / size(zz,1);
			xx = xx(round(size(xx,1)/2),:);		yy = yy(round(size(xx,1)/2),:);		% Mid track coordinates
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
			axSize = [axSize(1) (axSize(2)+axSize(4)/2) axSize(3) min(axSize(4),200)];
			hAxes = axes('Parent',hFig,'Units','pixels','Position',axSize, ...
				'Tag','axDynProf','Visible','on', 'Color','none');
			r = dist_along_profile(xx, yy);
			hLine = line('Parent',hAxes,'XData',r,'YData',zz,'Color',handles.DefLineColor,'LineWidth',2);

			% Estimate a right shift so that it can accommodate the Ylabels
			ylab = get((hAxes),'YTickLabel');
			nChars = max(numel(ylab(1,:)), numel(ylab(2,:)));
			dx = nChars * 7;
			axSize(1) = axSize(1) + dx;
			axSize(3) = axSize(3) - dx;
			set(hAxes,'Position',axSize)

			r_max =  max(r);
			if (r_max == 0),	r_max = 1;	end		% To not set a [0 0] xlim
			set(hAxes,'xlim', [0 r_max])
			setappdata(handles.axes1, 'dynProfile', hAxes)
			setappdata(hAxes,'theLine',hLine)
		else
			r = dist_along_profile(xx, yy);
			hLine = getappdata(hDynProfAx,'theLine');
			set(hLine,'XData',r,'YData',zz,'UserData',[xx(:) yy(:)])
			r_max =  max(r);
			if (r_max == 0),	r_max = 1;	end		% To not set a [0 0] xlim
			set(hDynProfAx,'xlim', [0 r_max])
		end
	end

% -------------------------------------------------------------------------------------
function rd = dist_along_profile(x, y)
    xd = diff(x);		yd = diff(y);
    tmp = sqrt(xd.*xd + yd.*yd);
	rd = [0; cumsum(tmp(:))];

% -------------------------------------------------------------------------------------
function [xx, yy] = make_track(x, y, dx, dy, do_stack)
% Here I don't realy know what is the good increment for interpolation, so I'll
% interpolate at the half grid spacing for each dimension.
% Construct the vectors with the points where to interpolate the profile

	if (~do_stack)				% Simple, most common and original case (one track only interpolation)
		xx = [];	yy = [];
		n_vertex = numel(x);
		for (i = 1:n_vertex-1)
			n_int = round( max( abs(x(i+1)-x(i))/(dx/2), abs(y(i+1)-y(i))/(dy/2) ) );			% find ...
			xx = [xx linspace(x(i),x(i+1),n_int)];     yy = [yy linspace(y(i),y(i+1),n_int)];	% at nodes, values are repeated
		end
	else
		ud = get(do_stack, 'UserData');		% See end of case 'thicken' in line_operations() for ud contents
		[n_vertex, n_lines] = size(ud{1});
		if ( isequal(x(:), ud{4}(:,1)) && isequal(y(:), ud{4}(:,2)) )	% Stack track line was not edited
			for (k = 1:n_lines)
				x = ud{1}(:,k);			y = ud{2}(:,k);
				x0 = [];	y0 = [];
				for (i = 1:n_vertex-1)
					n_int = round( max( abs(x(i+1)-x(i))/(dx/2), abs(y(i+1)-y(i))/(dy/2) ) );
					x0 = [x0 linspace(x(i),x(i+1),n_int)];     y0 = [y0 linspace(y(i),y(i+1),n_int)];
				end
				xx(k,:) = x0;			yy(k,:) = y0;
				if (k == 1)		% Now that we know the final array size, we can pre-allocate
					xx = zeros(n_lines,numel(x0));				yy = zeros(n_lines,numel(x0));
					xx(k,:) = x0;		yy(k,:) = y0;
				end
			end

		else								% Shit. Line was edited. Need to recompute the new tracks in it

			% Compute line angle 
			if (ud{5})				% Is geog?
				[dumb, az] = vdist(y(1:end-1),x(1:end-1), y(2:end),x(2:end));
				az = 90 - az;				% Make it trigonometric
				co = cos(az * pi / 180);	si = sin(az * pi / 180);
			else
				dx = diff(x);				dy = diff(y);
				hy = (dx.^2 + dy.^2)^.5;
				co =  dx ./ hy;				si =  dy ./ hy;			% calculate the cosine and sine
			end
			
			% Compute the new lines stack collection
			dl = ud{3} / n_lines;
			xL = zeros(n_vertex, n_lines+1);		yL = zeros(n_vertex, n_lines+1);
			if (size(x,1) == 1),		x = x(:);	y = y(:);		end
			for (k = 1:n_lines+1)
				th = ud{3}/2 - (k-1) * dl;
				foo = [co -si; si  co] * [0; th];
				xL(:,k) = x + repmat(foo(1), n_vertex, 1);		% One line per column
				yL(:,k) = y + repmat(foo(2), n_vertex, 1);
			end

			% And finaly interpoate them at inc/2
			for (k = 1:n_lines)
				x = xL(:,k);		y = yL(:,k);
				x0 = [];			y0 = [];
				for (i = 1:n_vertex-1)
					n_int = round( max( abs(x(i+1)-x(i))/(dx/2), abs(y(i+1)-y(i))/(dy/2) ) );
					x0 = [x0 linspace(x(i),x(i+1),n_int)];     y0 = [y0 linspace(y(i),y(i+1),n_int)];
				end
				xx(k,:) = x0;		yy(k,:) = y0;
				if (k == 1)		% Now that we know the final array size, we can pre-allocate
					xx = zeros(n_lines,numel(x0));				yy = zeros(n_lines,numel(x0));
					xx(k,:) = x0;	yy(k,:) = y0;
				end
			end
			
		end			% End of, stack track line edited?

	end

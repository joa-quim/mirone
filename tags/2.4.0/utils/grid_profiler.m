function [xx, yy, zz] = grid_profiler(hFig, xp, yp, point_int, do_dynamic, do_stack)
% Compute profiles on the grid displayed by the Mirone figure.
%
% POINT_INT  true|false		If true interpolate only on [xp, yp] vertex
% DO_DYNAMIC true|false 	If true do dynamic profiling
%
% When profiling an RGB image the output variable ZZ is a cell array.

%	Copyright (c) 2004-2012 by J. Luis
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

	handles = guidata(hFig);
	if (nargin == 2 && ishandle(xp))	% Normally, a call from "Radial Average"
		ud = getappdata(xp, 'donut');
		try						% Use a try to find if the handle is from a donutified circle
			rad_out = ud(3);	rad_in = ud(4);
			% ###################### DO THE RADIAL AVERAGE AND LEAVE $$$$$$$$$$$$$$$
			profile = do_radialAverage(hFig, ud(1), ud(2), rad_in, rad_out, ud(5));
			[pato,name,ext] = fileparts(get(handles.figure1,'Name'));
			ext = strtok(ext);		% Remove the "@ ??%" part
			%ecran(handles, profile(:,1),profile(:,2),['Radial average from ' name ext])
			x1 = profile(:,1);		y1 = profile(:,2)+profile(:,3);		y2 = profile(:,2)-profile(:,3);
			ecran(handles, [x1;x1(end:-1:1);NaN; x1],[y1;y2(end:-1:1);NaN;profile(:,2)],['Radial average (+/- STD) ' name ext])
			return
			% ####################################$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
		catch
			if (strcmp(get(xp,'Tag'), 'Donut'))
				errordlg(['Radial average trouble: ' lasterr],'Error'),		return
			end
		end
		% Reach here when not dealing with a donutified circle. So keep going.
		yp = get(xp, 'YData');	xp = get(xp, 'XData');
		point_int = false;		do_dynamic = false;		do_stack = false;
	elseif (nargin == 5)
		do_stack = false;
	end

	[X,Y,Z,head] = load_grd(handles,'silent');
	if (getappdata(handles.figure1,'PixelMode')),	do_stack = false;	end		% Case not programed
	
	if (~point_int)         % Profile interp
		[xx, yy] = make_track(xp, yp, handles.head(8), handles.head(9), do_stack);
	else					% Interpolation at line vertex (cannot be stacked)
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
		else												% NEARNEIGBOR mode (non stackable)
			[rows,cols] = size(Z);
			rp = getPixel_coords(rows, get(handles.hImg,'YData'),yy);
			cp = getPixel_coords(cols, get(handles.hImg,'XData'),xx);
			r = min(rows, max(1, round(rp)));	c = min(cols, max(1, round(cp)));
			rc = (c - 1) * rows + r;
			zz = double(Z(rc));
        end
	elseif (~handles.validGrid)								% An image, interp linear
		Z = get(handles.hImg,'CData');
		img_lims = getappdata(handles.axes1,'ThisImageLims');	% Get limits and correct them for the pix reg problem
		x_inc = (img_lims(2)-img_lims(1)) / size(Z,2);		y_inc = (img_lims(4)-img_lims(3)) / size(Z,1);
		img_lims = img_lims + [x_inc -x_inc y_inc -y_inc]/2;	% Remember that the Image is ALWAYS pix reg
		X = linspace(img_lims(1),img_lims(2),size(Z,2));	Y = linspace(img_lims(3),img_lims(4),size(Z,1));

		if (ndims(Z) == 2)
			zz = bi_linear(X,Y,Z,xx,yy);
			if ((size(zz,1) ~= 1) && (size(zz,2) ~= 1))		% A thickned line
				zz = stack_it(zz);
			end
		else
			zz{1} = bi_linear(X,Y,Z(:,:,1),xx,yy);
			zz{2} = bi_linear(X,Y,Z(:,:,2),xx,yy);
			zz{3} = bi_linear(X,Y,Z(:,:,3),xx,yy);
			if ((size(zz{1},1) ~= 1) && (size(zz{1},2) ~= 1))	% A thickned line
				zz{1} = stack_it(zz{1});
				zz{2} = stack_it(zz{2});
				zz{3} = stack_it(zz{3});
			end
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
function z = stack_it(z)
% Stack a set of profiles contained in the matriz z = [n_profiles x n_pts]
	zz = 0;
	for (k = 1:size(z,1))
		zz = zz + z(k,:);
	end
	z = zz / size(z,1);

% -------------------------------------------------------------------------------------
function profile = do_radialAverage(hFig, clon, clat, rad_in, rad_out, geog)
% ...
	handles = guidata(hFig);
	[X,Y,Z,head] = load_grd(handles,'silent');
	dr = (head(8) + head(9)) / 2;		% Don't really know what's the best so I'll use the plain mean
	profile = zeros( round((rad_out - rad_in)/dr)+2,3);		% pre-allocate for the result
	done = false;		k = 1;
	rad = rad_in;		% Start with inner radius
	while (~done)
		% Aproximate number of 'dr' chunks in this perimeter = round(2*pi*rad*D2R*R / (dr*D2R*R))
		n_pts = round(2*pi*rad / dr);		% n points for this circle. Variable with perimeter
		if (rad > 0)
			if (geog)
				[latc, lonc] = circ_geo(clat, clon, rad, [], n_pts);
			else	% Cartesian
				xx = linspace(-pi,pi,n_pts);	yy = sin(xx);		xx = cos(xx);
				lonc = clon + rad * xx;			latc = clat + rad * yy;
			end
		else
			latc = clat;	lonc = clon;
			if (handles.geog == 2),		lonc = lonc - 360;	end		% To be consistent with above cases.
			n_pts = 1;
		end
		if (handles.geog == 2),		lonc = lonc + 360;		end		% Longitudes in the [0 360] interval
		zz = grdtrack_m(Z,head,[lonc(:) latc(:)],'-Z')';
		profile(k,1) = rad;
		profile(k,2) = sum(zz) / n_pts;
		profile(k,3) = std(zz);
		rad = rad + dr;
		done = (rad > rad_out);
		if (done && (rad - rad_out) < (dr / 2))		% Condition to also do the outer circle radius
			rad = rad_out;		done = false;
		end
		k = k + 1;
	end
	profile(k-1:end,:) = [];			% Removed unused

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

% -------------------------------------------------------------------------------------
function pix_coords = getPixel_coords(img_length, XData, axes_coord)
% Convert coordinates from axes (real coords) to image (pixel) coordinates.
% IMG_LENGTH is the image width (n_columns)
% XDATA is the image's [x_min x_max] in axes coordinates
% AXES_COORD is the (x,y) coordinate of the point(s) to be converted

	slope = (img_length - 1) / (XData(end) - XData(1));
	if ((XData(1) == 1) && (slope == 1))
		pix_coords = axes_coord;
	else
		pix_coords = slope * (axes_coord - XData(1)) + 1;
	end

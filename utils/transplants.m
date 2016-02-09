function varargout = transplants(hLine, tipo, handles)
% Replace a zone in the host grid/image by the contents of an external grid/image
%
% In the IMPLANTGRID mode this function:
% Read an external grid, that can have most of the recognized formats, and "inserts" it
% on the mother grid. Grid resolutions do not need to be equal as over the modified zone
% the grid is reinterpolated with gmtmbgrid. A padding zone of 1 cell is left between
% the two datasets so that the interpolation can be made minimally smooth.
% The to-be-inserted grid may have NaNs but only the outside ones will be acknowledged.
% That is, inner holes are ignored.
%
% At 9-Oct-2015 I changed a bit the above behavior in result of a situation wehere roipoly_j
% was found to (sometimes?) fail when both grid limits and rectangle were equal (one line
% of zeros at Left and Top). So now, if implanting grid contains host grid, has no NaNs in
% common region and is of cruder resolution than we use grdsample to obtain the wished result.
% If resolution is finer we do as before to avoid aliasing that would be introduced by grdsample.
%
%  HLINE     -> a line or patch handle to rectangle. Currently it's only used to fetch
%               the Mirone handles in the case the HANDLES arg is not provided.
%  TIPO      -> Is either 'grid' or 'image' (case insensitive) and select the operation to do
%  HANDLES   -> is the (optional) Mirone handles. If not provided HLINE must be a valid line handle
%
%  VARARGOUT -> [Z_rect, r_c] Where Z_rect is a sub-region of the host grid comprised of the
%               BoundingBox of the imported grid plus a padding of 6 host grids cells. The padd
%               cells were used to guaranty a smooth transition between host and inplanted grid
%               This is the default for usage in Mirone when this function is called from ImageCrop_CB().
%  VARARGOUT -> [Z] Where Z is the host array (single) with the implant grid already embeded
%
% In the IMPLANTIMAGE mode this function calls a subfunction that:
% An external image will be inplanted inside the zone defined by the rectangle whose handle is HLINE.

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

% $Id$

	if (nargin < 3)		% NOTE: in the IMPLANTGRID mode I'm not really using (YET) hLine for anything else
		if (~ishandle(hLine) && ~strcmp(get(hLine,'type'), 'line') && ~strcmp(get(hLine,'type'), 'patch'))
			error('TRANSPLANTS First argument is not a valid line handle')
		end
		handles = guidata(hLine);
	end

	% ------------------- If implant Image, call subfunction and return ------------------
	if (tipo(1) == 'i' || tipo(1) == 'I')
		transplat_img(handles, hLine),		return
	end
	% ------------------------------------------------------------------------------------

	varargout = cell(1, nargout);

	str = {'*.grd;*.nc;*.tif;*.tiff;*.jpg;*.jp2;*.png;*.gif;*.mat;*.cpt;*.hdf;*.img', ...
			'Files (*.grd,*.nc,*.tif,*.tiff,*.jpg,*.jp2,*.png,*.gif,*.mat,*.cpt,*.hdf,*.img)'; '*.*', 'All Files (*.*)'};
	[FileName,PathName] = put_or_get_file(handles,str,'Select file','get');
	if (isequal(FileName,0)),	return,		end				% User gave up
	fname = [PathName FileName];
	drv = aux_funs('findFileType', fname);
	if (isempty(drv)),		errordlg('Sory, don''t know what type of grid is this','Error'),	return,		end
	switch drv
		case 'gmt',			tipo = 'GMT';
		case 'mola',		tipo = 'MOLA';
		case 'las',			read_las(handles, fname);
		case {'geotif' 'ecw' 'dono'},		tipo = 'whatever';
		otherwise
			errordlg('Sorry, this file is either not a grid or not supported here.','Error'),	return
	end
	% Minimalist handles to get send to read_grid
	handlesInner.grdMaxSize = handles.grdMaxSize;
	handlesInner.ForceInsitu = handles.ForceInsitu;
	handlesInner.IamCompiled = handles.IamCompiled;
	handlesInner.path_tmp = handles.path_tmp;
	[Zinner, Xinner, Yinner, srsWKT, handlesInner] = read_grid(handlesInner, fname, tipo);
	if (isempty(Zinner)),	set(handles.figure1,'pointer','arrow'),		return,		end

	% Check that both grids intersect
	rect = aux_funs('rectangle_and', handles.head, handlesInner.head);
	if (isempty(rect))
		errordlg('The to-be-transplanted grid does not have any region in common with the base grid. Bye, Bye.','Error')
		set(handles.figure1,'pointer','arrow'),		return
	else
		h = line('xdata',rect(:,1), 'ydata',rect(:,2), 'Parent', handles.axes1, 'LineWidth', handles.DefLineThick, ...
			'Color',handles.DefLineColor);
		draw_funs(h,'line_uicontext')
		drawnow
	end

	% Compute a rectangle = InnerGrid + pad, but that doesn't get out of OuterGrid
	pad = 6;
	x_min = max(handlesInner.head(1) - pad * handles.head(8), handles.head(1));
	x_max = min(handlesInner.head(2) + pad * handles.head(8), handles.head(2));
	y_min = max(handlesInner.head(3) - pad * handles.head(9), handles.head(3));
	y_max = min(handlesInner.head(4) + pad * handles.head(9), handles.head(4));
	rect_outer = [x_min y_min (x_max - x_min) (y_max - y_min)];
	[X,Y,Z,head] = load_grd(handles);
	[Z_rect,r_c] = cropimg(head(1:2), head(3:4), Z, rect_outer, 'out_grid');

	% ------------------ See if Inner grid actually contains Outer grid ----------------------
	t = (abs([x_min x_max y_min y_max] - handles.head(1:4)) < 1e-8);
	if (all(t))				% Yes it does
		hd = handlesInner.head;		% Short of
        rect_t = [0 0 0 0];
		rect_t(1) = max(rect_outer(1) - 2*hd(8), hd(1));		    rect_t(2) = max(rect_outer(2) - 2*hd(9), hd(2));
		rect_t(3) = min(rect_outer(3) + 4*hd(8), diff(hd(1:2)));    rect_t(4) = min(rect_outer(4) + 4*hd(9), diff(hd(3:4)));
		[Zinner,rc] = cropimg(hd(1:2), hd(3:4), Zinner, rect_t, 'out_grid');
		Xinner = linspace(hd(1) + (rc(3)-1)*hd(8), hd(1) + (rc(4)-1)*hd(8), rc(4) - rc(3) + 1 );
		Yinner = linspace(hd(3) + (rc(1)-1)*hd(9), hd(3) + (rc(2)-1)*hd(9), rc(2) - rc(1) + 1 );
		handlesInner.head(1) = Xinner(1);		handlesInner.head(2) = Xinner(end);
		handlesInner.head(3) = Yinner(1);		handlesInner.head(4) = Yinner(end);
		% When no NaNs and implanting grid has cruder resolution (no aliasing risk) use grdsample and return
		if (~handlesInner.have_nans && hd(8) >= head(8))
			opt_R = sprintf('-R%.12g/%.12g/%.12g/%.12g', head(1:4));
			opt_I = sprintf('-I%.14g/%.14g', head(8:9));
			varargout{1} = c_grdsample(Zinner, handlesInner.head, opt_R, opt_I);
			if (nargout == 2)
				varargout{2} = r_c;
			end
			return
		end
		% Else we continue but Zinner is now a potentially much smaller array that original
		handlesInner.have_nans = grdutils(Zinner,'-N');		% Must recheck this now.
	end
	% --------------------------------------------------------------------- ----------------------

	if (handlesInner.have_nans)				% Case in which the inner grid has NaNs
		indNotNaN = ~isnan(Zinner);
		mask = img_fun('bwmorph',indNotNaN,'dilate');
		B    = img_fun('bwboundaries',mask, 8, 'noholes');
		B{1} = cvlib_mex('dp', B{1}, 0.1);	% Simplify line
		Bh   = img_fun('bwboundaries',mask, 8, 'holes');
		if (~isempty(Bh))					% So we have holes. Must find only those inside B{1}
			c = false(numel(Bh),1);
			for (k = 1:numel(Bh))
				Bh{k} = cvlib_mex('dp', Bh{k}, 0.1);		% Simplify line
				if (size(Bh{k},1) == size(B{1},1))			% One of them comes out repeated
					c(k) = true;	continue
				end
				IN = inpolygon(Bh{k}(:,1),Bh{k}(:,2), B{1}(:,1),B{1}(:,2));
				if (~all(IN)),	c(k) = true;	end		% This one is out
			end
			Bh(c) = [];
			if (~isempty(Bh))				% If we have any left its because they are inside holes
				B = [B; Bh];	clear Bh;
			end
		end
		x_inc = handlesInner.head(8);		y_inc = handlesInner.head(9);
		x_min = handlesInner.head(1);		y_min = handlesInner.head(3);
		if (handlesInner.head(7))			% Work in pixel registration
			x_min = x_min + x_inc/2;	y_min = y_min + y_inc/2;
		end
		y = (B{1}(:,1)-1)*y_inc + y_min;
		x = (B{1}(:,2)-1)*x_inc + x_min;
		set(h, 'xdata',x, 'ydata',y)		% Convert the rectangle into the outer polygon

	else
		% Compute another rect but this time with pad = 1 that will be used to clip inside
		x1 = max(handlesInner.head(1) - handles.head(8), handles.head(1));
		x2 = min(handlesInner.head(2) + handles.head(8), handles.head(2));
		y1 = max(handlesInner.head(3) - handles.head(9), handles.head(3));
		y2 = min(handlesInner.head(4) + handles.head(9), handles.head(4));
		x = [x1 x1 x2 x2 x1];		y = [y1 y2 y2 y1 y1];
	end

	X = linspace(head(1) + (r_c(3)-1)*head(8), head(1) + (r_c(4)-1)*head(8), r_c(4) - r_c(3) + 1);
	Y = linspace(head(3) + (r_c(1)-1)*head(9), head(3) + (r_c(2)-1)*head(9), r_c(2) - r_c(1) + 1);
	opt_R = sprintf('-R%.10f/%.10f/%.10f/%.10f', X(1), X(end), Y(1), Y(end));
	opt_I = sprintf('-I%.10f/%.10f',head(8),head(9));
	mask = ~(img_fun('roipoly_j',[x_min x_max],[y_min y_max],Z_rect,x,y));		% Mask at the outer grid resolution

	if (handlesInner.have_nans && numel(B) > 1)		% If we have holes in Inner grid, let the outer grid values survive there
		for (k = 2:numel(B))
			y = (B{k}(:,1)-1)*y_inc + y_min;
			x = (B{k}(:,2)-1)*x_inc + x_min;
			mask_t = (img_fun('roipoly_j',[x_min x_max],[y_min y_max],Z_rect,x,y));
			mask = mask | mask_t;
		end
	end
	ZZ = double(Z_rect(mask));			% It F... has to be

	[X,Y] = meshgrid(X,Y);
	XX = X(mask(:));
	YY = Y(mask(:));
	indNaN = isnan(ZZ);					% But are those NaN free?
	if (any(indNaN))
		ZZ(indNaN) = [];	XX(indNaN) = [];	YY(indNaN) = [];
	end

	[X,Y] = meshgrid(Xinner,Yinner);
	if (handlesInner.have_nans)
		XX = [XX(:); X(indNotNaN(:))];		clear X
		YY = [YY(:); Y(indNotNaN(:))];		clear Y
		ZZ = [ZZ(:); double(Zinner(indNotNaN(:)))];	clear Zinner	% Join outer and inner Zs
		opt_C = '-C5';		% Just a heuristic value
		if (head(8) < handlesInner.head(8))		% If host grid is finer than implanting grid
			% Than we must increase -C because data to grid may be too sparse and gaps appear.
			opt_C = sprintf('-C%d', 2*round(handlesInner.head(8) / head(8)));
		end
	else
		XX = [XX(:); X(:)];			clear X
		YY = [YY(:); Y(:)];			clear Y
		ZZ = [ZZ(:); double(Zinner(:))];	clear Zinner
		opt_C = ' ';
	end

	Z_rect = gmtmbgrid_m(XX,YY,ZZ(:), opt_R, opt_I, '-T0.25', '-Mz', opt_C);
	clear XX YY ZZ;

	if (nargout == 1)
		if (isa(Z,'single')),		Z(r_c(1):r_c(2),r_c(3):r_c(4)) = single(Z_rect);
		elseif (isa(Z,'int16')),	Z(r_c(1):r_c(2),r_c(3):r_c(4)) = int16(Z_rect);
		elseif (isa(Z,'uint16')),	Z(r_c(1):r_c(2),r_c(3):r_c(4)) = uint16(Z_rect);
		else						Z(r_c(1):r_c(2),r_c(3):r_c(4)) = single(Z_rect);
		end
		varargout{1} = Z;
	else
		varargout{1} = Z_rect;
		varargout{2} = r_c;
	end
	
% ---------------------------------------------------------------------------------------
function transplat_img(handles, h)
% Cirurgy Imagery operation. An external image will be inplanted inside the
% rectangular zone defined by the rectangle whose handle is h.
% Notice that we have to forsee the possibility of transplanting RGB images
% into indexed bg images and vice-versa.
%
% This one day will be donne correctly by overlaying the implant image in a new axes 

	hAxes = handles.axes1;		hImg = handles.hImg;
	out = implanting_img(hImg, h, get(hAxes,'xlim'), get(hAxes,'ylim'));
	if isempty(out),	return,		end
	zz = get(hImg,'CData');

	% Find if Implanting image needs to be ud fliped
	if(strcmp(get(hAxes,'XDir'),'normal') && strcmp(get(hAxes,'YDir'),'reverse'))
			flip = 0;
	else	flip = 1;
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
			for (i = 1:3)
				%ZI(:,:,i) = interp2(double(out.ip_img(:,:,i)),X,Y,'*cubic');
				ZI(:,:,i) = c_grdsample(single(out.ip_img(:,:,i)),head,opt_N);
			end
		else
			if isempty(out.ip_cmap)
				errordlg('Implanting image has no colormap. Don''t know what to do.','Sorry');  return
			end
			%ZI = interp2(double(out.ip_img),X,Y,'*cubic');
			ZI = c_grdsample(single(out.ip_img),head,opt_N);
		end
		if (flip),   ZI = flipdim(ZI,1);    end
	elseif (out.resizeIP == 10) % So pra nao funcionar (da erro na penultima linha)
		%nl_new = linspace(1,nl_bg,(out.bg_size_updated(1)));
		%nc_new = linspace(1,nc_bg,(out.bg_size_updated(2)));
		%[X,Y] = meshgrid(nc_new,nl_new);
		%if (~indexed_bg)                            % Background image is of RGB type
			%for (i=1:3)
				%zz(:,:,i) = interp2(double(zz(:,:,i)),X,Y,'*cubic');
				%zz(:,:,i) = c_grdsample(zz(:,:,i),head,opt_N);
			%end
		%else
			%zz = interp2(double(zz),X,Y,'*cubic');
			%zz = c_grdsample(zz,head,opt_N);
		%end
		%if (flip)    out.ip_img = flipdim(out.ip_img,1);    end
	end

	if (indexed_ip && ~indexed_bg)				% Implanting indexed image on a RGB bg image
		%I = ind2rgb8(out.ip_img,out.ip_cmap);	% Transform implanting image to RGB
	elseif (indexed_ip && indexed_bg)			% Shit, both ip & bg images are indexed. We have to RGB them
		zz = ind2rgb8(zz,colormap);
		%I = ind2rgb8(out.ip_img,out.ip_cmap);
	elseif (~indexed_ip && ~indexed_bg)			% Nice, nothing to do
	elseif (~indexed_ip && indexed_bg)			% Implanting RGB image on a indexed bg image.
		zz = ind2rgb8(zz,colormap);				% Transform bg image to RGB
	end

	zz(out.r_c(1):out.r_c(2), out.r_c(3):out.r_c(4), :) = uint8(ZI);
	set(hImg,'CData',zz)

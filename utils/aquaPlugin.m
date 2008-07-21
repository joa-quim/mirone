function aquaPlugin(handles)
% Plugin function that is called by Aquamoto. Use this function to write custom code
% to solve particular problems taking advantage from the fact that a LOT of information
% about the netCDF files is accessible here. There are no instructions/manual but you
% can learn by studing the functions on aquamoto.m file or the (not so clean) working
% examples below. 

	if ( isempty(handles.fname) )
		errordlg('Hey Lou. What about a walk on the Wild Side? Maybe you''ll find a little file there that you can use here!','Chico clever')
		return
	end

	do_return = true;
	casos = {'zonal' ...			% 1 - Compute zonal means
			'tvar' ...				% 2 - Compute the Temp time rate of a file with anual means by fit of a straight line (Load entire file in memory)
			'yearMean' ...			% 3 - Compute yearly averages from monthly data
			'yearMeanFlag' ...		% 4 - Compute yearly averages from monthly data but checked against a quality flag file
			'polygAVG' ...			% 5 - Compute averages of whatever inside polygons (if any)
			'flagsStats' ...		% 6 - Compute per/pixel anual or mounth counts of pixel values with a quality >= flag
			'do_math' ...			% 7 - Perform some basic agebraic operations with the 3D planes
			};

	qual = casos{4};		% <== Active selection

	switch qual
		case 'zonal'				% case 1
			integ_lon = true;
			dlat = 0.5;
			have_polygon = true;
			if (have_polygon);
				fname = 'plataforma_poly.dat';
				S = load(fname);
				x = S(:,1)+360;		y = S(:,2);
			else
				x = [];		y = [];
			end

			zonal(handles, dlat, integ_lon, have_polygon, x, y)
		case 'tvar'					% case 2
			calcGrad(handles) 
		case 'yearMean'				% case 3
			ano = 1:12;				% Compute yearly means
			calc_yearMean(handles, ano)
		case 'yearMeanFlag'			% case 4
			ano = 1:12;				% Compute yearly means
			fname = 'C:\SVN\mironeWC\qual_85_07.nc';
			quality = 7;			% Retain only values of quality >= this
			pintAnoes = true;		% If true instruct to fill holes <= nCells
			nCells = 200;			% Holes (bad data) smaller than this are filled by interpolation
			calc_yearMean(handles, ano, fname, quality, pintAnoes, nCells)
		case 'polygAVG'				% case 5
			calc_polygAVG(handles)
		case 'flagsStats'			% case 6
			ano = 1:12;				% Compute yearly stats
			opt = 'per_year';		% Make the counting on a per year basis
			opt = '';				% Make the counting on a per month basis
			calc_flagsStats(handles, ano, 7, opt)
		case 'do_math'				% case 7
			opt = 'sum';			% Sum all layers (only operation for the time beeing)
			do_math(handles, opt)
		otherwise
			do_return = false;		% Not finish (run the code below)
	end
	
	if (do_return),		return,		end		% Finish work here
	
	% Old part of code that computes yearly rate of change from monthly data

	z_id = handles.netcdf_z_id;
	s = handles.nc_info;			% Retrieve the .nc info struct
	rows = s.Dataset(z_id).Size(end-1);
	cols = s.Dataset(z_id).Size(end);
	
	anos = handles.number_of_timesteps / 12;
	
	aguentabar(0,'title','Calcula as medias anuais ')
	Tmed = zeros([rows, cols, anos]);	% Temp media para cada um dos anos
	mn = 0;
	warning off MATLAB:divideByZero
	for (m = 1:anos)
		contanoes = zeros(rows, cols);
		for (n = 1:12)
			Z = nc_funs('varget', handles.fname, s.Dataset(z_id).Name, [mn 0 0], [1 rows cols]);
			ind = isnan(Z);
			Z(ind) = 0;				% transmuta os anoes
			contanoes = contanoes + ~ind;
			Tmed(:,:,m) = Tmed(:,:,m) + double(Z);
			mn = mn + 1;
		end
		Tmed(:,:,m) = Tmed(:,:,m) ./ contanoes;
	end

	clear Z
	Tmed(Tmed == 0) = nan;			% Repoe os anoes
	Tvar = zeros(rows,cols);

	aguentabar(0,'title','Calcula a variaçao temporal','CreateCancelBtn')
	
	x = (0:anos-1)';
	for (m = 1:rows)
		for (n = 1:cols)
			p = polyfit(x,squeeze(Tmed(m,n,:)),1);
			Tvar(m,n) = p(1);
		end
		h = aguentabar(m/rows);
		if (isnan(h)),	break,	end
	end
	if (isnan(h)),	return,		end
	
	clear Tmed
	Tvar = single(Tvar);
	
	tmp.head = handles.head;
	zz = grdutils(Tvar,'-L');  tmp.head(5:6) = [zz(1) zz(2)];
	tmp.X = linspace(tmp.head(1),tmp.head(2),cols);
	tmp.Y = linspace(tmp.head(3),tmp.head(4),rows);
	tmp.name = 'Time rate (deg/year)';
	mirone(Tvar, tmp)
	
% ----------------------------------------------------------------------
function zonal(handles, dlat, integ_lon, have_polygon, x, y)
% Compute zonal mean from a multi-layer file

	z_id = handles.netcdf_z_id;
	s = handles.nc_info;			% Retrieve the .nc info struct
	rows = s.Dataset(z_id).Size(end-1);
	cols = s.Dataset(z_id).Size(end);

	% Build the vectors to deal with the zonal integration
	if (integ_lon)
		N_spatialSize = rows;		% Number of points in the spatial dim
		integDim = 2;							% Dimension along which we are going to integrate
		vecD = floor(handles.head(3)):dlat:(handles.head(4));
		Y = linspace(handles.head(3),handles.head(4), N_spatialSize);
	else
		N_spatialSize = cols;
		integDim = 1;
		vecD = floor(handles.head(1)):dlat:ceil(handles.head(2));
		Y = linspace(handles.head(1),handles.head(2), N_spatialSize);
	end
	nStripes = numel(vecD) - 1;
	indStripe = ones(numel(vecD),1);

	for (k = 2:nStripes)
		ind = find(Y >= vecD(k));
		if (~isempty(ind)),		indStripe(k) = ind(1);	end
	end
	indStripe(end) = N_spatialSize;

	aguentabar(0,'title','Computing zonal means','CreateCancelBtn')

	nSeries = handles.number_of_timesteps;
	allSeries = zeros(nStripes, nSeries);
	mask = [];
	for (k = 1:nSeries)
		Z = nc_funs('varget', handles.fname, s.Dataset(z_id).Name, [k-1 0 0], [1 rows cols]);
		this_has_nans = false;

		% NaNify polygon exterior points?
		if (have_polygon && k == 1)
			mask = img_fun('roipoly_j',handles.head(1:2),handles.head(3:4),double(Z),x,y);
		end
		if (~isempty(mask))
			Z(~mask) = NaN;
		end

		ind = isnan(Z);
		if (any(ind(:))),		Z(ind) = 0;		this_has_nans = true;		end

		tmp = sum(Z,integDim);			% Add along integration dim
		if (integ_lon),		N = size(Z,2);
		else				N = size(Z,1);
		end
		if (this_has_nans)
			tmp2 = sum(ind,integDim);
			tmp = tmp ./ (N - tmp2);
		else
			tmp = tmp / N;
		end
		% Now add all inside each stripe
		for (m = 1:nStripes)
			tmp2 = tmp( indStripe(m):indStripe(m+1) );
			allSeries(m,k) = sum(tmp2) / numel(tmp2);
		end

		h = aguentabar(k/nSeries);
		if (isnan(h)),	break,	end
	end
	if (isnan(h)),	return,		end

	allSeries = single(allSeries);
	zz = grdutils(allSeries,'-L');
	head = [1 nSeries vecD(1) vecD(end) zz(1) zz(2) 0 1 dlat];
	tmp.X = 1:nSeries;		tmp.Y = linspace( (vecD(1)+dlat/2), (vecD(end)-dlat/2), nStripes );
	tmp.head = [head(1:2) tmp.Y(1) tmp.Y(end) head(5:end)];
	tmp.geo = 0;			tmp.name = 'Zonal integration';
	mirone(allSeries, tmp)
	
% ----------------------------------------------------------------------
function calcGrad(handles)
% Calcula o gradiente de um fiche ja com as medias anuais por ajuste de um recta (Loada o fiche todo na memoria)

	z_id = handles.netcdf_z_id;
	s = handles.nc_info;			% Retrieve the .nc info struct
	rows = s.Dataset(z_id).Size(end-1);
	cols = s.Dataset(z_id).Size(end);

	anos = handles.number_of_timesteps;
	Tmed = zeros([rows, cols, anos]);	% Temp media para cada um dos anos
	for (m = 1:anos)
		Tmed(:,:,m) = nc_funs('varget', handles.fname, s.Dataset(z_id).Name, [m-1 0 0], [1 rows cols]);
	end

	aguentabar(0,'title','Calcula o gradiente','CreateCancelBtn')
	
	Tvar = zeros(rows,cols);
	x = (0:anos-1)';
	for (m = 1:rows)
		for (n = 1:cols)
			p = polyfit(x,squeeze(Tmed(m,n,:)),1);
			Tvar(m,n) = p(1);
		end
		h = aguentabar(m/rows);
		if (isnan(h)),	break,	end
	end
	if (isnan(h)),	return,		end

	clear Tmed
	Tvar = single(Tvar);
	
	tmp.head = handles.head;
	zz = grdutils(Tvar,'-L');  tmp.head(5:6) = [zz(1) zz(2)];
	tmp.X = linspace(tmp.head(1),tmp.head(2),cols);
	tmp.Y = linspace(tmp.head(3),tmp.head(4),rows);
	tmp.name = 'Time gradient (deg/year)';
	mirone(Tvar, tmp)

% ----------------------------------------------------------------------
function calc_yearMean(handles, months, fname2, flag, pintAnoes, nCells)
% Calcula media anuais a partir de dados mensais
% MONTHS 	is a vector with the months uppon which the mean is to be computed
%		example: 	months = 1:12		==> Computes yearly mean
%					months = 6:8		==> Computes June-July-August seazonal means
%
% OPTIONS:
% FNAME2 	name of a netCDF file with quality flags. Obviously this file must be of
% 			the same size as the series under analysis.
% FLAG		Threshold quality value. Only values of quality >= FLAG will be taken into account
% PINTANOES	Logical that if true instruct to fill holes <= NCELLS

	txt1 = 'netCDF grid format (*.nc,*.grd)';
	txt2 = 'Select output netCDF grid';
	[FileName,PathName] = put_or_get_file(handles,{'*.nc;*.grd',txt1; '*.*', 'All Files (*.*)'},txt2,'put');
	if isequal(FileName,0),		return,		end
	[pato, fname, EXT] = fileparts(FileName);
	if (isempty(EXT)),		FileName = [fname '.nc'];	end
	grd_out = [PathName FileName];

	z_id = handles.netcdf_z_id;
	s = handles.nc_info;				% Retrieve the .nc info struct
	rows = s.Dataset(z_id).Size(end-1);
	cols = s.Dataset(z_id).Size(end);
	do_flags = false;

	if (nargin == 1),		months = 1:12;		end		% Default to yearly means

	if (nargin > 2)			% We have a quality-flag ghost file to check
		s_flags = nc_funs('info',fname2);
		[X,Y,Z,head,misc] = nc_io(fname2,'R');
		z_id_flags = misc.z_id;
		if ~(numel(head) == 9 && isfield(misc,'z_id'))
			errordlg(['Blheak!! ' fname2 ' is is not a file with presumably with quality flags. By'],'Error'),	return
		end
		if (numel(misc.z_dim) <= 2)
			errordlg(['Ghrrr!! The ' fname2 ' is is not a 3D file. By'],'Error'),		return
		end
		if (misc.z_dim(1) < handles.number_of_timesteps)
			errordlg(['Buhhuu!! The quality flags file has less "planes" than the-to-be-flagged-file. By'],'Error'),	return
		end
		if (~isequal([rows cols], [s_flags.Dataset(z_id_flags).Size(end-1) s_flags.Dataset(z_id_flags).Size(end)]))
			errordlg(['Buhhuu!! quality flags and the-to-be-flagged-file have not the same size. By'],'Error'),		return
		end
		do_flags = true;
		
		if (nargin == 3),		flag = 7;		end			% If not provided, defaults to best quality
	end

	handles.geog = 1;		handles.was_int16 = 0;		handles.computed_grid = 0;

	anos = handles.number_of_timesteps / 12;
	
	aguentabar(0,'title','Calcula as medias anuais','CreateCancelBtn')
	Tmed = zeros([rows, cols]);			% Temp media para cada um dos anos
	warning off MATLAB:divideByZero

	anos = 2;
	for (m = 1:anos)
		contanoes = zeros(rows, cols);
		for (n = months)
			mn = (m - 1)*12 + n - 1;
			Z = nc_funs('varget', handles.fname, s.Dataset(z_id).Name, [mn 0 0], [1 rows cols]);

			if (do_flags)
				Z_flags = nc_funs('varget', fname2, s_flags.Dataset(z_id_flags).Name, [mn 0 0], [1 rows cols]);
				Z(Z_flags < flag) = NaN;
			end

			ind = isnan(Z);

			if (pintAnoes && any(ind(:)))
				Z = inpaint_nans(handles, Z, ind, nCells);		% Select interp method inside inpaint_nans()
				ind = isnan(Z);
			end

			Z(ind) = 0;				% transmuta os anoes
			contanoes = contanoes + ~ind;
			Tmed(:,:) = Tmed(:,:) + double(Z);
		end
		Tmed(:,:) = Tmed(:,:) ./ contanoes;
		tmp = single(Tmed(:,:));
		
		if (m == 1),		nc_io(grd_out,sprintf('w%d/time',anos), handles, reshape(tmp,[1 size(tmp)]))
		else				nc_io(grd_out, sprintf('w%d', m-1), handles, tmp)
		end
		
		h = aguentabar(m/anos);
		if (isnan(h)),	break,	end
		
		Tmed = Tmed * 0;			% Reset it to zeros
	end

% ----------------------------------------------------------------------
function Z = inpaint_nans(handles, Z, bw, nCells)
% Interpolate holes in Z that are smaller than NCELLS in size
%
% BW is a logicall array which maps where Z has NaNs
	
	if (nargin == 3),	nCells = 100;	end

	use_surface = true;
	use_bicubic = false;
	pad = 4;				% Number of cells to increase the rectangle that encloses each hole
	head = handles.head;
	[rows, cols] = size(Z);
	
	% Retain only <= handles.nCells sized of connected groups
	bw2 = img_fun('bwareaopen', bw, nCells);
	bw = xor(bw, bw2);
	clear bw2;
	
	B = img_fun('find_holes',bw);

	if (use_surface)
		opt_I = sprintf('-I%.10f/%.10f',head(8),head(9));
	else
		opt_I = ' ';
	end

	n_buracos = numel(B);
	for (i = 1:n_buracos)
		% Get rectangles arround each hole
		x_min = min(B{i}(:,2));			x_max = max(B{i}(:,2));
		y_min = min(B{i}(:,1));			y_max = max(B{i}(:,1));
		x_min = max(1,x_min-pad);		x_max = min(x_max+pad,cols);
		y_min = max(1,y_min-pad);		y_max = min(y_max+pad,rows);
		x_min = head(1) + (x_min-1)*head(8);    x_max = head(1) + (x_max-1)*head(8);
		y_min = head(3) + (y_min-1)*head(9);    y_max = head(3) + (y_max-1)*head(9);

		rect_crop = [x_min y_min (x_max-x_min) (y_max-y_min)];
		[Z_rect, r_c]  = cropimg(head(1:2),head(3:4),Z,rect_crop,'out_grid');
		bw_rect = cropimg(head(1:2),head(3:4),bw,rect_crop,'out_grid');
		Z_rect = double(Z_rect);      % It has to be (GHRRRRRRRRRRRRR)

		X = x_min:head(8):x_max;	Y = y_min:head(9):y_max;
		[XX,YY] = meshgrid(X,Y);
		XX(bw_rect) = [];			YY(bw_rect) = [];		Z_rect(bw_rect) = [];

		if (use_surface)
			opt_R = sprintf('-R%.10f/%.10f/%.10f/%.10f', X(1), X(end), Y(1), Y(end));
			%Z_rect = surface_m( XX(:), YY(:), Z_rect(:), opt_R, opt_I, '-T.25' );
			Z_rect = gmtmbgrid_m( XX(:), YY(:), Z_rect(:), opt_R, opt_I, '-T.25', '-Mz' );
		elseif (use_bicubic)
			Z_rect = griddata_j(XX(:), YY(:), Z_rect(:), X, Y', 'cubic');
		else
			Z_rect = griddata_j(XX(:), YY(:), Z_rect(:), X, Y', 'linear');
		end

		% Inprint the processed rectangle back into orig array
		if (isa(Z,'single')),		Z(r_c(1):r_c(2),r_c(3):r_c(4)) = single(Z_rect);
		elseif (isa(Z,'int16')),	Z(r_c(1):r_c(2),r_c(3):r_c(4)) = int16(Z_rect);
		elseif (isa(Z,'uint16')),	Z(r_c(1):r_c(2),r_c(3):r_c(4)) = uint16(Z_rect);
		else						Z(r_c(1):r_c(2),r_c(3):r_c(4)) = single(Z_rect);
		end
	end
	%clear gmtmbgrid_m

% ----------------------------------------------------------------------
function calc_polygAVG(handles)
	
	if (~ishandle(handles.handMir.figure1))		% No insult. Just quit
		return
	end

	hLine = findobj(handles.handMir.axes1,'Type','line');
	hLine = [hLine; findobj(handles.handMir.axes1,'Type','patch')];

	N = 1;
	for (k = 1:numel(hLine))
		x = get(hLine(k),'XData');   y = get(hLine(k),'YData');
		if (numel(x) >= 3 && x(1) == x(end) && y(1) == y(end) )
			polys{N} = [x(:) y(:)];
			N = N + 1;
		end
	end
	N = N - 1;		% There was one too much incement above

	if (N == 0)
		errordlg('Fiu Fiu! No closed polygons to compute whaterver average value inside. Bye.','Error')
		return
	end
	
	z_id = handles.netcdf_z_id;
	s = handles.nc_info;			% Retrieve the .nc info struct
	rows = s.Dataset(z_id).Size(end-1);
	cols = s.Dataset(z_id).Size(end);
	nLayers = handles.number_of_timesteps;
	avg = zeros(nLayers,N);
	THRESH = 0.5;					% Minimum percentage of valid points inside poly

	aguentabar(0,'title','Calcula as medias poligonais','CreateCancelBtn')

	for (m = 1:nLayers)				% Loop over layers ensemble
		Z = nc_funs('varget', handles.fname, s.Dataset(z_id).Name, [m-1 0 0], [1 rows cols]);

		for (k = 1:N)				% Loop over polygons
			x = polys{k}(:,1);				y = polys{k}(:,2);
			xp(1) = min(x);     xp(2) = max(x);
			yp(1) = min(y);     yp(2) = max(y);
			rect_crop = [xp(1) yp(1) (xp(2) - xp(1)) (yp(2) - yp(1))];
			x_lim = [xp(1) xp(2)];		y_lim = [yp(1) yp(2)];

			% Extrai um rect que englobe o poligono para poupar na conta da mascara
			[Z_rect,r_c] = cropimg(handles.head(1:2),handles.head(3:4),Z,rect_crop,'out_grid');
			mask = img_fun('roipoly_j',x_lim,y_lim,Z_rect,x,y);

			% Test for a minimum of valid elements inside polygon
			zz = Z_rect(mask);
			zz = zz(:);
			ind = isnan(zz);
			if (~any(ind))
				avg(m,k) = sum(double(zz)) / numel(zz);
			else			% Accept/Reject based on % of valid numbers
				nAnoes = sum(ind);		nInPoly = numel(zz);
				if ( nAnoes / nInPoly < THRESH )
					zz = zz(~ind);
					avg(m,k) = sum(double(zz)) / numel(zz);
				end
			end
		end
		h = aguentabar(m/nLayers);
		if (isnan(h)),	break,	end
		
	end

	if (isnan(h)),	return,		end
	
% 	% -------------- We still need to determine polygon's area to compute the areal average
% 	for (k = 1:N)
% 		x = polys{k}(:,1);		y = polys{k}(:,2);
% 		if (handles.geog)
% 			area = area_geo(y,x);    % Area is reported on the unit sphere
% 			area = area * 4 * pi * (6371005^2);
% 		else
% 			area = polyarea(x,y);   % Area is reported in map user unites
% 		end
% 		avg(:,k) = avg(:,k) / (area * 1e-6);	% per km^2
% 	end
	
	% --------------- Now finaly save the result in a file	------------------
	[FileName,PathName] = put_or_get_file(handles,{'*.dat;*.DAT','ASCII file'; '*.*', 'All Files (*.*)'},'Output file','put');
	if isequal(FileName,0),		return,		end
	[PATH,FNAME,EXT] = fileparts([PathName FileName]);
	if isempty(EXT),	fname = [PathName FNAME '.dat'];
	else				fname = [PathName FNAME EXT];
	end
	
	%Open and write to ASCII file
	if (ispc),		fid = fopen(fname,'wt');
	elseif (isunix),fid = fopen(fname,'w');
	else			error('aquamoto: Unknown platform.');
	end

	% Calculate a a rough polygon centroid
	centro = zeros(N,2);
	for (k = 1:N)				% Loop over polygons
		x = polys{k}(:,1);				y = polys{k}(:,2);
		centro(k,:) = mean(polys{k});
	end
	
	fprintf(fid, ['#  \t', repmat('%g(X)\t', [1,N]) '\n'], centro(:,1));
	fprintf(fid, ['# T\t', repmat('%g(Y)\t', [1,N]) '\n'], centro(:,2));

	try			t = handles.time;		t = t(:);		% Layers's times
	catch		t = (1:size(avg,1))';
	end

	fprintf(fid,['%.2f\t' repmat('%f\t',[1,N]) '\n'], [t avg]');
	fclose(fid);

% ----------------------------------------------------------------------
function calc_flagsStats(handles, months, flag, opt)
% Compute per/pixel anual counts of pixel values with a quality >= flag
% Perfect locations will have a count of 12. Completely cloudy => count = 0.
%
% MONTHS 	is a vector with the months uppon which the mean is to be computed
%		example: 	months = 1:12		==> Computes yearly mean
%					months = 6:8		==> Computes June-July-August seazonal means
% FLAG		Threshold quality value. Only values of quality >= FLAG will be taken into account
% OPT		== 'per_year' output a file with N years planes (count per year)
%			otherwise ouputs a file with N months planes (each plane has a mounthly count)

	if (nargin < 3),		flag = 7;	opt = 'per_year';		end			% If not provided, defaults to best quality
	if (nargin < 4),		opt = 'per_year';		end

	txt1 = 'netCDF grid format (*.nc,*.grd)';
	txt2 = 'Select output netCDF grid';
	[FileName,PathName] = put_or_get_file(handles,{'*.nc;*.grd',txt1; '*.*', 'All Files (*.*)'},txt2,'put');
	if isequal(FileName,0),		return,		end
	[pato, fname, EXT] = fileparts(FileName);
	if (isempty(EXT)),		FileName = [fname '.nc'];	end
	grd_out = [PathName FileName];

	s = handles.nc_info;				% Retrieve the .nc info struct
	rows = s.Dataset(handles.netcdf_z_id).Size(end-1);
	cols = s.Dataset(handles.netcdf_z_id).Size(end);

	handles.geog = 1;		handles.was_int16 = 0;		handles.computed_grid = 0;

	anos = handles.number_of_timesteps / 12;

	aguentabar(0,'title','Compute flag quality counts','CreateCancelBtn')
	goodCount = zeros([rows, cols]);			% 

	if (strcmp(opt, 'per_year'))		% Make the counting on a per year basis
		for (m = 1:anos)
			contanoes = zeros(rows, cols);
			for (n = months)
				mn = (m - 1)*12 + n - 1;
				Z = nc_funs('varget', handles.fname, s.Dataset(handles.netcdf_z_id).Name, [mn 0 0], [1 rows cols]);
				ind = Z < flag;
				Z(ind) = 0;
				Z(~ind) = 1;				% Bellow threshold quality are set to zero
				contanoes = contanoes + ~ind;
				goodCount(:,:) = goodCount(:,:) + double(Z);
			end
			tmp = int16(goodCount(:,:));
	
			if (m == 1),	nc_io(grd_out,sprintf('w%d/time',anos), handles, reshape(tmp,[1 size(tmp)]))
			else			nc_io(grd_out, sprintf('w%d', m-1), handles, tmp)
			end
			
			h = aguentabar(m/anos);
			if (isnan(h)),	break,	end
			goodCount = goodCount * 0;			% Reset it to zeros
		end

	else			% Make the counting on a per month basis

		for (n = months)
			contanoes = zeros(rows, cols);
			for (m = 1:anos)
				mn = (m - 1)*12 + n - 1;
				Z = nc_funs('varget', handles.fname, s.Dataset(handles.netcdf_z_id).Name, [mn 0 0], [1 rows cols]);
				ind = Z < flag;
				Z(ind) = 0;
				Z(~ind) = 1;				% Bellow threshold quality are set to zero
				contanoes = contanoes + ~ind;
				goodCount(:,:) = goodCount(:,:) + double(Z);
			end
			tmp = int16(goodCount(:,:));
	
			if (n == 1),	nc_io(grd_out,sprintf('w%d/time',numel(months)), handles, reshape(tmp,[1 size(tmp)]))
			else			nc_io(grd_out, sprintf('w%d', n-1), handles, tmp)
			end
			
			h = aguentabar(n/numel(months));
			if (isnan(h)),	break,	end	
			goodCount = goodCount * 0;			% Reset it to zeros
		end
	end

% ----------------------------------------------------------------------
function do_math(handles, opt)
% Perform some basic agebraic operations on 3D planes
%
% OPT = 'sum'	=> add all layers

	s = handles.nc_info;				% Retrieve the .nc info struct
	rows = s.Dataset(handles.netcdf_z_id).Size(end-1);
	cols = s.Dataset(handles.netcdf_z_id).Size(end);
	nLayers = handles.number_of_timesteps;

	handles.geog = 1;		handles.was_int16 = 0;		handles.computed_grid = 0;

	if (strcmp(opt, 'sum'))
		soma = nc_funs('varget', handles.fname, s.Dataset(handles.netcdf_z_id).Name, [0 0 0], [1 rows cols]);
		is_int8 = isa(soma, 'int8');		is_uint8 = isa(soma, 'uint8');
		if (is_int8 || is_uint8),		soma = int16(soma);		end			% To avoid ovelflows
		for (m = 2:nLayers)
			Z = nc_funs('varget', handles.fname, s.Dataset(handles.netcdf_z_id).Name, [m-1 0 0], [1 rows cols]);
			if (is_int8 || is_uint8),		Z = int16(Z);		end			% To avoid ovelflows
			cvlib_mex('add',soma, Z);
		end
	end
	
	tmp.head = handles.head;
	if (isa(soma,'single'))
		zz = grdutils(soma,'-L');  tmp.head(5:6) = [zz(1) zz(2)];		% Singles & NaNs = BUGs in R13
	else
		tmp.head(5:6) = [double(min(soma(:))) double(max(soma(:)))];
	end
	tmp.X = linspace(tmp.head(1),tmp.head(2),cols);
	tmp.Y = linspace(tmp.head(3),tmp.head(4),rows);
	tmp.name = 'Computed grid';
	mirone(soma, tmp)
	
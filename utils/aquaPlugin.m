function aquaPlugin(handles)
% Plugin function that is called by Aquamoto. Use this function to write custom code
% to solve particular problems taking advantage from the fact that a LOT of information
% about the netCDF files is accessible here. There are no instructions/manual but you
% can learn by studing the functions on aquamoto.m file or the (not very clean) working
% examples below. 

% Compute the mean anual temperature gradient
	if ( isempty(handles.fname) )
		errordlg('Hey Lou. What about a walk on the Wild Side? Maybe you''ll find a little file there that you can use here!','Chico clever')
		return
	end

	qual = 'tvar';
	do_return = true;
	switch qual
		case 'zonal'
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
		case 'tvar'
			calcGrad(handles)		% Compute the Temp time rate of a file with anual means by fit of a straight line (Load entire file in memory)
		case 'yearMean'
			calc_yearMean(handles)	% Compute yearly averages from monthly data
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
function calc_yearMean(handles)
% Calcula media anuais a partir de dados mensais

	txt1 = 'netCDF grid format (*.nc,*.grd)';
	txt2 = 'Select output netCDF grid';
	[FileName,PathName] = put_or_get_file(handles,{'*.nc;*.grd',txt1; '*.*', 'All Files (*.*)'},txt2,'put');
	if isequal(FileName,0),		return,		end
	[pato, fname, EXT] = fileparts(FileName);
	if (isempty(EXT)),		FileName = [fname '.nc'];	end
	grd_out = [PathName FileName];
% 	grd_out = 'PF_avgYear_85-06.nc';

	z_id = handles.netcdf_z_id;
	s = handles.nc_info;				% Retrieve the .nc info struct
	rows = s.Dataset(z_id).Size(end-1);
	cols = s.Dataset(z_id).Size(end);

	handles.geog = 1;
	handles.was_int16 = 0;
	handles.computed_grid = 0;

	anos = handles.number_of_timesteps / 12;
	
	aguentabar(0,'title','Calcula as medias anuais','CreateCancelBtn')
	Tmed = zeros([rows, cols]);			% Temp media para cada um dos anos
	mn = 0;
	warning off MATLAB:divideByZero
	for (m = 1:anos)
		contanoes = zeros(rows, cols);
		for (n = 1:12)
			Z = nc_funs('varget', handles.fname, s.Dataset(z_id).Name, [mn 0 0], [1 rows cols]);
			ind = isnan(Z);
			Z(ind) = 0;				% transmuta os anoes
			contanoes = contanoes + ~ind;
			Tmed(:,:) = Tmed(:,:) + double(Z);
			mn = mn + 1;
		end
		Tmed(:,:) = Tmed(:,:) ./ contanoes;
		tmp = single(Tmed(:,:));
		
		if (m == 1)
			nc_io(grd_out,sprintf('w%d/time',anos), handles, reshape(tmp,[1 size(tmp)]))
		else
			nc_io(grd_out, sprintf('w%d', m-1), handles, tmp)
		end
		
		h = aguentabar(m/anos);
		if (isnan(h)),	break,	end
		
		Tmed = Tmed * 0;			% Reset it to zeros
	end

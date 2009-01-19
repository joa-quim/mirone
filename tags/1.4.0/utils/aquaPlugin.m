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

	casos = {'zonal' ...			% 1 - Compute zonal means
			'tvar' ...				% 2 - Compute the Temp time rate of a file with anual means by fit of a straight line (Load entire file in memory)
			'yearMean' ...			% 3 - Compute yearly averages from monthly data
			'yearMeanFlag' ...		% 4 - Compute yearly averages from monthly data but checked against a quality flag file
			'polygAVG' ...			% 5 - Compute averages of whatever inside polygons (if any)
			'flagsStats' ...		% 6 - Compute per/pixel anual or month counts of pixel values with a quality >= flag
			'pass_by_count' ...		% 7 - Check the curently active 3D file against a count file
			'do_math' ...			% 8 - Perform some basic agebraic operations with the 3D planes
			'conv2vtk' ...			% 9 - Convert a 3D netCDF file into a VTK format
			};

	qual = casos{2};		% <== Active selection

	switch qual
		case 'zonal'				% case 1
			integ_lon = true;
			dlat = 1.0;
 			have_polygon = true;
			trends = true;			% If true compute the trends (per stripe) of the zonal integration
			fname = 'plataforma_poly.dat';		% If this name is uncorrect, another will be asked
			fnam2 = 'plataforma_offset_poly.dat';		% If it exists, compute difference of zonal integrations
			fnam2 = 'poly_largo.dat';	%fnam2= [];
			zonal(handles, dlat, integ_lon, trends, have_polygon, fname, fnam2)
		case 'tvar'					% case 2
			calcGrad(handles) 
		case 'yearMean'				% case 3
			ano = 1:12;				% Compute yearly means
			calc_yearMean(handles, ano)
		case 'yearMeanFlag'			% case 4
			ano = [1:6 10:12];			% Compute yearly (ano = 1:12) or seasonal means (ano = start_month:end_month)
			fname  = 'C:\SVN\mironeWC\qual_85_07.nc';
			quality = 6;			% Retain only values of quality >= this
			pintAnoes = true;		% If true instruct to fill holes <= nCells
			nCells = 1000;			% Holes (bad data) smaller than this are filled by interpolation
			splina = true;			% Fill missing monthly data by a spline interpolation taken over two years
			% Where to save track of filled holes. Ignored if pintAnoes = false OR fname3 = []
			fname3 = 'C:\SVN\mironeWC\qual7_85_07_Interp200_Q6.nc';
			fname3 = [];
			calc_yearMean(handles, ano, fname, quality, pintAnoes, nCells, fname3, splina)
		case 'polygAVG'				% case 5
			calc_polygAVG(handles)
		case 'flagsStats'			% case 6
			ano = 1:12;				% Compute yearly stats
			opt = '';				% Make the counting on a per month basis
			opt = 'per_year';		% Make the counting on a per year basis
			calc_flagsStats(handles, ano, 7, opt)
		case 'pass_by_count'		% case 7
			count = 11;
			fname = 'C:\SVN\mironeWC\countPerYear_flag7_Interp200.nc';
			pass_by_count(handles, count, fname)
		case 'do_math'				% case 8
			opt = 'sum';			% Sum all layers (only operation for the time beeing)
			do_math(handles, opt)
		case 'conv2vtk'				% case 9
			write_vtk(handles)
	end

% ----------------------------------------------------------------------
function out = zonal(handles, dlat, integ_lon, do_trends, have_polygon, fname, fname2)
% Compute zonal means from a multi-layer file
% DLAT 			width of the box in the direction orthogonal to INTEG_LON
% INTEG_LON 	If true, integration is done along longitude
% DO_TRENDS		If false compute zonal integrations. Otherwise compute trends of the zonal integrations (per DLAT)
% HAVE_POLYGON	If true limit the analisys to the are delimited by the polygon stored in file FNAME 
% FNAME2		Optional polygon file. If it points to a vald file. This function is called twice and results are subtracted

	if (have_polygon)
		if (exist(fname,'file') ~= 2)
			[FileName,PathName] = put_or_get_file(handles, ...
				{'*.dat;*.DAT', 'Data files (*.dat)';'*.*', 'All Files (*.*)'},'Enter polygon file','get');
			if (isequal(FileName,0)),		return,		end
			fname = [PathName FileName];
		end
		S = load(fname);
		x = S(:,1);		y = S(:,2);
	else
		x = [];			y = [];
	end

	[z_id, s, rows, cols] = get_ncInfos(handles);

	% Build the vectors to deal with the zonal integration
	if (integ_lon)
		N_spatialSize = rows;		% Number of points in the spatial dim
		integDim = 2;							% Dimension along which we are going to integrate
		ini = fix(handles.head(3) / dlat) * dlat;
		fim = fix(handles.head(4) / dlat) * dlat + dlat;
		vecD = (ini:dlat:fim);
		Y = linspace(handles.head(3),handles.head(4), N_spatialSize);
	else
		N_spatialSize = cols;
		integDim = 1;
		ini = fix(handles.head(1) / dlat) * dlat;
		fim = fix(handles.head(2) / dlat) * dlat + dlat;
		vecD = (ini:dlat:fim);
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
	if (integ_lon),		N_tot = cols + 1e-10;		% Add eps so that we never have divisions by zero
	else				N_tot = rows + 1e-10;
	end
	mask = [];
	for (k = 1:nSeries)
		Z = nc_funs('varget', handles.fname, s.Dataset(z_id).Name, [k-1 0 0], [1 rows cols]);
		this_has_nans = false;

		% NaNify polygon exterior points?
		if (have_polygon && k == 1)
			mask = img_fun('roipoly_j',handles.head(1:2),handles.head(3:4),double(Z),x,y);
		end
		if (~isempty(mask)),	Z(~mask) = NaN;		end

		ind = isnan(Z);				% This may, or may not, be equal to 'mask'
		if (any(ind(:))),		Z(ind) = 0;		this_has_nans = true;		end

		tmp = sum(Z,integDim);				% Add along integration dim
		if (this_has_nans)
			tmp2 = sum(ind,integDim);		% Get total number of NaNs along interp dim
			tmp = tmp ./ (N_tot - tmp2);	% Now get the number of valid values along interp dim
		else
			tmp = tmp / N_tot;
		end
		% Now add all inside each stripe
		for (m = 1:nStripes)
			tmp2 = tmp( indStripe(m):indStripe(m+1) );
			tmp2(tmp2 == 0) = [];			% Not so unlikely
			allSeries(m,k) = sum(tmp2) / numel(tmp2);
		end

		h = aguentabar(k/nSeries);
		if (isnan(h)),	break,	end
	end
	if (isnan(h)),	return,		end	

	allSeries(allSeries == 0) = nan;		% NaN is more reasonable to denote data absence

	if ( nargin == 7 && have_polygon  && exist(fname2,'file') == 2 )
		out2 = zonal(handles, dlat, integ_lon, false, have_polygon, fname2);
		allSeries = double(out2) - allSeries;
	end
	allSeries = single(allSeries);

	if (~nargout)			% If no argout, show result in a Mirone/Ecran window
		zz = grdutils(allSeries,'-L');
		head = [1 nSeries vecD(1) vecD(end) zz(1) zz(2) 0 1 dlat];
		tmp.X = 1:nSeries;		tmp.Y = linspace( (vecD(1)+dlat/2), (vecD(end)-dlat/2), nStripes );
		if (~do_trends)		% 2D, Mirone
			tmp.head = [head(1:2) tmp.Y(1) tmp.Y(end) head(5:end)];
			tmp.geo = 0;			tmp.name = 'Zonal integration';
			mirone(allSeries, tmp)
		else				% 1D, Ecran
			trend = zeros(1,nStripes);
			for (k = 1:nStripes)
				p = polyfit(tmp.X, double(allSeries(k,:)), 1);
				trend(k) = p(1);
			end
			ind = find(~isnan(trend));			% Remove any eventual leading or trailing NsNs
			trend = trend(ind(1):ind(end) );
			tmp.Y = tmp.Y(ind(1):ind(end) );
 			ecran(handles, tmp.Y, trend, 'Slope of line fit')
		end
	else
		out = allSeries;
	end
	
% ----------------------------------------------------------------------
function calcGrad(handles)
% Calcula o gradiente de um fiche ja com as medias anuais por ajuste de um recta (Loada o fiche todo na memoria)

	[z_id, s, rows, cols] = get_ncInfos(handles);

	anos = handles.number_of_timesteps;
% 	anos = anos - 1;
	Tmed = zeros([rows, cols, anos]);	% Temp media para cada um dos anos
	for (m = 1:anos)
 		Tmed(:,:,m) = nc_funs('varget', handles.fname, s.Dataset(z_id).Name, [m-1 0 0], [1 rows cols]);
% 		Tmed(:,:,m) = nc_funs('varget', handles.fname, s.Dataset(z_id).Name, [m 0 0], [1 rows cols]);
	end

	aguentabar(0,'title','Compute the Time rate','CreateCancelBtn')

	Tvar = zeros(rows,cols);
	x = (0:anos-1)';
	for (m = 1:rows)
		for (n = 1:cols)
			y = squeeze(Tmed(m,n,:));
			ind = isnan(y);
			y(ind) = [];
			if (numel(y) < 10)				% Completely ad-hoc test
				Tvar(m,n) = NaN;
				continue
			end
%  			p = polyfit(x(~ind),y,1);
 			p = trend1d_m([x(~ind) y],'-L','-N2r');
%  			p = trend1d_m([x(~ind) y],'-L','-N2r','-R','-P');
% 			z=[xvalues(1:4);ones(1,4)]'\yvalues';
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

% ------------------------------------------------------------------------------
function calc_yearMean(handles, months, fname2, flag, pintAnoes, nCells, fname3, splina)
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
% FNAME3 	Optional name of a netCDF file where interpolated nodes will be set to FLAG
%			and the others retain their FNAME2 value. This corresponds to the promotion
%			of interpolated nodes to quality FLAG. Only used if PINTANOES == TRUE
% SPLINA	Logical that if true instruct to spline interpolate the missing monthly values
%			before computing the yearly mean. This option acumulates with that of PINTANOES

	txt1 = 'netCDF grid format (*.nc,*.grd)';	txt2 = 'Select output netCDF grid';
	[FileName,PathName] = put_or_get_file(handles,{'*.nc;*.grd',txt1; '*.*', 'All Files (*.*)'},txt2,'put','.nc');
	if isequal(FileName,0),		return,		end
	grd_out = [PathName FileName];

	[z_id, s, rows, cols] = get_ncInfos(handles);
	do_flags = false;
	track_filled = false;

	if (nargin == 1),		months = 1:12;		end		% Default to yearly means
	if (nargin >= 7 && ~isempty(fname3)),		track_filled = true;	end		% Qeep track of interpolated nodes

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
			errordlg('Buhhuu!! The quality flags file has less "planes" than the-to-be-flagged-file. By','Error'),	return
		end
		if (~isequal([rows cols], [s_flags.Dataset(z_id_flags).Size(end-1) s_flags.Dataset(z_id_flags).Size(end)]))
			errordlg('Buhhuu!! quality flags and the-to-be-flagged-file have not the same size. By','Error'),		return
		end
		do_flags = true;
		
		if (nargin == 3),		flag = 7;		end			% If not provided, defaults to best quality
	end

	if (nargin < 8),		splina = false;		end
	if (splina)
		n_pad_months = 7;		% Example: if months = 7:9 interpolation domain is 7-n_pad_months-1:9+n_pad_months
		total_months = numel(months) + 2*n_pad_months;
		ZtoSpline = alloc_mex(rows, cols, total_months, 'single', NaN);
	end

	handles.geog = 1;		handles.was_int16 = 0;		handles.computed_grid = 0;
	anos = handles.number_of_timesteps / 12;

	h = aguentabar(0,'title','Computing anual means.','CreateCancelBtn');
	Tmed = zeros([rows, cols]);			% Temp media para cada um dos anos
	in_break = false;					% Inner loop cancel option
	last_processed_month = 0;		already_processed = 0;
	warning off MATLAB:divideByZero
	%anos = 3;	%pintAnoes =0;

	for (m = 1:anos)
		if (splina)
% 			if (m == 1)
% 				this_months = 1:(months(end)+n_pad_months);		past_months = max(1,months(1)-n_pad_months) - 1;
% 			elseif (m == anos)
% 				this_months = 1:( n_pad_months + numel(months) + min(n_pad_months, 12-months(end)) );
% 				past_months = (m - 1)*12 + months(1) - n_pad_months - 1;
% 			else
% 				this_months = 1:(numel(months)+2*n_pad_months);		past_months = (m - 1)*12 + months(1) - n_pad_months - 1;
% 			end
			if (m == 1)
				past_months = max(1,months(1)-n_pad_months) - 1;
				this_months = past_months+1 : past_months+(months(end)+n_pad_months);
			elseif (m == anos)
				past_months = (m - 1)*12 + months(1) - n_pad_months - 1;
				this_months = past_months+1 : past_months+( n_pad_months + numel(months) + min(n_pad_months, 12-months(end)) );
			else
				past_months = (m - 1)*12 + months(1) - n_pad_months - 1;
				this_months = past_months+1 : past_months+(numel(months)+2*n_pad_months);
			end
		else
			this_months = months;				past_months = (m - 1)*12;
			contanoes = zeros(rows, cols);
		end

		counter = 0;
		for (n = this_months)
% 			mn = past_months + n - 1;
% 			Z = nc_funs('varget', handles.fname, s.Dataset(z_id).Name, [mn 0 0], [1 rows cols]);
			counter = counter + 1;
			if (m == 1)
				Z = nc_funs('varget', handles.fname, s.Dataset(z_id).Name, [n-1 0 0], [1 rows cols]);
			else
				if (n <= last_processed_month)
					already_processed = already_processed + 1;
					offset = 1;
					if (m == 2),	offset = months(1)-n_pad_months;		end		% Because for the 1st year the series may be shorter
					offset = total_months - (last_processed_month - n) + offset - 1;
					ZtoSpline(:,:,already_processed) = ZtoSpline(:,:,offset);
				else
					Z = nc_funs('varget', handles.fname, s.Dataset(z_id).Name, [n-1 0 0], [1 rows cols]);
					already_processed = 0;
				end
			end

% 			if (do_flags)
% 				Z_flags = nc_funs('varget', fname2, s_flags.Dataset(z_id_flags).Name, [mn 0 0], [1 rows cols]);
			if (do_flags && ~already_processed)
				Z_flags = nc_funs('varget', fname2, s_flags.Dataset(z_id_flags).Name, [n-1 0 0], [1 rows cols]);
				Z(Z_flags < flag) = NaN;
			end

			ind = isnan(Z);

 			if (pintAnoes && ~already_processed && any(ind(:)))			% If fill spatial holes is requested
%  			if (pintAnoes && any(ind(:)))			% If fill spatial holes is requested
				if (track_filled),		ind0 = ind;		end			% Get this Z level original NaNs mask
				Z = inpaint_nans(handles, Z, ind, nCells);			% Select interp method inside inpaint_nans()
				ind = isnan(Z);
				if (track_filled)
					Z_flags(ind0 & ~ind) = flag;	% Promote interpolated pixels to quality 'flag'
					% Write updated quality file
					if (mn == 0),		nc_io(fname3, sprintf('w%d/time',anos*numel(months)), handles, reshape(Z_flags,[1 size(Z_flags)]))
					else				nc_io(fname3, sprintf('w%d', mn), handles, Z_flags)
					end
				end
			end

			if (~splina)				% Do not interpolate along time (months)
				Z(ind) = 0;				% Transmutate the anoes
				contanoes = contanoes + ~ind;
				Tmed(:,:) = Tmed(:,:) + double(Z);
			else						% Pack this year into a 3D temporary variable, to be processed later.
%  				ZtoSpline(:,:,n) = Z;
 				if (~already_processed),	ZtoSpline(:,:,counter) = Z;		end
			end
		end			% End loop over months
		last_processed_month = this_months(end);

		if (~splina)				% Do not interpolate along time. Compute averages with all non NaNs
			Tmed(:,:) = Tmed(:,:) ./ contanoes;
			tmp = single(Tmed(:,:));

		else						% Fill missing month data by interpolation based on non-NaN data
			hh = aguentabar(eps,'title','Splining it.');
			if (isnan(hh)),		in_break = true;	break,		end		% Over time loop said: break
			n_meses = numel(this_months);

			if (m == 1),	first_wanted_month = months(1);				% First year in the stack
			else			first_wanted_month = n_pad_months + 1;
			end
			last_wanted_month = first_wanted_month + numel(months) - 1;
			
			for (i = 1:cols)
				for (j = 1:rows)
					if (already_processed),		break,	end
					y = double( squeeze( ZtoSpline(j,i,1:n_meses) ) );
					ind = ~isnan(y);
					% If have NaNs inside the months of interest and the overall series has enough points, interp in the missing positions
					if ( all( ind(first_wanted_month:last_wanted_month) ) )		% We have them all, so nothing to interp
						ZtoSpline(j,i,1:n_meses) = single(y);
					elseif ( numel(ind(ind)) >= (fix(numel(months)/2) + n_pad_months - 0) )
						x = this_months(ind);			y0 = y(ind);
						%yy = gmtmbgrid_m(x(:), zeros(numel(x),1)+2, y0(:), '-I1', sprintf('-R%d/%d/0/4', x(1), x(end)), '-Mz', '-W-2', '-T100');
						%yy = yy(3,:);
						try
							if ((m == 1 || m == anos))
% 								yy = interp1(x, y0, this_months, 'spline', nan);
%  								yy = cubicspline(x, y0, this_months);
 								yy = akimaspline(x, y0, this_months);
% 								yy = spline1d(this_months, x, y0, [], [], 0.0);
							else
% 								yy = spline(x, y0, this_months);
%  								yy = cubicspline(x, y0, this_months);
								yy = akimaspline(x, y0, this_months);
% 								yy = spline1d(this_months, x, y0, [], [], 0.0);
							end
							y(first_wanted_month:last_wanted_month)= yy(first_wanted_month:last_wanted_month);
							ZtoSpline(j,i,1:n_meses) = single(y);
						end
					else
						ZtoSpline(j,i,1:n_meses) = single((1:n_meses)*nan);
					end
				end
				hh = aguentabar(i/(cols+1));
				if (isnan(hh)),		in_break = true;	break,		end		% Over time loop said: break (comment: FCK ML SHIT)
			end				% end loops over this 2D layer

			if (in_break),		break,		end		% Fck no gotos paranoia obliges to this recursive break

			% Now we can finaly compute the season mean
			tmp = ZtoSpline(:,:,first_wanted_month);
			for (n = (first_wanted_month+1):last_wanted_month)
				tmp = cvlib_mex('add',tmp,ZtoSpline(:,:,n));
			end
			cvlib_mex('CvtScale',tmp,1/numel(months),0);
		end				% End interpolate along time

		if (in_break),		break,		end		% Fck no gotos paranoia obliges to this recursive break

		% Write this layer to file
		if (m == 1),		nc_io(grd_out, sprintf('w%d/time',anos), handles, reshape(tmp,[1 size(tmp)]))
		else				nc_io(grd_out, sprintf('w%d', m-1), handles, tmp)
		end

		h = aguentabar(m/anos,'title','Computing anual means.');
		if (isnan(h)),	break,	end

		if (~splina),	Tmed = Tmed * 0;	end			% Reset it to zeros
	end

% ----------------------------------------------------------------------
function pass_by_count(handles, count, fname2)
% Check the curently active 3D file against a count file
%
% COUNT		Threshold count value. Nodes on in-memory file that have a count on the
% 			corresponding node of FNAME2 < COUNT are set to NaN
%
% OPTIONS:
% FNAME2 	name of a netCDF file with the count quality flags. Asked if not provided
%			Obviously this file must be of the same size as the series under analysis.

	txt1 = 'netCDF grid format (*.nc,*.grd)';		txt2 = 'Select output netCDF grid';
	[FileName,PathName] = put_or_get_file(handles,{'*.nc;*.grd',txt1; '*.*', 'All Files (*.*)'},txt2,'put','.nc');
	if isequal(FileName,0),		return,		end
	grd_out = [PathName FileName];
	
	if (nargin == 2)		% No count grid transmitted. Ask for it
		[FileName,PathName] = put_or_get_file(handles,{'*.nc;*.grd',txt1; '*.*', 'All Files (*.*)'},'Select input netCDF file','get');
		if isequal(FileName,0),		return,		end
		fname2 = [PathName FileName];
	else					% Got a name. Check that it exists
		if (exist(fname2,'file') ~= 2)
			errordlg(['Blheak!! ' fname2 ' does not exist (even if you think so). Bye Bye'],'Error'),	return
		end
	end

	z_id = handles.netcdf_z_id;
	s = handles.nc_info;				% Retrieve the .nc info struct
	rows = s.Dataset(z_id).Size(end-1);
	cols = s.Dataset(z_id).Size(end);

	s_flags = nc_funs('info',fname2);
	[X,Y,Z,head,misc] = nc_io(fname2,'R');
	z_id_flags = misc.z_id;
	if ~(numel(head) == 9 && isfield(misc,'z_id'))
		errordlg(['Blheak!! ' fname2 ' is is not a file with presumably with a count of quality flags. By'],'Error'),	return
	end
	if (misc.z_dim(1) < handles.number_of_timesteps)
		errordlg('Buhhuu!! The count flags file has less "planes" than the-to-be-counted-file. By','Error'),	return
	end
	if (~isequal([rows cols], [s_flags.Dataset(z_id_flags).Size(end-1) s_flags.Dataset(z_id_flags).Size(end)]))
		errordlg('Buhhuu!! quality flags and the-to-be-counted-file have not the same size. By','Error'),		return
	end

	handles.geog = 1;		handles.was_int16 = 0;		handles.computed_grid = 0;
	
	aguentabar(0,'title',['NaNify countings < ' sprintf('%d',count)],'CreateCancelBtn')

	n_layers = handles.number_of_timesteps;
	for (m = 1:n_layers)

		Z = nc_funs('varget', handles.fname, s.Dataset(z_id).Name, [m-1 0 0], [1 rows cols]);
		Z_flags = nc_funs('varget', fname2, s_flags.Dataset(z_id_flags).Name, [m-1 0 0], [1 rows cols]);
		Z(Z_flags < count) = NaN;

		if (m == 1),		nc_io(grd_out, sprintf('w%d/time',n_layers), handles, reshape(Z,[1 size(Z)]))
		else				nc_io(grd_out, sprintf('w%d', m-1), handles, Z)
		end

		h = aguentabar(m/n_layers);
		if (isnan(h)),	break,	end
		
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
			Z_rect = cropimg(handles.head(1:2),handles.head(3:4),Z,rect_crop,'out_grid');
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
	[FileName,PathName] = put_or_get_file(handles,{'*.nc;*.grd',txt1; '*.*', 'All Files (*.*)'},txt2,'put','.nc');
	if isequal(FileName,0),		return,		end
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

% -----------------------------------------------------------------------------------------
function fid = write_vtk(handles)
% Write a 3D netCDF into a VTK format

	txt1 = 'VTK format (*.vtk)';
	txt2 = 'Select output VRT file';
	[FileName,PathName] = put_or_get_file(handles,{'*.vtk',txt1; '*.*', 'All Files (*.*)'},txt2,'put','.vtk');
	if isequal(FileName,0),		return,		end
	fname_out = [PathName FileName];

	s = handles.nc_info;				% Retrieve the .nc info struct
	rows = s.Dataset(handles.netcdf_z_id).Size(end-1);
	cols = s.Dataset(handles.netcdf_z_id).Size(end);
	nLayers = handles.number_of_timesteps;

	fid = fopen(fname_out, 'wb','b');
	fprintf(fid, '# vtk DataFile Version 2.0\n');
	fprintf(fid, 'converted from A B\n');
	fprintf(fid, 'BINARY\n');
	fprintf(fid, 'DATASET RECTILINEAR_GRID\n');
	fprintf(fid, 'DIMENSIONS %d %d %d\n', cols, rows, nLayers);
	fprintf(fid, 'X_COORDINATES %d float\n', cols);
	X = linspace(handles.head(1), handles.head(2), cols);
	fwrite(fid, X, 'real*4');
	fprintf(fid, 'Y_COORDINATES %d float\n', rows);
	X = linspace(handles.head(3), handles.head(4), rows);
	fwrite(fid, X, 'real*4');
	fprintf(fid, 'Z_COORDINATES %d float\n', nLayers);
	fwrite(fid, 1:nLayers, 'real*4');
	fprintf(fid, 'POINT_DATA %d\n', cols * rows * nLayers);
	fprintf(fid, 'SCALARS dono float 1\n');
	fprintf(fid, 'LOOKUP_TABLE default\n');

	for (m = 1:nLayers)
		Z = nc_funs('varget', handles.fname, s.Dataset(handles.netcdf_z_id).Name, [m-1 0 0], [1 rows cols]);
		Z = double(Z');
		fwrite(fid, Z(:), 'real*4');
	end
		
% ----------------------------------------------------------------------
function [z_id, s, rows, cols] = get_ncInfos(handles)
% Since this is done in nearly all functions, centralize it here
	z_id = handles.netcdf_z_id;
	s = handles.nc_info;			% Retrieve the .nc info struct
	rows = s.Dataset(z_id).Size(end-1);
	cols = s.Dataset(z_id).Size(end);

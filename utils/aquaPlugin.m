function aquaPlugin(handles, auto)
% Plugin function that is called by Aquamoto. Use this function to write custom code
% to solve particular problems taking advantage from the fact that a LOT of information
% about the netCDF files is stored in HANDLES.
%
% OPTIONS
%
% AUTO	- If logical and TRUE, search the OPTcontrol.txt file for the MIR_AQUAPLUG key that
%		  should point into a file name of a control script.
%		- If it is a string, than that is interpreted as the name of the control script.
%
%		A "control script" is a file with the EXACT arguments to select and run one of
%		main functions here as pointed by the CASOS cell array below.
%
%		One way of executing this functionality is to check the "Seek OPTcontrol.txt" checkbox
%		In which case the OPTcontrol.txt will be scanned for the name of the control script.
%		This works both for the ML and the standalone version.
%		The other way, restricted to the ML version, is to run in the Matlab command line:
%				aquamoto file.nc 'file_name_of_control_script'
%		OR
%				aquamoto('file.nc', 0)
%		In the later case the control script name is searched in the OPTcontrol.txt file

%	Copyright (c) 2004-2020 by J. Luis
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

% $Id: aquaPlugin.m 11412 2019-03-05 19:46:53Z j $

	if (isempty(handles.fname))
		errordlg('Fast trigger, you probably killed my previous encarnation. Now you have to start again. Bye.','Error')
		return
	end
	internal_master = true;		% To know if flow control is determined by the contents of an external file (def NO).

	casos = {'zonal' ...			% 1 - Compute zonal means
			'tvar' ...				% 2 - Compute the Temp time rate of a file with annual means by fit of a straight line (Load entire file in memory)
			'applyFlags' ...		% 3 - Check against its 3D flags companion and replace values < FLAG to NaN
			'yearMeanFlag' ...		% 4 - Compute yearly averages from monthly data but checked against a quality flag file
			'polygAVG' ...			% 5 - Compute averages of whatever inside polygons (if any)
			'flagsStats' ...		% 6 - Compute per/pixel annual or month counts of pixel values with a quality >= flag
			'pass_by_count' ...		% 7 - Check the curently active 3D file against a count file
			'do_math' ...			% 8 - Perform some basic algebric operations with the 3D planes
			'conv2vtk' ...			% 9 - Convert a 3D netCDF file into a VTK format
			'L2_periods' ...		% 10 - Calculate composites of L2 products over fixed periods
			'corrcoef' ...			% 11 - Calculate correlation coefficients between 2 3D arrays
			'count_blooms' ...		% 12 - Count chlor_a blooms
			};

	qual = casos{11};			% <== Active by MANUAL selection. May be override by next section

	n_args = nargin;
	if (isfield(handles,'check_plugFun') && get(handles.check_plugFun, 'Val'))	% This way, the stand alone version can work too
		n_args = 2;		auto = true;
	end

	if (n_args == 2)				% Go figure out if we have a controlling script
		out = script_control(handles, auto);
		if (~isempty(out))
			qual = casos{out{1}};	% The rest will be applied blindly. If it screws, screws
			internal_master = false;
		else
			errordlg('You directed AQUAPLUGIN to work on script mode but OPTcontrol.txt misses necessary info','Error')
			return
		end
	end

	switch qual
		case 'zonal'				% CASE 1
			integ_lon = true;
			dlat = 1.0;
			trends = false;			% If true compute the trends (per stripe) of the zonal integration
			sub_set = [0 0];		% [jump_start stop_before_end], make it [] or [0 0] to be ignored
			fnamPoly1 = 'C:\a1\pathfinder\plataforma_poly.dat';		% If this name is uncorrect, another will be asked
			fnamPoly2 = 'poly_largo.dat';	%fnamPoly2= [];	% If it exists, compute difference of zonal integrations
			fnameFlag  = 'C:\a1\pathfinder\qual_82_09.nc';	% If not empty check againts this file (For monthly data)
			quality = 6;			% Retain only values of quality >= this (or <= abs(this) when MODIS). Ingored if fname = []
			if (internal_master)
				zonal(handles, dlat, integ_lon, trends, sub_set, fnamPoly1, fnamPoly2, fnameFlag, quality)
			else
				zonal(handles, out{2:end})
			end
		case 'tvar'					% CASE 2
			slope = true;			% TRUE to compute slope of linear fit, FALSE to compute p-value parameter
			sub_set = [0 0];		% [jump_start stop_before_end], make it [] or [0 0] to be ignored
			fname  = 'C:\a1\pathfinder\qual_82_09.nc';	% If not empty check againts this file (For monthly data)
			quality = 6;			% Retain only values of quality >= this (or <= abs(this) when MODIS). Ingored if fname = []
			splina = false;			% Fill missing monthly data by a spline interpolation. Ignored if fname = [].
			scale = 12;				% Scale rate of change by this value (useful when input data has monthly data).
			if (internal_master)
				calcGrad(handles, slope, sub_set, fname, quality, splina, scale)
			else
				calcGrad(handles, out{2:end})
			end
		case 'applyFlags'			% CASE 3
			fnameFlag  = 'C:\a1\pathfinder\qual_82_09.nc';
			flag = 6;				% Compute yearly means
			if (internal_master),	applyFlags(handles, fnameFlag, flag, 200)
			else,					applyFlags(handles, out{2:end})
			end
		case 'yearMeanFlag'			% CASE 4
			ano = 1:12;				% Compute yearly (ano = 1:12) or seasonal means (ano = start_month:end_month)
			fname  = 'C:\a1\MODIS\algas\Algas_qual_nsst_TERRA_00_10.nc';
			%fname  = 'C:\a1\pathfinder\qual_82_09.nc';
			quality = 0;			% Retain only values of quality >= this (or <= abs(this) when MODIS)
			nCells = 200;			% Holes (bad data) smaller than this are filled by interpolation
			% Where to save track of filled holes. Ignored if nCells = 0 OR fname3 = []
			fname3 = [];	%'C:\a1\pathfinder\qual7_85_07_Interp200_Q6.nc';
			%splina = true;
			splina = [12 30];		% Fill missing monthly data by a spline interpolation taken over two years (out limits set to NaN)
			tipoStat = 0;			% 0, Compute MEAN, 1 -> Median; 2 -> MINimum; 3 -> MAXimum; 4 -> STD of the ANO period
			% If not empty, it must contain the name of a Lon,Lat file with locations where to output time series
			chkPts_file = [];	%chkPts_file = 'C:\a1\pathfinder\chkPts.dat';
			if (internal_master)
				calc_yearMean(handles, ano, fname, quality, nCells, fname3, splina, tipoStat, chkPts_file)
			else
				calc_yearMean(handles, out{2:end})
			end
		case 'polygAVG'				% CASE 5
			fnameOut = [];			% If not empty, file name where to save the result (otherwise, asked at the end)
			op = [];				% Type average (or other). [] means doing average. 'median' will do a median
			fnamePolys = [];		% If not empty, file name of polygon or list of polygons (otherwise, fished from Mirone fig)
			sub_set = [0 0];		% [jump_start stop_before_end], make it [] or [0 0] to be ignored
			fnameFlag  = 'C:\a1\pathfinder\qual_82_09.nc';	% If not empty check againts this file
			quality = 6;			% Retain only values of quality >= this (or <= abs(this) when MODIS). Ingored if fnameFlag = []
			if (internal_master),	calc_polygAVG(handles, fnameOut, op, fnamePolys, sub_set, fnameFlag, quality)
			else,					calc_polygAVG(handles, out{2:end})
			end
		case 'flagsStats'			% CASE 6
			ano = 1:12;				% Compute yearly stats
			%opt = '';				% Make the counting on a per month basis
			opt = 'per_year';		% Make the counting on a per year basis
			if (internal_master),	calc_flagsStats(handles, ano, 7, opt)
			else,					calc_flagsStats(handles, out{2:end})
			end
		case 'pass_by_count'		% CASE 7
			count = 11;
			fname = 'C:\a1\pathfinder\countPerYear_flag7_Interp200.nc';
			if (internal_master),	pass_by_count(handles, count, fname)
			else,					pass_by_count(handles, out{2:end})
			end
		case 'do_math'				% CASE 8
			opt = 'diffstd';		% Sum all layers
			if (internal_master)
				do_math(handles, opt, [19 0], 'C:\a1\MODIS\mediaAnual_TERRA_NSST_Interp200_Q0.nc', [1 1])
			else
				do_math(handles, out{2:end})
			end
		case 'conv2vtk'				% CASE 9
			write_vtk(handles)
		case 'L2_periods'			% CASE 10
			period = 3;				% Number of days in each period
			regMinMax = [0 inf];	% If want to limit the admissible Z values ([min max])
			tipoStat = 0;			% 0, Compute MEAN, 1 compute MINimum and 2 compute MAXimum of the ANO period
			grd_out = 'C:\SVN\mironeWC\tmp\lixoL2.nc';
			if (internal_master)
				calc_L2_periods(handles, period, tipoStat, regMinMax, grd_out)
			else
				calc_L2_periods(handles, out{2:end})
			end
		case 'corrcoef'				% CASE 11
			secondArray = 'C:\a1\pathfinder\lixoCloros.nc';
			grd_out = 'C:\SVN\mironeWC\tmp\lixoR.nc';
			sub_set = [0 0];		% [jump_start stop_before_end], make it [] or [0 0] to be ignored
			if (internal_master),	calc_corrcoef(handles, secondArray, sub_set, false, grd_out)
			else,					calc_corrcoef(handles, out{2:end})
			end
		case 'count_blooms'			% CASE 12
			if (internal_master),	count_blooms(handles)
			else,					count_blooms(handles, out{2:end})
			end
	end

% ----------------------1-------2--------3---------4---------5---------6-----------7----------8---------9----
function out = zonal(handles, dlat, integ_lon, do_trends, sub_set, fnamePoly1, fnamePoly2, fnameFlag, quality)
% Compute zonal means from a multi-layer file
%
% DLAT 			width of the box in the direction orthogonal to INTEG_LON
% INTEG_LON 	(LOGICAL) If true, integration is done along longitude
% DO_TRENDS		(LOGICAL) If false compute zonal integrations, otherwise compute trends of the zonal integrations (per DLAT)
%				Note that X coords represent different things in the above two cases:
%				- Layer number for the zonal integration case and is up tp the user to make it correspond to a date.
%				- Spatial coordinate orthogonal to the integration direction for the DO_TRENDS == TRUE case
%
% OPTIONS: Since the number of options is variable some make mandatory that prev args exist. In that case use [] if needed
%
% SUB_SET	->  A two columns row vec with number of the offset of years where analysis start and stop.
%				For example [3 1] Starts analysis on forth year and stops on the before last year.
%				[0 0] Means using the all dataset.
%
% FNAMEPOLY1	Optional polygon file delimiting an area where the analysis will be carried on.
%
% FNAMEPOLY2	Optional second polygon file. If it points to a valid file. This function is called 
%				twice and results are subtracted
%
% FNAMEFLAG		name of a netCDF file with quality flags. Obviously this file must be of
%				the same size as the series under analysis. If not provided no quality check is done.
%
% QUALITY		Threshold quality value. Only values of quality >= FLAG will be taken into account
%				NOTE: For MODIS use negative FLAG. Than, values are retained if quality <= abs(FLAG)

	if (nargin < 4)
		errordlg('ZONAL: called with less than ninimum number of arguments', 'Error'),	return
	end
	if (nargin == 4)
		sub_set = [0 0];		fnamePoly1 = [];	fnamePoly2 = [];	fnameFlag = [];		quality = 7;
	elseif (nargin == 5)
		fnamePoly1 = [];		fnamePoly2 = [];	fnameFlag = [];		quality = 7;
	elseif (nargin == 6)
		fnamePoly2 = [];		fnameFlag = [];		quality = 7;
	elseif (nargin == 7)
		fnameFlag = [];			quality = 7;	% This 'quality' def is only to not error in one case.
	end

	if (numel(sub_set) == 2)
		jump_start = sub_set(1);		stop_before_end = sub_set(2);
	else
		jump_start = 0;					stop_before_end = 0;
	end

	if (~isempty(fnamePoly1))
		if (exist(fnamePoly1,'file') ~= 2)		% If given name does not exist, give another chance
			[FileName,PathName] = put_or_get_file(handles, ...
				{'*.dat;*.DAT', 'Data files (*.dat)';'*.*', 'All Files (*.*)'},'Enter polygon file','get');
			if (isequal(FileName,0)),		return,		end
			fnamePoly1 = [PathName FileName];
		end
		S = load(fnamePoly1);
		x = S(:,1);		y = S(:,2);
	else
		x = [];			y = [];
	end

	[z_id, s, rows, cols] = get_ncInfos(handles);

	%------------- Check for quality flags request -------------------
	if (~isempty(fnameFlag))
		[s_flags, z_id_flags, msg] = checkFlags_compat(fnameFlag, handles.number_of_timesteps, rows, cols);
		if (~isempty(msg)),	errordlg(msg, 'Error'),		return,		end
		do_flags = true;
		if (quality > 0 ),	growing_flag = true;		% PATHFINDER flags
		else,				growing_flag = false;		% MODIS flags
		end
	else
		do_flags = false;
	end

	% ------------- Build the vectors to deal with the zonal integration -------------
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

	nSeries = handles.number_of_timesteps - (jump_start + stop_before_end);	% Number of layers to be used in this run
	series_vec = (jump_start:(nSeries - 1 + jump_start)) + 1;		% Add 1 so it never starts at 0 (no good for indices)
	allSeries = zeros(nStripes, nSeries);
	if (integ_lon),		N_tot = cols + 1e-10;		% Add eps so that we never have divisions by zero
	else,				N_tot = rows + 1e-10;
	end
	mask = [];

	for (k = series_vec)			% Loop over time layers
		Z = nc_funs('varget', handles.fname, s.Dataset(z_id).Name, [k-1 0 0], [1 rows cols]);
		this_has_nans = false;

		if (do_flags)
			flags = nc_funs('varget', fnameFlag, s_flags.Dataset(z_id_flags).Name, [k-1 0 0], [1 rows cols]);
			if (growing_flag),		Z(flags < quality) = NaN;	% Pathfinder style (higher the best) quality flag
			else,					Z(flags > quality) = NaN;	% MODIS style (lower the best) quality flag
			end
		end

		% NaNify polygon exterior points?
		if (~isempty(fnamePoly1) && k == series_vec(1))
			mask = img_fun('roipoly_j',handles.head(1:2),handles.head(3:4),double(Z),x,y);
		end
		if (~isempty(mask)),	Z(~mask) = NaN;		end

		ind = isnan(Z);						% This may, or may not, be equal to 'mask'
		if (any(ind(:))),		Z(ind) = 0;		this_has_nans = true;		end

		tmp = sum(Z, integDim);				% Add along integration dim
		if (this_has_nans)
			tmp2 = sum(ind,integDim);		% Get total number of NaNs along interp dim
			tmp = tmp ./ (N_tot - tmp2);	% Now get the number of valid values along interp dim
		else
			tmp = tmp / N_tot;
		end
		% Now add all inside each stripe
		for (m = 1:nStripes)
			tmp2 = tmp(indStripe(m):indStripe(m+1));
			tmp2(tmp2 == 0) = [];			% Not so unlikely
			allSeries(m,k) = sum(tmp2) / numel(tmp2);
		end

		h = aguentabar(k/nSeries);
		if (isnan(h)),	break,	end
	end

	if (isnan(h)),	return,		end	

	allSeries(allSeries == 0) = nan;		% NaN is more reasonable to denote data absence

	if (~isempty(fnamePoly2) && exist(fnamePoly2,'file') == 2)
		out2 = zonal(handles, dlat, integ_lon, false, sub_set, fnamePoly2, [], fnameFlag, quality);
		allSeries = double(out2) - allSeries;
	end
	allSeries = single(allSeries);

	% ------------ If no argout, show result in a Mirone/Ecran window ------------
	if (~nargout)
		zz = grdutils(allSeries,'-L');
		head = [1 nSeries vecD(1) vecD(end) zz(1) zz(2) 0 1 dlat];
		clear tmp;				% To shut up a useless ML warning due to what will happen at next line
		tmp.X = 1:nSeries;		tmp.Y = linspace((vecD(1)+dlat/2), (vecD(end)-dlat/2), nStripes);
		if (~do_trends)			% 2D, Mirone
			tmp.head = [head(1:2) tmp.Y(1) tmp.Y(end) head(5:end)];
			tmp.geo = 0;		tmp.name = 'Zonal integration';
			mirone(allSeries, tmp)
		else					% 1D, Ecran
			trend = zeros(1,nStripes);
			for (k = 1:nStripes)
				p = polyfit(tmp.X, double(allSeries(k,:)), 1);
				trend(k) = p(1);
			end
			ind = find(~isnan(trend));			% Remove any eventual leading or trailing NaNs
			trend = trend(ind(1):ind(end));
			tmp.Y = tmp.Y(ind(1):ind(end));
 			ecran(handles, tmp.Y, trend, 'Slope of line fit')
		end
	else
		out = allSeries;
	end
	
% -------------------1-------2-------3---------4---------5-------6-------7-------8-------9--------10-----
function calcGrad(handles, slope, sub_set, fnameFlag, quality, splina, scale, grd_out, do_3x3, mask_file)
% Compute the rate of change of a file by fitting a LS straight line. The file can be one of the
% already computed yearly means, in which case last three input arguments do not apply,
% OR the full time series. In this case optional checking against quality flags and spline
% interpolation can be done using info transmitted via the last three args.
%
% NOTE: THIS IS A HIGHLY MEMORY CONSUMPTION ROUTINE AS ALL DATA IS LOADED IN MEMORY
%
% SLOPE		Logical indicating if compute slope of linear fit (TRUE) or the p parameter (FALSE)
%
% SUB_SET -> A two elements row vec with number of the offset of years where analysis start and stop.
%			For example [3 1] Starts analysis on forth year and stops on the before last year.
%			[0 0] Means using the all dataset.
%
% OPTIONS:
% FNAMEFLAG	name of a netCDF file with quality flags. Obviously this file must be of
% 			the same size as the series under analysis.
%
% QUALITY	Threshold quality value. Only values of quality >= FLAG will be taken into account
%			NOTE: For MODIS use negative FLAG. Than, values are retained if quality <= abs(FLAG)
%
% SPLINA	Logical that if true instruct to spline interpolate the missing monthly values
%			before computing the rate of change.
%
% SCALE		Scale the final rate by this value. The idea of all this is that input data
%			can have monthly means and we want to compute time rate of change per year.
%			In this case use SCALE = 12.
%			Ignored if SLOPE is FALSE
%
% GRD_OUT	Name of the netCDF file where to store the result. If not provided, open Mirone Fig.
%
% DO_3x3	Logical. If true, compute over 3x3 windows
%
% MASK_FILE	name of Land mask file. If it has NaNs it will be multiplied but if only 1/0's, 0's will become NaNs

	if (nargin <= 6),	scale = 1;		end
	if (nargin < 8),	grd_out = [];	end
	if (nargin < 9),	do_3x3 = false;	end
	if (nargin < 10),	mask_file = '';	end
	if (~slope),		scale = 1;		end		% Make sure to not scale p-values

	do_blockMean = false;	% Must be turned into an input parameter
	do_flags = false;		% Will be set to true if we do a checking against a quality flgas file
	get_profiles_in_polygon = false;			% Save all profiles (along third dim) located inside the polygonal area
	n_anos = handles.number_of_timesteps;
	[z_id, s, rows, cols] = get_ncInfos(handles);

	% Find if we are dealing with a Pathfinder V5.2 daily file
	[is_PFV52, tempos] = PFV52(handles, s);		% if not a PFV5.2, tempos = handles.time(:)

	if (nargin >= 3 && (numel(sub_set) == 2))
		jump_anos = sub_set(1);		stop_before_end_anos = sub_set(2);
	else
		jump_anos = 0;				stop_before_end_anos = 0;
	end

	% For PFV5.2 daily data, accept the start-stop as time in years. e.g. [1986 2011]
	if (is_PFV52 && jump_anos >= 1982 && stop_before_end_anos <= 2014)	% Though we only have them till 2011
		ind = find(tempos > jump_anos);
		jump_anos = ind(1);
		ind = find(tempos > stop_before_end_anos);
		stop_before_end_anos = numel(tempos) - ind(1);
	end

	n_anos = n_anos - (jump_anos + stop_before_end_anos);	% Number of layers to be used in this run

	% When files are too big (not difficult when number of layers is high) must do a patch processing
	one_Gb = 1024*1024*1024;		% Use patches of 1 Gb
	if (rows * cols * n_anos * 4 > one_Gb)
		slice_cols = round(one_Gb / (rows * n_anos * 4));
		slice_cols = 1:slice_cols:cols;
		if (slice_cols(end) < cols),	slice_cols(end+1) = cols;	end
		slicing = true;
	else
		slice_cols = [1 cols];
		slicing = false;
	end
	read_cols = cols;

	if (nargin >= 4 && ~isempty(fnameFlag))
		if (slicing && ~handles.IamCompiled)
			slice_cols = [1 cols];
			slicing = false;
			warndlg('File very big that would be processed by patches, but that is not implemented whid Flags. Trying wirh a single patch.', 'Warning')
		elseif (slicing && handles.IamCompiled)
			errordlg('File too big. Processing it with quality flags is not implemented','Error'),	return
		end
		[s_flags, z_id_flags, msg] = checkFlags_compat(fnameFlag, handles.number_of_timesteps, rows, cols);
		if (~isempty(msg)),	errordlg(msg, 'Error'),		return,		end
		flags = alloc_mex(rows, cols, n_anos, 'uint8');
		for (m = 1:n_anos)
			flags(:,:,m) = nc_funs('varget', fnameFlag, s_flags.Dataset(z_id_flags).Name, [(m - 1 + jump_anos) 0 0], [1 rows cols]);
		end
		do_flags = true;
		if (nargin == 4)			% Default to Pathfinder max quality
			quality = 7;	splina = false;
		elseif (nargin == 5)
			splina = false;
		end
		if (quality > 0 ),	growing_flag = true;		% PATHFINDER flags
		else,				growing_flag = false;		% MODIS flags
		end
	else
		splina = false;		flags = [];	growing_flag = false;	% Not used
	end

	if (is_PFV52)			% For Pathfinder V5.2 (daily) use the true time coordinates
		x = tempos;		yy = [];
		if (splina),	yy = tempos;			end
		threshold = 0.10;	% AD-HOC threshold percentage of the n points below which do NOT compute slope
	else
		x = (0:n_anos-1)';	yy = [];
		if (splina),	yy = (0:n_anos-1)';		end
		threshold = 0.66;	% 66%
	end

try
	nSlices = numel(slice_cols) - 1;
	for (ns = 1:nSlices)		% Loop over number of slices (which may be only one)
		msg = 'Loading data ...';
		if (slicing)
			read_cols = slice_cols(ns+1) - slice_cols(ns) + 1;
			if (ns == 1)
				Tvar = zeros(rows, cols);	% To gether the final result. Allocate on first usage
			end
			msg = sprintf('Loading (big) data. %d of %d ...', ns, nSlices);
		end
		Tmed = alloc_mex(rows, read_cols, n_anos, 'single');
		aguentabar(0,'title', msg)
		if (~slicing)
			for (m = 1:n_anos)
		 		Tmed(:,:,m) = nc_funs('varget', handles.fname, s.Dataset(z_id).Name, [(m - 1 + jump_anos) 0 0], [1 rows cols]);
			end
		else
			for (m = 1:n_anos)
		 		t = nc_funs('varget', handles.fname, s.Dataset(z_id).Name, [(m - 1 + jump_anos) 0 0], [1 rows cols]);
				Tmed(:, :, m) = t(:, slice_cols(ns):slice_cols(ns+1));
			end
			clear t
		end

		% ---- save profiles of points, located inside polygon of Mirone fig, as a multi-segment file
		if (get_profiles_in_polygon)
			profiles_in_polygon(handles, Tmed, n_anos)
			return
		end
		% -----------------------------------------------------------------------------------------

		if (do_blockMean)
			aguentabar(0,'title','Compute block means','CreateCancelBtn')
			for (m = 1:n_anos)
				if (do_flags)
					slice = Tmed(:,:,m);
					fslice = flags(:,:,m);
					if (growing_flag),	slice(fslice < quality) = NaN;
					else,				slice(fslice > quality) = NaN;
					end
					Tmed(:,:,m) = mirblock(slice, '-A3');
				else
					Tmed(:,:,m) = mirblock(Tmed(:,:,m), '-A3');
				end
				h = aguentabar(m/n_anos);
				if (isnan(h)),	return,		end
			end
			if (do_flags)
				clear slice fslice
			end
		end

		if (is_PFV52)
			aguentabar(0,'title','Computing and removing Seazonal cycle','CreateCancelBtn')
			Tavg = tideman(handles, Tmed, tempos, 52);
			Tmed = remove_seazon(handles, Tmed, Tavg, tempos);
			clear Tavg
		end

		if (~slicing)			%Do it all in one passage
			if (do_3x3)
				[Tvar, stopit] = get_slopes_conn8(Tmed, flags, x, rows, cols, slope, do_flags, quality, growing_flag, threshold);
			else
				[Tvar, stopit] = get_slopes(Tmed, flags, x, yy, rows, cols, slope, do_flags, quality, growing_flag, splina, threshold);
			end
		else
			if (do_3x3)
				[Tvar_slice, stopit] = get_slopes_conn8(Tmed, flags, x, rows, read_cols, slope, do_flags, quality, growing_flag, threshold);
			else
				[Tvar_slice, stopit] = get_slopes(Tmed, flags, x, yy, rows, read_cols, slope, do_flags, quality, growing_flag, splina, threshold);
			end
			if (stopit),	return,		end
			clear Tmed
			Tvar(:, slice_cols(ns):slice_cols(ns+1)) = Tvar_slice;
		end
	end
catch
	disp(lasterror)
end

	clear Tmed
	if (stopit),	return,		end

	if (scale ~= 1),		cvlib_mex('CvtScale', Tvar, double(scale),0);		end
	Tvar = single(Tvar);

	if (~isempty(mask_file))
		Tvar = apply_mask(mask_file, [], Tvar);		% Apply the Land mask file
	end

	zz = grdutils(Tvar,'-L');  handles.head(5:6) = [zz(1) zz(2)];
	tmp.head = handles.head;
	if (isempty(grd_out))			% Show result in a Mirone figure
		tmp.X = linspace(tmp.head(1),tmp.head(2),cols);
		tmp.Y = linspace(tmp.head(3),tmp.head(4),rows);
		tmp.name = 'Time gradient (deg/year)';
		mirone(Tvar, tmp)
	else							% Got output name from input arg
		handles.was_int16 = false;
		handles.computed_grid = true;
		handles.geog = 1;
		nc_io(grd_out, 'w', handles, Tvar)
	end

% --------------------------------------------------------------------------------------
function [Tvar, stopit] = get_slopes(Tmed, flags, x, yy, rows, cols, slope, do_flags, quality, growing_flag, splina, threshold)
% Compute the best fit slopes (or p-values) on each node
	Tvar = zeros(rows, cols) * NaN;
	stopit = false;
	n_anos = numel(x);
	atLeast = round(n_anos * threshold);
	aguentabar(0,'title','Computing the Time rate','CreateCancelBtn')

	try
	for (n = 1:cols)
		for (m = 1:rows)
			y = double(squeeze(Tmed(m,n,:)));
			if (do_flags)
				this_flag = squeeze(flags(m,n,:));
				if (growing_flag),		y(this_flag < quality) = NaN;	% Pathfinder style (higher the best) quality flag
				else,					y(this_flag > quality) = NaN;	% MODIS style (lower the best) quality flag
				end
			end
			ind = isnan(y);
			y(ind) = [];
			if (numel(y) < atLeast),	continue,	end	% Completely ad-hoc test (it also jumps land cells)

			if (splina)
				if (~all(ind))		% Otherwise we have them all and so nothing to interp
					akimaspline(x(~ind), y, x, yy);
					y = yy;										
					if (ind(1) || ind(end))				% Cases when where we would have extrapolations
						if (ind(1))
							ki = 1;
							while (ind(ki)),	ki = ki + 1;	end
							y(1:ki) = NaN;				% Reset extraped values to NaN
						end
						if (ind(end))
							kf = numel(ind);
							while (ind(kf)),	kf = kf - 1;	end
							y(kf:numel(ind)) = NaN;		% Reset extraped values to NaN
						end
						ind = isnan(y);					% Get new nan indices again
						y(ind) = [];
					else
						ind = false(n_anos,1);			% Pretend no NaNs for the rest of the code below
					end
				end
			end

 			%p = polyfit(x(~ind),y,1);
			%z=[xvalues(1:4);ones(1,4)]'\yvalues';
			if (slope)		% Compute sople of linear fit
				p = trend1d_m([x(~ind) y],'-L','-N2r');
				Tvar(m,n) = p(1);
			else			% Compute p value
	 			p = trend1d_m([x(~ind) y],'-L','-N2r','-R','-P');
				if (p(1) < -0.5 || p(1) > 1),	continue,	end		% Another ad-hoc (CLIPPING)
				Tvar(m,n) = p(4);
			end
		end
		h = aguentabar(n/cols);
		if (isnan(h)),	stopit = true;	break,	end
	end
	aguentabar(1)
	catch
		disp(lasterror)
	end

% --------------------------------------------------------------------------------------
function [Tvar, stopit] = get_slopes_conn8(Tmed, flags, x, rows, cols, slope, do_flags, quality, growing_flag, threshold)
% Compute the best fit slopes (or p-values) centered on each node but using a 8-connection.
% That is, on each node the value will be an average of all the 8 data points arround, plut itself.
	Tvar = zeros(rows, cols) * NaN;
	stopit = false;
	aguentabar(0,'title','Computing the Time rate','CreateCancelBtn')

	% First do the left and right columns with the connection 1 algo. Won't do the Top & Bottom rows
	atLeast = round(numel(x) * threshold);		% Minimum number of points needed to do the fit
	for (n = [1 cols])
		for (m = 2:rows-1)
			y = double(Tmed(m,n,:));
			y = y(:);
			if (do_flags)
				this_flag = flags(m,n,:);
				if (growing_flag),		y(this_flag < quality) = NaN;	% Pathfinder style (higher the best) quality flag
				else,					y(this_flag > quality) = NaN;	% MODIS style (lower the best) quality flag
				end
			end
			ind = isnan(y);
			y(ind) = [];
			if (numel(y) < atLeast),	continue,	end				% Completely ad-hoc test (it also jumps land cells)

			if (slope)		% Compute sople of linear fit
				p = trend1d_m([x(~ind) y],'-L','-N2r');
				Tvar(m,n) = p(1);
			else			% Compute p value
	 			p = trend1d_m([x(~ind) y],'-L','-N2r','-R','-P');
				if (p(1) < -0.5 || p(1) > 1),	continue,	end		% Another ad-hoc (CLIPPING)
				Tvar(m,n) = p(4);
			end
		end
	end
	
	atLeast = round(numel(x) * threshold * 9);		% 9 because we are doing  a 3x3 window
	x = repmat(x, 1, 9)';
	x = x(:);
	for (n = 2:cols-1)
		for (m = 2:rows-1)
			y = double(Tmed([m-1 m m+1],[n-1 n n+1],:));
			y = y(:);
			if (do_flags)
				this_flag =flags([m-1 m m+1],[n-1 n n+1],:);
				if (growing_flag),		y(this_flag < quality) = NaN;	% Pathfinder style (higher the best) quality flag
				else,					y(this_flag > quality) = NaN;	% MODIS style (lower the best) quality flag
				end
			end
			ind = isnan(y);
			y(ind) = [];
			if (numel(y) < atLeast),	continue,	end				% Completely ad-hoc test (it also jumps land cells)

			if (slope)		% Compute sople of linear fit
				p = trend1d_m([x(~ind) y],'-L','-N2r');
				Tvar(m,n) = p(1);
			else			% Compute p value
	 			p = trend1d_m([x(~ind) y],'-L','-N2r','-R','-P');
				if (p(1) < -0.5 || p(1) > 1),	continue,	end		% Another ad-hoc (CLIPPING)
				Tvar(m,n) = p(4);
			end
		end
		h = aguentabar(n/cols);
		if (isnan(h)),	stopit = true;	break,	end
	end
	aguentabar(1)

% --------------------------------------------------------------------------------------
function profiles_in_polygon(handles, Tmed, n_anos)
% Save profiles of points, located inside polygon of Mirone fig, as a multi-segment file

	hFigs = findobj(0,'type','figure');						% Fish all figures
	IAmAMir = zeros(1, numel(hFigs));
	for (k = 1:numel(hFigs))								% Get the first Mirone figure with something in it
		if (~isempty(getappdata(hFigs(k), 'IAmAMirone')))
			handMir = guidata(hFigs(k));
			if (handMir.no_file),	continue,	end			% A virgin Mirone bar figure
			IAmAMir(k) = 1;		break,	
		end
	end
	if (sum(IAmAMir) ~= 1)
		errordlg('Did not find any valid Mirone figure with data displayed.','Error'),	return
	end

	hLine = findobj(handMir.axes1,'Type','line');
	if (isempty(hLine)),	hLine = findobj(handMir.axes1,'Type','patch');	end		% Try once more
	x = get(hLine,'XData');		y = get(hLine,'YData');
	if (isempty(x))
		errordlg('The Mirone figure needs to have at least one polygon loaded.','Error'),	return
	end
	mask = img_fun('roipoly_j',handles.head(1:2),handles.head(3:4),get(handMir.hImg,'CData'),x,y);
	B = img_fun('find_holes',mask);
	col_min = min(B{1}(:,2));		col_max = max(B{1}(:,2));
	row_min = min(B{1}(:,1));		row_max = max(B{1}(:,1));
	row_vec = row_min:row_max;
	col_vec = col_min:col_max;
	x = (0:n_anos-1)';
	k = 1;
	stack = cell(1, 3);		% Obviously not enough but will shut up MLint
	for (n = col_vec)
		IN = inpolygon(repmat(n, numel(row_vec), 1), row_vec, B{1}(:,2), B{1}(:,1));		% See if ...
		this_row = 1;
		for (m = row_vec)
			if (~IN(this_row)),		continue,	end				% This pixel is outside polygon POI
			this_row = this_row + 1;
			y = double(squeeze(Tmed(m,n,:)));
			ind = isnan(y);
			y(ind) = [];
			if (numel(y) < n_anos/2),		continue,	end		% Completely ad-hoc test
			p = trend1d_m([x(~ind) y],'-L','-N2r','-R','-P');
			stack{k,1} = [x(~ind)+1 y];	% x,temp
			stack{k,2} = p(1);			% Slope
			stack{k,3} = p(4);			% p-value
			k = k + 1;
		end
	end

	str1 = {'*.dat;*.DAT', 'Symbol file (*.dat,*.DAT)'; '*.*', 'All Files (*.*)'};
	[FileName,PathName] = put_or_get_file(handles,str1,'Select Output File name','put','.dat');
	if isequal(FileName,0),		return,		end
	f_name = [PathName FileName];
	double2ascii(f_name, stack, '%.0f\t%.3f', 'maybeMultis');
	[PATH, FNAME, EXT] = fileparts(f_name);
	f_name = [PATH filesep FNAME '_mp' EXT];		% Write a second file with 2 columns where 1st col is slope and 2nth is p-value
	xy = [cat(1,stack{:,2}) cat(1, stack{:,3})];
	double2ascii(f_name, xy);

% ------------------------------------------------------------------------------
function applyFlags(handles, fname, flag, nCells, grd_out)
% Check a 3D file against its 3D flags companion and replace values < FLAG to NaN.
% Optionaly interpolate to fill gaps smaller than NCELLS.
%
% FNAME 	name of a netCDF file with quality flags. Obviously this file must be of
% 			the same size as the series under analysis.
%
% FLAG		Threshold quality value. Only values of quality >= FLAG will be taken into account
%			NOTE: For MODIS use negative FLAG. Than, values are retained if quality <= abs(FLAG)
%
% OPTIONS:
% NCELLS	Holes (bad data) groups smaller than this are filled by interpolation (default = 0). 
%			For example if NCELL = 200 groups of equal or less than a total of 200 will be filled.
%
% GRD_OUT	Name of the netCDF file where to store the result. If not provided, it will be asked here.

	if (nargin < 3),	error('calc_yearMean:Not enough input arguments'),	end

	[z_id, s, rows, cols] = get_ncInfos(handles);

	[s_flags, z_id_flags, msg] = checkFlags_compat(fname, handles.number_of_timesteps, rows, cols);
	if (~isempty(msg)),	errordlg(msg, 'Error'),		return,		end
	if (nargin == 3),	nCells = 0;		end			% If not provided, defaults to no gaps fill

	if (nargin < 5)
		txt1 = 'netCDF grid format (*.nc,*.grd)';	txt2 = 'Select output netCDF grid';
		[FileName,PathName] = put_or_get_file(handles,{'*.nc;*.grd',txt1; '*.*', 'All Files (*.*)'},txt2,'put','.nc');
		if isequal(FileName,0),		return,		end
		grd_out = [PathName FileName];
	end

	% -------------------------------------- END PARSING SECTION --------------------------------------------

	if (flag > 0),		growing_flag = true;				% Pathfinder style (higher the best) quality flag
	else,				growing_flag = false;	flag = -flag;	% MODIS style (lower the best) quality flag
	end
	pintAnoes = (nCells > 0);

	handles.geog = 1;		handles.was_int16 = 0;		handles.computed_grid = 0;
	n_layers = handles.number_of_timesteps;

	aguentabar(0,'title','Applying flags.','CreateCancelBtn');

	for (m = 1:n_layers)
		if (m == 1)						% First layer
			Z = nc_funs('varget', handles.fname, s.Dataset(z_id).Name, [m-1 0 0], [1 rows cols]);
		else
			Z = nc_funs('varget', handles.fname, s.Dataset(z_id).Name, [m-1 0 0], [1 rows cols]);
		end

		Z_flags = nc_funs('varget', fname, s_flags.Dataset(z_id_flags).Name, [m-1 0 0], [1 rows cols]);
		if (growing_flag),		Z(Z_flags < flag) = NaN;	% Pathfinder style (higher the best) quality flag
		else,					Z(Z_flags > flag) = NaN;	% MODIS style (lower the best) quality flag
		end

		ind = isnan(Z);
		if (pintAnoes && any(ind(:)))		% If fill spatial holes is requested
			Z = inpaint_nans(handles, Z, ind, nCells);			% Select interp method inside inpaint_nans()
		end

		% Write this layer to file
		if (m == 1),		nc_io(grd_out, sprintf('w%d/time',n_layers), handles, reshape(Z,[1 size(Z)]))
		else,				nc_io(grd_out, sprintf('w%d', m-1), handles, Z)
		end

		h = aguentabar(m/n_layers,'title','Applying flags.');	drawnow
		if (isnan(h)),	break,	end
	end

% ------------------------1-------2-------3------4------5-------6--------7-------8---------9-----------10--------11----
function calc_yearMean(handles, months, fname2, flag, nCells, fname3, splina, tipoStat, chkPts_file, grd_out, mask_file)
% Compute anual means or climatologies from monthly data (well, and daily too).
%
% MONTHS 	a vector with the months uppon which the mean is to be computed
%		example: 	months = 1:12		==> Computes yearly mean
%					months = 6:8		==> Computes June-July-August seazonal means
%			Default (if months = []) 1:12
%
% OPTIONS:
% FNAME2 	name of a netCDF file with quality flags. Obviously this file must be of
% 			the same size as the series under analysis.
%
% FLAG		Threshold quality value. Only values of quality >= FLAG will be taken into account
%			NOTE: For MODIS use negative FLAG. Than, values are retained if quality <= abs(FLAG)
%
% NCELLS	Holes (bad data) groups smaller than this are filled by interpolation. 
%			For example if NCELL = 200 groups of equal or less than a total of 200 will be filled.
%
% FNAME3 	Optional name of a netCDF file where interpolated nodes will be set to FLAG
%			and the others retain their FNAME2 value. This corresponds to the promotion
%			of interpolated nodes to quality FLAG.
%
% SPLINA	Logical that if true instructs to spline interpolate the missing monthly values
%			before computing the yearly mean. Optionaly, it may be a 2 elements vector with
%			the MIN and MAX values allowed on the Z function (default [0 32]).
%			----------------------- OR (to CLIMA) --------------------
%			A string containing 'CLIMA' to instruct to compute climatologies from montly data
%			This is a uggly hack but I don't want to add yet another input argument.
%
% TIPOSTAT	Variable to control what statistic to compute.
%			0 Compute MEAN of MONTHS period. 1 Compute MINimum and 2 compute MAXimum
%
% CHKPTS_FILE	(Optional)
%			Name of a file with Lon,Lat locations where to output the entire time series.
%			Output name file is constructed from the input name and appended '_tseries'.
%			Fitst column has the month number, even columns the original data and odd columns
%			the data with the holes (NaNs) interpolated with an Akima spline function.
%
% GRD_OUT	Name of the netCDF file where to store the result. If not provided, it will be asked here.
%
% MASK_FILE	name of Land mask file. If it has NaNs it will be multiplied but if only 1/0's, 0's will become NaNs

	% Variables that are not always used but need to exist
	Tmed = [];	ZtoSpline = [];	contanoes = [];		total_months = [];	n_pad_months = [];
	z_id_flags = [];	s_flags= [];

	do_flags = false;		track_filled = false;		do_saveSeries = false;	do_climatologies = false;
	[z_id, s, rows, cols] = get_ncInfos(handles);

	% Find if we are dealing with a Pathfinder V5.2 daily file
	is_PFV52 = PFV52(handles, s);

	if (isempty(months)),	months = 1:12;		end

	if (nargin >= 6 && ~isempty(fname3)),		track_filled = true;	end		% Keep track of interpolated nodes

	if (nargin >= 3 && ~isempty(fname2))				% We have a quality-flag ghost file to check
		[s_flags, z_id_flags, msg] = checkFlags_compat(fname2, handles.number_of_timesteps, rows, cols);
		if (~isempty(msg)),	errordlg(msg, 'Error'),		return,		end
		if (nargin == 3),	nCells = 0;		flag = 7;	end			% If not provided, defaults to best quality
		do_flags = true;
	end

	if (nargin < 7)
		splina = false;		tipoStat = 0;	chkPts_file = [];	grd_out = [];
	elseif (nargin < 8)
		tipoStat = 0;	chkPts_file = [];	grd_out = [];
	end
	if (nargin < 11),	mask_file = '';		end			% Need to test this before the "chkPts_file"
	if (~isempty(mask_file)),	chkPts_file = [];	end

	% -------------- Test if output time series at locations provided in the CHKPTS_FILE --------------------
	if (nargin >= 9 && ~isempty(chkPts_file))
		if (exist(chkPts_file,'file') == 2)
			pts_pos = text_read(chkPts_file);
			indTimeSeries = zeros(size(pts_pos,1), 2);
			for (k = 1:size(pts_pos,1))
				indTimeSeries(k,1) = round((pts_pos(k,1) - handles.head(1)) / handles.head(8)) + 1;
				indTimeSeries(k,2) = round((pts_pos(k,2) - handles.head(3)) / handles.head(9)) + 1;
			end
			%timeSeries = zeros(s.Dataset(3).Size, k);		% Total number of layers
			timeSeries = [(1:s.Dataset(3).Size)' zeros(s.Dataset(3).Size, 2*k)];	% Total N of layers + N of check pts
			indTSCurr_o = 1;		indTSCurr_s = 1;

			if (rem(size(timeSeries,1), 12) ~= 0)
				warndlg('Output time series works only with complete years of monthly data. Ignoring request','Warning')
			else
				do_saveSeries = true;
			end
		end
	end
	if (nargin < 10 || isempty(grd_out))	% Note: old and simple CASE 3 in main cannot send here the output name 
		txt1 = 'netCDF grid format (*.nc,*.grd)';	txt2 = 'Select output netCDF grid';
		[FileName,PathName] = put_or_get_file(handles,{'*.nc;*.grd',txt1; '*.*', 'All Files (*.*)'},txt2,'put','.nc');
		if isequal(FileName,0),		return,		end
		grd_out = [PathName FileName];
	end
	mask = [];		% If needed more than once, this var will hold the amsking array

	% -------------------------------------------------------------------------------------------------------
	% -------------------------------------- END PARSING SECTION --------------------------------------------
	% -------------------------------------------------------------------------------------------------------

	if (flag > 0),		growing_flag = true;				% Pathfinder style (higher the best) quality flag
	else,				growing_flag = false;	flag = -flag;	% MODIS style (lower the best) quality flag
	end
	pintAnoes = (nCells > 0);

	% The following limits are used to clip unreasonable temperatures computed during the spline interpolation
	% When no time (spline) interpolation is used, they are simply ignored
	if (isa(splina, 'char'))
		if (strcmpi(splina,'CLIMA'))
			do_climatologies = true;
		end
		splina = false;					% Make it logic again for use in IF tests
	elseif (numel(splina) == 2)
		regionalMIN = splina(1);		regionalMAX = splina(2);
		splina = true;					% Make it logic again for use in IF tests
	else
		regionalMIN = 0;				regionalMAX = 32;
	end

	handles.geog = 1;		handles.was_int16 = 0;		handles.computed_grid = 0;

	if (rem(handles.number_of_timesteps, 12) == 0 && handles.number_of_timesteps < 336)	% SHIT
		n_anos = handles.number_of_timesteps / 12;
	else
		n_anos = 1;		% TEMPORARY	FOR COMPUTING YEARLY MEANS FROM DAILY DATA --- NON SECURED AND NON DOCUMENTED
	end
	if (is_PFV52)		% This relies on the fact that the apropriate 'description' global attribute has been set
		anos = fix(handles.time);
		n_anos = max(anos) - min(anos) + 1;
		if (n_anos == 2 && anos(1) ~= anos(2))		% Crazy NOAAs have year of 2006 start a 2005.99xxx
			n_anos = 1;
		end

		if (months(end) > n_anos && months(end) < 20)	% A crude test to detect a year-request-overflow (it may screw)
			months = months(1):n_anos;
			if (isempty(months))
				errordlg('Nonsense: Requested period is not covered by the input data','Error'),		return
			else
				warndlg('You requested more years than the data actually has. Reseting max to max years in data.','Warning')
			end
		end
		if (months(end) <= n_anos)					% Convert a period given in number of years to true years.
			months = (anos(1) - months(1) + 1):(anos(end) - months(end) + 1);
		end
		if (months(1) >= anos(1) && months(end) <= anos(end))	% In this case 'months' are actually the years
			new_months = [];
			for (k = 1:numel(months))			% So compute the new 'months' vector
				this_year = months(k);
				ind = find(anos == this_year);
				new_months = [new_months(1:end) ind(:)'];
			end
			months = new_months;
		end
	end

	% Take care of the case when output file has no extension
	[p,f,e] = fileparts(grd_out);
	if (isempty(e))
		if (n_anos == 1),	e = '.grd';
		else,				e = '.nc';
		end
		grd_out = [p filesep f e];
	end

	if (splina)
		n_pad_months = 7;		% Example: if months = 7:9 interpolation domain is 7-n_pad_months-1:9+n_pad_months
		if (n_anos == 1)
			n_pad_months = 0;
			warndlg('Time series is too short to do the spline interpolation option.','WARNING')
		end
		total_months = numel(months) + 2*n_pad_months;
		ZtoSpline = alloc_mex(rows, cols, total_months, 'single', NaN);
	end

	aguentabar(0,'title','Computing means.','CreateCancelBtn');
	if (~splina)
		Tmed = zeros(rows, cols);		% Temp media para cada um dos anos
	end
	in_break = false;					% Inner loop cancel option
	last_processed_month = 0;		already_processed = 0;

	if (do_climatologies)
		multi_climat = false;		n_outer_loop = 1;
		if (months(1) > 100)			% Compute a climat for each element in 'months'
			multi_climat = true;	n_outer_loop = numel(months);
		end
	else
		n_outer_loop = n_anos;
	end

	for (m = 1:n_outer_loop)

		if (splina)
			if (m == 1)
				past_months = max(1,months(1)-n_pad_months) - 1;
				this_months = past_months+1 : (months(end)+n_pad_months);
			elseif (m == n_anos)
				past_months = (m - 1)*12 + months(1) - n_pad_months - 1;
				this_months = past_months+1 : past_months+(n_pad_months + numel(months) + min(n_pad_months, 12-months(end)));
			else
				past_months = (m - 1)*12 + months(1) - n_pad_months - 1;
				this_months = past_months+1 : past_months+(numel(months)+2*n_pad_months);
			end
		elseif (do_climatologies)
			% Litle manip to create a row vector with all months of interest across the years
			if (multi_climat)
				this_months = repmat(months(m)-100, n_anos, 1);
			else
				this_months = repmat(months, n_anos, 1);
			end
			v = ((0:n_anos) * 12)';
			for (k = 1:n_anos)
				this_months(k,:) = this_months(k,:) + v(k);
			end
			this_months = this_months';
			this_months = this_months(:)';
			contanoes = zeros(rows, cols);
		else
			this_months = (m - 1) * 12 + months;
			contanoes = zeros(rows, cols);
		end

		if (tipoStat > 0)
			l_inc = this_months(2) - this_months(1);
			s2 = struct('fname',handles.fname, 'info',s, 'n_layers', handles.number_of_timesteps, 'rows',rows, ...
				 'cols',cols, 'layerOI',this_months(1), 'layer_inc',l_inc, 'z_id',z_id);
			tmp = doM_or_M_or_M([], this_months(1), l_inc, handles.number_of_timesteps, Inf, Inf, tipoStat, s2);
		else
			% For averages and for the time being (not break compat), continue to use the old code in form of function
			[Tmed, contanoes, ZtoSpline, already_processed] = ...
				calc_average_old(handles, s, s_flags, Tmed, ZtoSpline, contanoes, this_months, splina, last_processed_month, ...
				nCells, n_anos, already_processed, total_months, months, n_pad_months, do_flags, flag, ...
				fname2, fname3, pintAnoes, track_filled, rows, cols, z_id, z_id_flags, growing_flag, m);
			if (~splina)					% Do not interpolate along time. Compute averages with all non NaNs
				cvlib_mex('div', Tmed, contanoes);		% The mean for current year
				tmp = single(Tmed);
			end
		end

		last_processed_month = this_months(end);

		if (tipoStat == 0 && splina)					% Fill missing month data by interpolation based on non-NaN data
			hh = aguentabar(eps,'title','Splining it.');	drawnow
			if (isnan(hh)),		break,		end			% Over time loop said: break
			n_meses = numel(this_months);

			if (m == 1),	first_wanted_month = months(1);				% First year in the stack
			else,			first_wanted_month = n_pad_months + 1;
			end
			last_wanted_month = first_wanted_month + numel(months) - 1;

			yy = zeros(1, n_meses)';				% A kind of pre-allocation to be used by the akimaspline MEX

			% ----------- Test if we are saving time series in array for later saving ----------------------
			if (do_saveSeries)
				[timeSeries, indTSCurr_o] = getTimeSeries(ZtoSpline, timeSeries, indTimeSeries, indTSCurr_o, ...
											true, first_wanted_month, last_wanted_month);
			end
			% ----------------------------------------------------------------------------------------------

			for (i = 1:cols)
				for (j = 1:rows)
					if (already_processed),		break,	end
					y = double(squeeze(ZtoSpline(j,i,1:n_meses)));
					ind = ~isnan(y);
					% If have NaNs inside months of interest and the overall series has enough points, interp in the missing positions
					if (all(ind(first_wanted_month:last_wanted_month)))		% We have them all, so nothing to interp
						ZtoSpline(j,i,1:n_meses) = single(y);
					elseif (~any(ind))		% They are all NaNs -- Almost sure a land pixel
						continue
					elseif (numel(ind(~ind)) <= round((numel(this_months) - (n_pad_months * (m ~= 1))) / 2))
						% At least > 1/2 number of valid pts not counting first n_pad_months that were already interpolated
						x = this_months(ind);			y0 = y(ind);
 						akimaspline(x, y0, this_months, yy);
						y(first_wanted_month:last_wanted_month) = yy(first_wanted_month:last_wanted_month);
						if (~ind(1) || ~ind(end))			% Cases when we would have extrapolations
							if (~ind(1))
								ki = 1;
								while (~ind(ki)),	ki = ki + 1;	end
								y(1:ki) = NaN;				% Reset extraped values to NaN
							end
							if (~ind(end))
								kf = numel(ind);
								while (~ind(kf)),	kf = kf - 1;	end
								y(kf:numel(ind)) = NaN;		% Reset extraped values to NaN
							end
						end
						ZtoSpline(j,i,1:n_meses) = single(y);
					else									% Less than half valid points. We'll make them be all NaNs
						if (~tipoStat),	 ZtoSpline(j,i,1:n_meses) = y * NaN;	end		% For MIN & MAX we want to retain data
					end
				end
				hh = aguentabar(i/(cols+1));	drawnow
				if (isnan(hh)),		in_break = true;	break,		end		% Over time loop said: break (comment: FCK ML SHIT)
			end				% End loops over this 2D layer

			if (in_break),		break,		end			% Fck no gotos paranoia obliges to this recursive break

			% ----------- Test if we are saving time series in array for later saving ----------------------
			if (do_saveSeries)
				[timeSeries, indTSCurr_s] = getTimeSeries(ZtoSpline, timeSeries, indTimeSeries, indTSCurr_s, ...
											false, first_wanted_month, last_wanted_month);
			end
			% ----------------------------------------------------------------------------------------------

			% Now we can finaly compute the season MEAN or MIN or MAX
			tmp = doM_or_M_or_M(ZtoSpline, first_wanted_month, 1, last_wanted_month, regionalMIN, regionalMAX, tipoStat);
		end							% End interpolate along time (SPLINA)

		tmp(tmp == 0) = NaN;		% Reset the NaNs

		if (~isempty(mask_file))
			[tmp, mask] = apply_mask(mask_file, mask, tmp);		% First time reads from MASK_FILE, second on uses MASK
		end

		if (in_break),		break,		end		% Fckng no gotos paranoia obliges to this recursive break

% 		% Clip obvious bad data based on cheap median statistics
% 		aguentabar(0.5,'title','Filtering obvious bad data based on cheap statistics.');	drawnow
% 		medianas = c_grdfilter(tmp,[1 size(tmp,2) 1 size(tmp,1) 0 50 0 1 1], '-D0', '-Fm11');
% 		difa = cvlib_mex('absDiff', tmp, medianas);
% 		tmp(difa > 0.75) = NaN;					% 0.75 is probably still too permissive

		% Write this layer to file
		if (m == 1)
			if (n_outer_loop == 1)		% If one single layer don't make it 3D
				% Defaults and srsWKT fishing are set in nc_io
				misc = struct('x_units',[],'y_units',[],'z_units',[],'z_name',[],'desc',[], ...
					'title','Yearly or Seazonal average','history',[],'srsWKT',[], 'strPROJ4',[]);
				zz = grdutils(tmp,'-L');
				handles.head(5:6) = [double(zz(1)) double(zz(2))];
				nc_io(grd_out, 'w', handles, tmp, misc)
			else
				nc_io(grd_out, sprintf('w%d/time',n_outer_loop), handles, reshape(tmp,[1 size(tmp)]))
			end
		else
			nc_io(grd_out, sprintf('w%d', m-1), handles, tmp)
		end

		h = aguentabar(m/n_outer_loop,'title','Computing means.');	drawnow
		if (isnan(h)),	break,	end

		if (~splina && m < n_outer_loop),	cvlib_mex('CvtScale', Tmed, 0.0, 0.0);	end			% Reset it to zeros
	end

	if (do_saveSeries)			% Save the time series file. The name is build from that of locations file
		[pato, fname, ext] = fileparts(chkPts_file);
		fname = [fname '_tseries' ext];
		if (~isempty(pato)),	fname = [pato filesep fname];	end
		double2ascii(fname, timeSeries, ['%d' repmat('\t%.4f',[1 size(timeSeries,2)-1])]);
	end

% ------------------------------------------------------------------------------
function [Tmed, contanoes, ZtoSpline, already_processed] = ...
		calc_average_old(handles, s, s_flags, Tmed, ZtoSpline, contanoes, this_months, splina, last_processed_month, ...
		nCells, n_anos, already_processed, total_months, months, n_pad_months, do_flags, flag, ...
		fname2, fname3, pintAnoes, track_filled, rows, cols, z_id, z_id_flags, growing_flag, m)
% Chunk of code that calculates the average in the old way and converted to a function.

	counter = 0;
	for (n = this_months)
		counter = counter + 1;
		if (m == 1)						% First year
			Z = nc_funs('varget', handles.fname, s.Dataset(z_id).Name, [n-1 0 0], [1 rows cols]);
		else
			if (splina && n <= last_processed_month)
				already_processed = already_processed + 1;
				offset = 1;
				if (m == 2),	offset = min(1, months(1)-n_pad_months);	end		% Because for the 1st year the series may be shorter
				offset = total_months - (last_processed_month - n) + offset - 1;
				ZtoSpline(:,:,already_processed) = ZtoSpline(:,:,offset);
			else
				Z = nc_funs('varget', handles.fname, s.Dataset(z_id).Name, [n-1 0 0], [1 rows cols]);
				already_processed = 0;
			end
		end

		if (do_flags && ~already_processed)
			Z_flags = nc_funs('varget', fname2, s_flags.Dataset(z_id_flags).Name, [n-1 0 0], [1 rows cols]);
			if (growing_flag),		Z(Z_flags < flag) = NaN;	% Pathfinder style (higher the best) quality flag
			else,					Z(Z_flags > flag) = NaN;	% MODIS style (lower the best) quality flag
			end
		end

		ind = isnan(Z);

		if (pintAnoes && ~already_processed && any(ind(:)))		% If fill spatial holes is requested
			if (track_filled),		ind0 = ind;		end			% Get this Z level original NaNs mask
			Z = inpaint_nans(handles, Z, ind, nCells);			% Select interp method inside inpaint_nans()
			ind = isnan(Z);
			if (track_filled && counter <= 12)					% Write updated quality file (The 'splina' case has counters >> 12)
				mn = (m - 1)*12 + counter - 1;					% This will work only for entire years (not seasons)
				Z_flags(ind0 & ~ind) = flag;					% Promote interpolated pixels to quality 'flag'
				grdutils(Z_flags,'-c');							% Shift by -128 so it goes well with the uint8 add_off elsewere
				if (mn == 0),		nc_io(fname3, sprintf('w%d/time',n_anos*numel(months)), handles, reshape(Z_flags,[1 size(Z_flags)]))
				else,				nc_io(fname3, sprintf('w%d', mn), handles, Z_flags)
				end
			end
		end

		if (~splina)				% Do not interpolate along time (months)
			Z(ind) = 0;				% Transmutate the Anoes
			contanoes = contanoes + ~ind;
			cvlib_mex('add', Tmed, double(Z));
		else						% Pack this year into a 3D temporary variable, to be processed later.
			if (~already_processed),	ZtoSpline(:,:,counter) = Z;		end
		end

		if (n_anos == 1 && numel(this_months) > 12)			% For the secret daily data case
			aguentabar(n / (numel(this_months) + 1)),		drawnow
		end
	end								% End loop over months

% ------------------------1----------2-------3----------4---------5---------6-------7---
function calc_L2_periods(handles, period, tipoStat, regMinMax, grd_out, mask_file, filt)
% Compute averages for 1, 3, 8, month periods of L2 data processed by empilhador
%
% PERIOD	Number of days of the composit period (e.g. 3, 8, 30). A variable number
%			of layers per day is allowed.
%			Alternatively PERIOD can be a vector with the periods to compute, and they don't
%			need to have all the same interval. But ATTENTION that the elements of this
%			vector need to cover the dates in the 3rth dim of the input nc file.
%			This allows selecting a fix start. E.G this will compute 3 days composits of
%			the first 29 days of January 2012.
%				PERIOD=2012.0:3/365:(2012+29/365);
%			Off course dates in input must fall in this period, otherwise output is NaNs only
%
% TIPOSTAT	Variable to control what statistic to compute.
%			0 Compute MEAN of MONTHS period. 1 Compute MINimum, 2 compute MAXimum, 3 compute STD
%
% REGMINMAX	A 2 elements vector with the MIN and MAX values allowed on the Z function (default [0 inf])
%
% GRD_OUT	Name of the netCDF file where to store the result. If not provided, it will be asked here.
%
% MASK_FILE	name of Land mask file. If it has NaNs it will be multiplied but if only 1/0's, 0's will become NaNs
%
% FILT      If provided apply a gaussian filtering with FILT filtering full width

	if (nargin < 5 || isempty(grd_out))			% Note: old and simple CASE 3 in main cannot send here the output name 
		txt1 = 'netCDF grid format (*.nc,*.grd)';	txt2 = 'Select output netCDF grid';
		[FileName,PathName] = put_or_get_file(handles,{'*.nc;*.grd',txt1; '*.*', 'All Files (*.*)'},txt2,'put','.nc');
		if isequal(FileName,0),		return,		end
		grd_out = [PathName FileName];
	end
	if (nargin < 4)
		regionalMIN = -2;	regionalMAX = inf;
	else
		regionalMIN = regMinMax(1);		regionalMAX = regMinMax(2);
		if (regionalMIN == 0 && regionalMAX == 0),	regionalMIN = -Inf;	regionalMAX = Inf;	end
	end
	if (nargin < 6),	mask_file = '';		end
	mask = [];			% If needed more than once, this var will hold the masking array
	if (nargin < 7),	filt = 0;			end

	[z_id, s, rows, cols] = get_ncInfos(handles);

	% Find if we are dealing with a Pathfinder V5.2 daily file
	[is_PFV52, tempos] = PFV52(handles, s);	% if not a PFV5.2, tempos = handles.time(:)

	if (numel(period) == 1)
		periods = fix(tempos(1)) : period : fix(tempos(end))+period-1;	% Don't risk to loose an incomplete last interval
		half_period = repmat(period/2, 1, numel(periods));		% for naming layers in nc file
	else
		periods = period;
		% Need to add 1 because HISTC counts in the interval edge(k) <= x < edge(k+1) and, e.g., day 3.7 is still day 3
		periods(2:end) = periods(2:end) + 1;
		half_period = [diff(periods(:)')/2 (periods(end) - periods(end-1))/2];	% Repeat last value
	end

	N = histc(fix(tempos), periods);
	N(end) = [];		periods(end) = [];		half_period(end) = [];	% Last N is for >= edge(end) and we don't want it
	if (isempty(periods))
		warndlg('There is nothing inside the period(s) you have requested. Bye.', 'Warning'),	return
	end

	handles.was_int16 = false;		% I have to get rid of the need to set this

	aguentabar(0,'title','Computing period means.','CreateCancelBtn');

	% C is the counter to the current layer number being processed.
	c = find(fix(tempos) < periods(1));		% Find the starting layer number
	if (isempty(c)),	c = 1;				% We start at the begining of file.
	else,				c = c(end) + 1;		% We start somewhere at the middle of file.
	end

	n_periods = numel(periods);
	for (m = 1:n_periods)
		if (N(m) ~= 0)
% 			Z = alloc_mex(rows, cols, N(m), 'single', NaN);
% 			for (n = 1:N(m))			% Loop over the days in current period
% 				Z(:,:,n) = nc_funs('varget', handles.fname, s.Dataset(z_id).Name, [c-1 0 0], [1 rows cols]);
% 				c = c + 1;
% 			end
% 			tmp = doM_or_M_or_M(Z, 1, size(Z,3), regionalMIN, regionalMAX, tipoStat);
% 			clear Z;					% Free memory (need because at the alloc time above, two of them would exist)
			s2 = struct('fname',handles.fname, 'info',s, 'n_layers', N(m), 'rows',rows, 'cols',cols, ...
				'layerOI',c, 'layer_inc',1, 'z_id',z_id);
			[tmp, s2] = doM_or_M_or_M([], 1, 1, N(m), regionalMIN, regionalMAX, tipoStat, s2);
			c = s2.layerOI;				% This is crutial because it tells us where we are in the layer stack
			tmp(tmp == 0) = NaN;		% Reset the NaNs
			if (~isempty(mask_file))
				[tmp, mask] = apply_mask(mask_file, mask, tmp);		% First time reads from MASK_FILE, second on uses MASK
			end

			if (filt)
				tmp = c_grdfilter(tmp, handles.head, sprintf('-Fg%f -Nr', filt), '-D3');
			end

			zzz = grdutils(tmp,'-L');
			handles.head(5) = min(handles.head(5), zzz(1));		handles.head(6) = max(handles.head(6), zzz(2));
		else
			tmp = alloc_mex(rows, cols, 1, 'single', NaN);
		end

		% Compute the mean time of this bin and use it to name the layer
		thisLevel = periods(m) + half_period(m);

		% Write this layer to file, but must treate compiled version differently since
		% it is not able to write UNLIMITED files
		if (true || ~handles.IamCompiled)		% TEMP. Later, if it works, we'll simply delete the other branch
			if (m == 1),	nc_io(grd_out, sprintf('w-%f/time',thisLevel), handles, reshape(tmp,[1 size(tmp)]))
			else,			nc_io(grd_out, sprintf('w%d\\%f', m-1, thisLevel), handles, tmp)
			end
		else
			if (m == 1)
				handles.levelVec = periods + half_period;
     			nc_io(grd_out,sprintf('w%d/time',n_periods), handles, reshape(tmp,[1 size(tmp)]))
			else
				nc_io(grd_out, sprintf('w%d', m-1), handles, tmp)
			end
		end

		h = aguentabar(m/n_periods,'title','Computing period means.');	drawnow
		if (isnan(h)),	break,	end
	end
	
% --------------------------------------------------------------------------------------
function [Z, mask] = apply_mask(mask_file, mask, Z)
% When MASK is empty, read data from MASK_FILE otherwise just use MASK to mask out Z.
% Masking is done either by multiplication by MASK or by logical op when MASK is a logical (or only 1 & 0's)
	if (isempty(mask))
		try
			G = gmtmex(['read -Tg ' mask_file]);
			if ~((size(G.z,1) == size(Z,1)) && (size(G.z,2) == size(Z,2)))
				errordlg('The sizes of the Land mask grid and the input array differ. Ignoring masking request.', 'Error')
				return
			end
		catch
			errordlg(sprintf('Error reading file %s\n%s', mask_file, lasterror),'Error')
			return
		end
		if (G.range(5) == 0 && G.range(6) == 1)		% A 1/0's mask
			if (nargout == 2)			% If it's going to be reused, better convert it to logicals right away
				mask = logical(G.z);
				clear G
			end
			Z(mask) = NaN;
		else
			cvlib_mex('mul', Z, G.z)
		end
	else
		if (isa(mask, 'logical'))
			Z(mask) = NaN;
		else
			cvlib_mex('mul', Z, mask)
		end
	end

% ------------------------------------------------------------------------------------
function Tavg = tideman(handles, Temps, tempos, nPeriods)
% Apply the method suggested by TideMan in this post
% http://www.mathworks.com/matlabcentral/newsreader/view_thread/292502#887770
% but use weeks instead of months which would introduce a significant phase
% offset because a month is too long period for this method

	if (nargin == 3),	nPeriods = 10;	end
	rPeriod = 1 / nPeriods;
	
	decTime = tempos - fix(tempos);				% The decimal part only
	Tavg = zeros(size(Temps,1), size(Temps,2), nPeriods);
	for (k = 1:nPeriods)
		s = (k - 1) * rPeriod;		e = k * rPeriod;
		ind = (decTime >= s & decTime < e);
		thisStackPeriod = Temps(:,:,ind);		% All layers that fall inside this (stacking) period
		Tavg(:,:,k) = doM_or_M_or_M(thisStackPeriod);
	end

% ------------------------------------------------------------------------------------
function T_measured = remove_seazon(handles, T_measured, Tavg, tempos)
% T_MEASURED is the array with measured temperatures (or something else)
% TAVG is the mean cycle computed by the tideman function
% Subtract both (after reinterpolations) and return result in inplace modified T_MEASURED

	nYears = numel(find(diff(fix(tempos)) ~= 0)) + 1;
	nPeriods = size(Tavg,3);			% Number of intervals in which the seazon model has been computed
	rPeriod = 1 / nPeriods;
	ty = linspace(0+rPeriod/2, 1-rPeriod/2, nPeriods)';	% Decimal time centered in the middle of each interval
	t = zeros(nPeriods * nYears, 1);
	for (k = 1:nYears)
		t((k-1)*nPeriods+1 : k*nPeriods) = fix(tempos(1)) + (k-1) + ty;
	end

	y = zeros(1, numel(tempos));
	for (m = 1:size(T_measured,1))
		for (n = 1:size(T_measured,2))
			T_seazon = double(repmat(squeeze(Tavg(m,n,:)), nYears, 1));		% Replicate seazonal cycle
			tmp = double(squeeze(T_measured(m,n,:)));
			indNan = isnan(tmp);
			%y = interp1(t, T_seazon, tempos, 'linear', 'extrap');
			akimaspline(t, T_seazon, tempos, y);		% we don't need a spline but it's much faster
			y(indNan) = NaN;
			T_measured(m,n,:) = single(tmp - y(:));
		end
	end

% ----------------------------------------------------------------------
function [out, s] = get_layer(Z, layer, s)
% Get the layer in one of the following two instances
%	1 - S is empty and Z is [m n p]
%	2 - S is a structure with info on how to read the layer directly from the netCDF file

	if (isempty(s))
		out = Z(:,:,layer);
	else
		out = nc_funs('varget', s.fname, s.info.Dataset(s.z_id).Name, [s.layerOI-1 0 0], [1 s.rows s.cols]);
		s.layerOI = s.layerOI + s.layer_inc;
	end

% ----------------------------------------------------------------------
function [out, s] = doM_or_M_or_M(Z, first_level, lev_inc, last_level, regionalMIN, regionalMAX, tipo, s)
% Compute either the MEAN (TIPO = 0) or the MIN (TIPO = 2), MAX (3) or STD (4) of the period selected
% by the first_level:lev_inc:last_level vector. Normaly a year but can be a season as well.
% NOTE1: This function was only used when SPLINA (see above in calc_yearMean()) up to Mirone 2.2.0
% NOTE2: It is now (2.5.0dev) used again by the tideman function (and other calls)
%
% Because of the potentially very large memory comsumption, we can do the data file reading from
% within this function. In that case Z can be empty and S must be a struct with:
%	S = struct('fname',handles.fname, 'info',s, 'n_layers', N(m), 'rows',rows, 'cols',cols, ...
%	           'layerOI',c, layer_inc,n, 'z_id',z_id);
% whre N_LAYERS is the numbers of layers to be read. 'layerOI' flags the layer to be read and must be
% incremented after each layer reading (which is done by the GET_LAYER() function).
% LAYER_INC is the increment between layers. It will be = 1 for consecutive layers, or 12 for climatologies
% We also return S because it keeps the trace of the current layer (the layerOI field), which is needed
% when there are multiple calls to this function.

	if (nargin == 1)			% Compute average of all layers without any constraint
		first_level = 1;		lev_inc = 1;			last_level = size(Z,3);
		regionalMIN = Inf;		regionalMAX = Inf;		tipo = 0;
	end
	if (nargin < 8),	s = [];	end

	if (tipo == 0)				% Compute the MEAN of the considered period
		% We don't use nanmean_j here because of the regionalMIN|MAX
		[out, s] = get_layer(Z, first_level, s);
		out(out < regionalMIN | out > regionalMAX) = NaN;
		ind = isnan(out);
		contanoes = alloc_mex(size(ind,1), size(ind,2), 'single');
		cvlib_mex('add', contanoes, single(~ind));
		out(ind) = 0;						% Mutate NaNs to 0 so that they don't screw the adition
		for (n = (first_level+lev_inc):lev_inc:last_level)
			[tmp, s] = get_layer(Z, n, s);
			if (~isinf(regionalMIN)),	tmp(tmp < regionalMIN) = NaN;	end
			if (~isinf(regionalMAX)),	tmp(tmp > regionalMAX) = NaN;	end
			ind = isnan(tmp);
			tmp(ind) = 0;
			cvlib_mex('add', contanoes, single(~ind));
			cvlib_mex('add', out, tmp);
		end
		cvlib_mex('div', out, contanoes);			% The mean
	elseif (tipo == 1)						% Median
		% too complicated if the entire dataset is not on memory 
	elseif (tipo == 2 || tipo == 3)			% ...
		v = version;
		if (str2double(v(1)) > 6)			% With R14 and above we use the built in min and max
			if (tipo == 2),		fh = @min;	% Minimum of the selected period
			else,				fh = @max;	% Maximum of the selected period
			end
		else								% But for R13 and compiled we must avoid the BUGGY NaNs comparisons
			if (tipo == 2),		fh = @min_nan;		test = 1e10;
			else,				fh = @max_nan;		test = -1e10;
			end
		end
		if (isempty(s))
			out = feval(fh, Z(:,:,first_level:lev_inc:last_level),[],3);
		else
			[out, s] = get_layer(Z, first_level, s);
			for (k = first_level+lev_inc:lev_inc:s.n_layers)
				[tmp, s] = get_layer(Z, k, s);
				out = feval(fh, out, tmp);
			end
			if (str2double(v(1)) <= 6)
				out(out == test) = NaN;		% Reset the ests values back to NaNs
			end
		end
	elseif (tipo == 4)			% STD
		out = nanstd_j(Z, first_level, lev_inc, last_level, s);
	end

% ----------------------------------------------------------------------
function out = min_nan(A, B)
% Compute minimum of a 3D array or two 2D arrays that have NaNs.
% This function is used only by R13 and compiled versionn that very bugged in whta concerns NaNs comparisons
% Fot the two arguments case, the caller function is responsible to reset the 1e10 to NaNs again.
% This to allow calls in a loop and avoid intermediate wasting replacements.
	if (nargin == 1)
		A(isnan(A)) = 1e10;
		out = min(A, [], 3);
		out(out == 1e10) = nan;
	else
		A(isnan(A)) = 1e10;		B(isnan(B)) = 1e10;
		out = min(A, B);
	end

% ----------------------------------------------------------------------
function out = max_nan(A, B)
% Compute maximum of a 3D array or two 2D arrays that have NaNs.
% This function is used only by R13 and compiled versionn that very bugged in whta concerns NaNs comparisons
% Fot the two arguments case, the caller function is responsible to reset the -1e10 to NaNs again.
% This to allow calls in a loop and avoid intermediate wasting replacements.
	if (nargin == 1)
		A(isnan(A)) = -1e10;
		out = max(A, [], 3);
		out(out == -1e10) = nan;
	else
		A(isnan(A)) = -1e10;	B(isnan(B)) = -1e10;
		out = max(A, B);
	end

% ----------------------------------------------------------------------
function out = nanmean_j(Z, first_level, lev_inc, last_level, s)
% ...
	if (nargin == 1)
		first_level = 1;	last_level = size(Z,3);		lev_inc = 1;	s = [];
	end
	[out, s] = get_layer(Z, first_level, s);
	ind = isnan(out);
	contanoes = alloc_mex(size(ind,1), size(ind,2), 'single');
	cvlib_mex('add', contanoes, single(~ind));
	out(ind) = 0;						% Mutate NaNs to 0 so that they don't screw the adition
	for (n = (first_level+lev_inc):lev_inc:last_level)
		[tmp, s] = get_layer(Z, n, s);
		ind = isnan(tmp);
		tmp(ind) = 0;
		cvlib_mex('add', contanoes, single(~ind));
		cvlib_mex('add', out, tmp);
	end
	cvlib_mex('div', out, contanoes);			% The mean

% ----------------------------------------------------------------------
function out = nanstd_j(Z, first_level, lev_inc, last_level, s)
% Compute the STD taking into account the presence of NaNs
% This is a bit more convoluted for memory efficiency concearns (somethig TMW does not care)
%
% S is a structure with info to read the layers from file instead of relying in Z (that hould be [] than)
% For further info, see help section of the doM_or_M_or_M() function

	if (nargin == 1)
		first_level = 1;	last_level = size(Z,3);		lev_inc = 1;	s = [];
	end

	this_mean = nanmean_j(Z, first_level, lev_inc, last_level, s);
	if (isa(this_mean, 'single'))		% Need this gimnastic because cvlib_mex screws if types are different
		out = alloc_mex(size(this_mean,1), size(this_mean,2), 'single');
	else
		out = alloc_mex(size(this_mean,1), size(this_mean,2), 'double');
	end
	denom = zeros(size(this_mean));			% Swallow the thing but this one has to be done with doubles
	for (n = first_level:lev_inc:last_level)
		[t, s] = get_layer(Z, n, s);
		denom = denom + ~isnan(t);
		cvlib_mex('sub', t, this_mean);
		cvlib_mex('mul', t, t);				% The squares of (xi - xm)
		t(isnan(t)) = 0;
		cvlib_mex('add', out, t)
	end

	denom = max(denom-1, 1);		% divide by (n-1). But when n == 0 or 1, we'll return ones
	denom(denom == 0) = NaN;		% When all NaNs return NaN, and thus avoid a divide by 0
	if (isa(this_mean, 'single')),	denom = single(denom);	end

	cvlib_mex('div', out, denom)	
	cvlib_mex('pow', out, 0.5)

% ----------------------------------------------------------------------
function calc_corrcoef(handles, secondArray, sub_set, pValue, grd_out)
% Compute the correlation coefficien between loaded array and 'secondArray'
% 
% NOTE: THIS IS A HIGHLY MEMORY CONSUMPTION ROUTINE AS ALL DATA IS LOADED IN MEMORY
%
% SECONDARRAY	name of the other netCDF file whose correlation with loaded array will be estimated.
% 			Obviously this file must be of the same size as the series under analysis.
%
% OPTIONS:
% SUB_SET -> A two columns row vec with number of the offset of years where analysis start and stop.
%			For example [3 1] Starts analysis on forth year and stops on the before last year.
%			[0 0] Means using the all dataset.
%
% PVALUE	Logical that if true instruct to compute also the p-value of the estimate
%
% GRD_OUT	Name of the netCDF file where to store the result. If not provided, open Mirone Fig.

	n_anos = handles.number_of_timesteps;
	[z_idA, sA, rows, cols] = get_ncInfos(handles);
	if (nargin < 5),	grd_out = [];	end

	if (nargin >= 3 && (numel(sub_set) == 2))
		jump_anos = sub_set(1);		stop_before_end_anos = sub_set(2);
	else
		jump_anos = 0;				stop_before_end_anos = 0;
	end

	n_anos = n_anos - (jump_anos + stop_before_end_anos);	% Number of layers to be used in this run

	[sB, z_idB, msg] = checkFlags_compat(secondArray, handles.number_of_timesteps, rows, cols);
	if (~isempty(msg)),	errordlg(msg, 'Error'),		return,		end
	
	arrayB = alloc_mex(rows, cols, n_anos, 'single');
	arrayA = alloc_mex(rows, cols, n_anos, 'single');
	for (m = 1:n_anos)
 		arrayA(:,:,m) = nc_funs('varget', handles.fname, sA.Dataset(z_idA).Name, [(m - 1 + jump_anos) 0 0], [1 rows cols]);
		arrayB(:,:,m) = nc_funs('varget', secondArray, sB.Dataset(z_idB).Name, [(m - 1 + jump_anos) 0 0], [1 rows cols]);
	end

	h = aguentabar(0,'title','Computing correlation coefficients','CreateCancelBtn');

	try
	Rgrid = alloc_mex(rows, cols, 'single', NaN);	% Works on R13 as well
	if (pValue),	Pgrid = alloc_mex(rows, cols, 'single', NaN);	end
	for (m = 1:rows)
		for (n = 1:cols)
			Var1 = double(squeeze(arrayA(m,n,:)));
			Var2 = double(squeeze(arrayB(m,n,:)));
			
			ind = isnan(Var1);
			Var1(ind) = [];		Var2(ind) = [];
			ind = isnan(Var2);
			Var1(ind) = [];		Var2(ind) = [];
			% Completely ad-hoc test (it also jumps land cells)
			if ((numel(Var1) < n_anos*0.66)),		continue,	end

%  			%[R_, p] = corrcoef([Var1 Var2]);	%
			thisN = numel(Var1);
			VarAB = double([Var1 Var2]);
			VarAB = VarAB - repmat(sum(VarAB,1) / thisN, thisN, 1);	% Remove mean
			R = (VarAB' * VarAB) / (thisN - 1);		% Covariance matrix
			d = sqrt(diag(R));		% sqrt first to avoid under/overflow
			R = R ./ (d*d');		% Correlation matrix

			Rgrid(m,n) = single(R(2));

			if (pValue)
	 			p = trend1d_m(VarAB,'-L','-N2r','-R','-P');
				if (p(1) < -0.5 || p(1) > 1),	continue,	end		% Another ad-hoc (CLIPPING)
				Pgrid(m,n) = p(4);

				% To compute p, the Matlab way (but too recent Matlab)
				% Tstat = +/-Inf and p = 0 if abs(r) == 1, NaN if r == NaN.
% 				Tstat = R(2) .* sqrt((thisN-2) ./ (1 - R(2).^2));
% 				p = zeros(2);
% 				p(2) = 2*tpvalue(-abs(Tstat), thisN-2);
% 				p = p + p' + diag(diag(R)); % Preserve NaNs on diag.
% 				Pgrid(m,n) = single(p(2));
			end

		end
		if (rem(m, 10) == 0),	h = aguentabar(m/rows);		end
		if (isnan(h)),	break,	end
	end
	if (isnan(h)),	return,		end
	if (ishandle(h)),	aguentabar(1),		end		% Make sure it goes away.

catch
	disp(lasterror)
end
	
	clear arrayA arrayB

	zz = grdutils(Rgrid,'-L');  handles.head(5:6) = [zz(1) zz(2)];
	tmp.head = handles.head;
	if (isempty(grd_out))	% Show result in a Mirone figure
		tmp.X = linspace(tmp.head(1),tmp.head(2),cols);
		tmp.Y = linspace(tmp.head(3),tmp.head(4),rows);
		tmp.name = 'Correlation';
		mirone(Rgrid, tmp)
	else					% Got output name from input arg
		handles.was_int16 = false;
		handles.computed_grid = true;
		handles.geog = 1;
		nc_io(grd_out, 'w', handles, Rgrid)
	end

	if (pValue)
		zz = grdutils(Pgrid,'-L');  handles.head(5:6) = [zz(1) zz(2)];
		tmp.head = handles.head;
		tmp.X = linspace(tmp.head(1),tmp.head(2),cols);
		tmp.Y = linspace(tmp.head(3),tmp.head(4),rows);
		tmp.name = 'p-values';
		mirone(Pgrid, tmp)
	end

% ----------------------------------------------------------------------
% function p = tpvalue(x,v)
% %TPVALUE Compute p-value for t statistic.

% normcutoff = 1e7;
% if length(x)~=1 && length(v)==1
%    v = repmat(v,size(x));
% end
% 
% % Initialize P.
% p = NaN(size(x));
% nans = (isnan(x) | ~(0<v)); % v == NaN ==> (0<v) == false
% 
% % First compute F(-|x|).
% %
% % Cauchy distribution.  See Devroye pages 29 and 450.
% cauchy = (v == 1);
% p(cauchy) = .5 + atan(x(cauchy))/pi;
% 
% % Normal Approximation.
% normal = (v > normcutoff);
% p(normal) = 0.5 * erfc(-x(normal) ./ sqrt(2));
% 
% % See Abramowitz and Stegun, formulas 26.5.27 and 26.7.1.
% gen = ~(cauchy | normal | nans);
% p(gen) = betainc(v(gen) ./ (v(gen) + x(gen).^2), v(gen)/2, 0.5)/2;
% 
% % Adjust for x>0.  Right now p<0.5, so this is numerically safe.
% reflect = gen & (x > 0);
% p(reflect) = 1 - p(reflect);
% 
% % Make the result exact for the median.
% p(x == 0 & ~nans) = 0.5;

% ----------------------------------------------------------------------
function [tSeries, indTSCurr] = getTimeSeries(ZtoSpline, tSeries, indTS, indTSCurr, orig, first_month, last_month)
% Fill the TSERIES array with the original as well as the spline interpolated time series.
% This is intended mostly to provide a way to check the goodness (or not) of the interp mechanism
%
% ZtoSpline The 3D array with yearly (+- pad months) time series -- before (orig = true) and after splining
% tSeries	a Mx(2*n+1) array to store the original data (even columns) and the spline interpolated (odd coluns)
% indTS		array indices of where to extract the time series (a Mx2 array)
% indTSCurr Start index to write on the tSeries vector. Needs to be incremented by 12 at the end of this function
% orig		Logical meaning that if true we are writting original series or, otherwise interpolated, where NANs, values

	if (orig),		iStart = 2;			% Remember that first column is always the time
	else,			iStart = 3;
	end
	for (k = 1:size(indTS,1))	% Loop over number of check points
		tSeries(indTSCurr:(indTSCurr+11), iStart + 2*(k-1)) = ZtoSpline(indTS(k,2), indTS(k,1), first_month:last_month);
	end
	indTSCurr = indTSCurr + 12;			% Remember, this works only for yearly means from month data

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
		else,				nc_io(grd_out, sprintf('w%d', m-1), handles, Z)
		end

		h = aguentabar(m/n_layers);
		if (isnan(h)),	break,	end
	end

% ------------------------------------------------------------------------------------------------
function [is_PFV52, times, anos] = PFV52(handles, s)
% Find if we are dealing with one Pathfinder V5.2 daily file.
% If yes and user request it, convert and return also the time vector in the form of decimal years
% In case file is not a PFV5.5, TIMES = HANDLES.TIME

	is_PFV52 = false;
	for (k = 1:numel(s.Attribute))
		if (strcmp(s.Attribute(k).Name, 'description') && ~isempty(strfind(s.Attribute(k).Value,'Pathfinder 5.2 daily')))
			is_PFV52 = true;
			break
		end
	end

	if (nargout > 1)		% This relies on the fact that the apropriate 'description' global attribute has been set
		anos = fix(handles.time(:));
		if (is_PFV52)
			n_days_in_year = datenummx(anos+1,1,1) - datenummx(anos,1,1);
			times = (handles.time - anos) * 1000;		% Decimal day of the year (wrongly aka julian day)
			times = anos + times ./ n_days_in_year;		% Now we have decimal years
		else
			times = handles.time(:);
		end
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

	opt_I = ' ';
	if (use_surface),		opt_I = sprintf('-I%.10f/%.10f',head(8),head(9));	end

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
		[Z_rect, r_c] = cropimg(head(1:2),head(3:4),Z,rect_crop,'out_grid');
		[bw_rect, l.lixo] = cropimg(head(1:2),head(3:4),bw,rect_crop,'out_grid');
		Z_rect = double(Z_rect);      % It has to be (GHRRRRRRRRRRRRR)

		%X = x_min:head(8):x_max;	Y = y_min:head(9):y_max;
		X = linspace(x_min, x_max, size(Z_rect, 2));		% Safer against round off errors 
		Y = linspace(y_min, y_max, size(Z_rect, 1));
		[XX,YY] = meshgrid(X,Y);
		XX(bw_rect) = [];			YY(bw_rect) = [];		Z_rect(bw_rect) = [];

		if (use_surface)
			opt_R = sprintf('-R%.10f/%.10f/%.10f/%.10f', X(1), X(end), Y(1), Y(end));
			%Z_rect = c_surface(XX(:), YY(:), Z_rect(:), opt_R, opt_I, '-T.25');
			Z_rect = gmtmbgrid_m(XX(:), YY(:), Z_rect(:), opt_R, opt_I, '-T.25', '-Mz');
		elseif (use_bicubic)
			Z_rect = griddata_j(XX(:), YY(:), Z_rect(:), X, Y', 'cubic');
		else
			Z_rect = griddata_j(XX(:), YY(:), Z_rect(:), X, Y', 'linear');
		end

		% Inprint the processed rectangle back into orig array
		if (isa(Z,'single')),		Z(r_c(1):r_c(2),r_c(3):r_c(4)) = single(Z_rect);
		elseif (isa(Z,'int16')),	Z(r_c(1):r_c(2),r_c(3):r_c(4)) = int16(Z_rect);
		elseif (isa(Z,'uint16')),	Z(r_c(1):r_c(2),r_c(3):r_c(4)) = uint16(Z_rect);
		else,						Z(r_c(1):r_c(2),r_c(3):r_c(4)) = single(Z_rect);
		end
	end

% ------------------------1--------2-------3------4----------5--------6---------8----
function calc_polygAVG(handles, fnameOut, op, fnamePolys, sub_set, fnameFlag, quality)
% This function search for polygons (patches or closed lines) and computes averages
% of whatever quantity is respresented inside those polygones. The result is saved
% in an ASCII file whose first two columns contain the the polygons (x,y) centroid.
%
% OPTIONS:
%
% FNAMEOUT		Name of the ouput file. If not provided, it will be asked for here.
%
% OP			Operation to apply to the data inside each polygon. If not provided or is [],
%				defaults to 'mean'. That is compute the mean of all values inside polygon.
%				Otherwise it must be the name of a Matlab function that can executed via
%				the function handles mechanism. For example 'median'.
%
% FNAMEPOLYS	A file name of a (x,y) polygon or the name of a list of polygons (one per line).
%				If given no attempt is made to fish the polygons from line handles.
%
% SUB_SET		A two columns row vec with number of the offset of years where analysis start and stop.
%				For example [3 1] Starts analysis on forth year and stops on the before last year.
%				[0 0] Means using the all dataset.
%
% FNAMEFLAG		name of a netCDF file with quality flags. Obviously this file must be of
%				the same size as the series under analysis. If not provided no quality check is done.
%
% QUALITY		Threshold quality value. Only values of quality >= FLAG will be taken into account
%				NOTE: For MODIS use negative FLAG. Than, values are retained if quality <= abs(FLAG)

	% -------------------------------- Options parsing -------------------------------------
	if (nargin < 4)
		fnamePolys = [];		sub_set = [0 0];	fnameFlag = [];
	end
	if (nargin == 1)
		fnameOut = [];			fhandle = @local_avg;
	elseif (nargin == 2)
		fhandle = @local_avg;	fnamePolys = [];
	elseif (nargin >= 3)
		if (ischar(op)),		fhandle = str2func(op);
		else,					fhandle = @local_avg;
		end
		if (nargin == 4)
			sub_set = [0 0];	fnameFlag = [];
		elseif (nargin == 5)
			fnameFlag = [];
		end
	end
	% -------------------------------- END of options parsing -----------------------------------

	if (numel(sub_set) == 2)
		jump_start = sub_set(1);		stop_before_end = sub_set(2);
	else
		jump_start = 0;					stop_before_end = 0;
	end

	[z_id, s, rows, cols] = get_ncInfos(handles);

	%------------- Check for quality flags request -------------------
	if (~isempty(fnameFlag))
		[s_flags, z_id_flags, msg] = checkFlags_compat(fnameFlag, handles.number_of_timesteps, rows, cols);
		if (~isempty(msg)),	errordlg(msg, 'Error'),		return,		end
		do_flags = true;
		if (quality > 0 ),	growing_flag = true;		% PATHFINDER flags
		else,				growing_flag = false;		% MODIS flags
		end
	else
		do_flags = false;
	end

	% --------------- Fish the polygon handles from figure (if not provided) ------------
	if (isempty(fnamePolys))

		handMir = guidata(handles.hMirFig);
		hLine = findobj(handMir.axes1,'Type','line');
		hLine = [hLine; findobj(handMir.axes1,'Type','patch')];
		if (isempty(hLine))
			errordlg('polygAVG: No polygon file provided and no polygon in fig either. Bye Bye.','Error')
			return
		end

		polys = cell(1, numel(hLine));
		N = 1;
		for (k = 1:numel(hLine))
			x = get(hLine(k),'XData');   y = get(hLine(k),'YData');
			if (numel(x) >= 3 && x(1) == x(end) && y(1) == y(end) )
				polys{N} = [x(:) y(:)];
				N = N + 1;
			end
		end
		N = N - 1;					% There was one too much incement above
		polys(N+1:end) = [];		% Remove unused

	else							% Got a filename or a list of polygons files in input
		[bin, n_column, multi_seg, n_headers] = guess_file(fnamePolys);
		if (isempty(bin))
			errordlg(['Error reading file (probaby empty)' fnamePolys],'Error'),	return
		end
		if (n_column == 1 && multi_seg == 0)			% Take it as a file names list
			fid = fopen(fnamePolys);
			c = fread(fid,'*char')';	fclose(fid);
			names = strread(c,'%s','delimiter','\n');   clear c fid;
		else
			names = {fnamePolys};
		end

		for (k = 1:numel(names))
			fname = names{k};
			if (isempty(n_headers)),    n_headers = NaN;    end
			if (multi_seg)
				polys = text_read(fname,NaN,n_headers,'>');
			else
				polys = {text_read(fname,NaN,n_headers)};
			end
		end

		N = numel(polys);
	end

	if (N == 0)
		errordlg('Fiu Fiu! No closed polygons to compute whaterver average value inside. Bye Bye.','Error')
		return
	end
	% --------------------------- END of polygons fishing section -----------------------------------

	nLayers = handles.number_of_timesteps - (jump_start + stop_before_end);	% Number of layers to be used in this run
	series_vec = (jump_start:(nLayers - 1 + jump_start)) + 1;		% Add 1 so it never starts at 0 (no good for indices)
	avg = zeros(nLayers,N) * NaN;
	THRESH = 0.25;						% Minimum percentage of valid points inside poly

	aguentabar(0,'title','Compute poligonal averages','CreateCancelBtn')

	for (m = series_vec)				% Loop over layers ensemble
		Z = nc_funs('varget', handles.fname, s.Dataset(z_id).Name, [m-1 0 0], [1 rows cols]);

		if (do_flags)
			flags = nc_funs('varget', fnameFlag, s_flags.Dataset(z_id_flags).Name, [m-1 0 0], [1 rows cols]);
			if (growing_flag),		Z(flags < quality) = NaN;	% Pathfinder style (higher the best) quality flag
			else,					Z(flags > quality) = NaN;	% MODIS style (lower the best) quality flag
			end
		end

		for (k = 1:N)				% Loop over all polygons
			x = polys{k}(:,1);			y = polys{k}(:,2);
			xp(1) = min(x);				xp(2) = max(x);
			yp(1) = min(y);				yp(2) = max(y);
			rect_crop = [xp(1) yp(1) (xp(2) - xp(1)) (yp(2) - yp(1))];
			x_lim = [xp(1) xp(2)];		y_lim = [yp(1) yp(2)];

			% Get a BoundingBox rect to save work in mask computing
			[Z_rect, l.lixo] = cropimg(handles.head(1:2),handles.head(3:4),Z,rect_crop,'out_grid');
			mask = img_fun('roipoly_j',x_lim,y_lim,Z_rect,x,y);

			% Test for a minimum of valid elements inside polygon
			zz = Z_rect(mask);
			zz = zz(:);
			ind = isnan(zz);
			if (~any(ind))
				%avg(m,k) = sum(double(zz)) / numel(zz);
				avg(m,k) = feval(fhandle, zz);
			else			% Accept/Reject based on % of valid numbers
				nAnoes = sum(ind);		nInPoly = numel(zz);
				if (nAnoes / nInPoly < THRESH)
					zz = zz(~ind);
					%avg(m,k) = sum(double(zz)) / numel(zz);
					avg(m,k) = feval(fhandle, zz);
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
	if (isempty(fnameOut))
		[FileName,PathName] = put_or_get_file(...
			handles,{'*.dat;*.DAT','ASCII file'; '*.*', 'All Files (*.*)'},'Output file','put');
		if isequal(FileName,0),		return,		end
		[PATH,FNAME,EXT] = fileparts([PathName FileName]);
		if isempty(EXT),	fname = [PathName FNAME '.dat'];
		else,				fname = [PathName FNAME EXT];
		end
	else
		fname = fnameOut;
	end
	
	% Open and write to ASCII file
	if (ispc),		fid = fopen(fname,'wt');
	elseif (isunix),fid = fopen(fname,'w');
	else,			error('aquamoto: Unknown platform.');
	end

	% Calculate a a rough polygon centroid
	centro = zeros(N,2);
	for (k = 1:N)				% Loop over polygons
		centro(k,:) = mean(polys{k});
	end

	if (strcmpi(handles.z_units, 'Decimal year'))
		fprintf(fid, '# DATENUM YYYY.XX\n');	fmt_t = '%.12f\t';
	elseif (strncmpi(handles.z_units, 'Seconds since 0000-01-01', 24))
		fprintf(fid, '# DATENUM SECONDS\n');	fmt_t = '%.3f\t';
	elseif (strncmpi(handles.z_units, 'Milliseconds since 0000-01-01', 29))
		fprintf(fid, '# DATENUM MILLISECONDS\n');	fmt_t = '%.0f\t';
	elseif (strncmpi(handles.z_units, 'UNIX', 4))
		fprintf(fid, '# DATENUM UNIX\n');		fmt_t = '%.3f\t';
	elseif (strncmpi(handles.z_units, 'RATA', 4))
		fprintf(fid, '# DATENUM RATA DIE\n');	fmt_t = '%.8f\t';
	elseif (strncmpi(handles.z_units, 'Days since 0000-01-01', 21))
		fprintf(fid, '# DATENUM DAYS\n');		fmt_t = '%.8f\t';
	else
		fmt_t = '%.5f\t';
	end

	fprintf(fid, ['#  \t', repmat('%g(X)\t', [1,N]) '\n'], centro(:,1));
	fprintf(fid, ['# T\t', repmat('%g(Y)\t', [1,N]) '\n'], centro(:,2));

	try			t = handles.time;		t = t(:);		% Layers's times
	catch,		t = (1:size(avg,1))';
	end

	fprintf(fid,[fmt_t repmat('%f\t',[1,N]) '\n'], [t avg]');
	fclose(fid);

% ----------------------------------------------------------------------
function out = local_avg(x)
% A simple average function, but that convert X do doubles and will be assigned a function handle
	out = sum(double(x)) / numel(x);

% ----------------------------------------------------------------------
function calc_flagsStats(handles, months, flag, opt)
% Compute per/pixel annual counts of pixel values with a quality >= flag
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

	n_anos = handles.number_of_timesteps / 12;

	aguentabar(0,'title','Compute flag quality counts','CreateCancelBtn')
	goodCount = zeros([rows, cols]);			% 

	if (strcmp(opt, 'per_year'))		% Make the counting on a per year basis
		for (m = 1:n_anos)
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
	
			if (m == 1),	nc_io(grd_out,sprintf('w%d/time',n_anos), handles, reshape(tmp,[1 size(tmp)]))
			else,			nc_io(grd_out, sprintf('w%d', m-1), handles, tmp)
			end
			
			h = aguentabar(m/n_anos);
			if (isnan(h)),	break,	end
			goodCount = goodCount * 0;			% Reset it to zeros
		end

	else			% Make the counting on a per month basis
		for (n = months)
			contanoes = zeros(rows, cols);
			for (m = 1:n_anos)
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
			else,			nc_io(grd_out, sprintf('w%d', n-1), handles, tmp)
			end
			
			h = aguentabar(n/numel(months));
			if (isnan(h)),	break,	end	
			goodCount = goodCount * 0;			% Reset it to zeros
		end
	end

% ----------------------------------------------------------------------
function do_math(handles, opt, subSet1, fname2, subSet2)
% Perform some basic algebraic operations on 3D planes
%
% OPT can be one of:
% 'sum'       => add all layers
% 'count'     => count the acumulated number of non-NaNs of all layers
% 'diffstd'   => Compute mean and sdt of A - B (inmemory - fname2)
%                The array read from FNAME2 does not need to be exactly comptible with A.
%                When it is not a call to grdsample will take care of compatibilization.
%
% SUBSET1 & SUBSET2 are (optional) two columns row vec with number of the offset layers
% where analysis starts and stop.
% For example [3 1] Starts analysis on forth layer and stops on the before last year.
% [0 0] Means using the all dataset (the default)
%
%	NOT FINISHED, NEED INPUT PARSING

	grid1 = [];		grid2 = [];
	s1 = handles.nc_info;				% Retrieve the .nc info struct
	rows1 = s1.Dataset(handles.netcdf_z_id).Size(end-1);
	cols1 = s1.Dataset(handles.netcdf_z_id).Size(end);
	nLayers = handles.number_of_timesteps;

	handles.geog = 1;		handles.was_int16 = 0;		handles.computed_grid = 0;

	% Find if we are dealing with a Pathfinder V5.2 daily file
	[is_PFV52, tempos] = PFV52(handles, handles.nc_info);		% if not a PFV5.2, tempos = handles.time(:)

	if (strcmp(opt, 'sum'))
		% WARNING: SUBSET1 not implemented
		grid1 = nc_funs('varget', handles.fname, s1.Dataset(handles.netcdf_z_id).Name, [0 0 0], [1 rows1 cols1]);
		is_int8 = isa(grid1, 'int8');		is_uint8 = isa(grid1, 'uint8');
		if (is_int8 || is_uint8),			grid1 = int16(grid1);		end		% To avoid ovelflows
		for (m = 2:nLayers)
			Z = nc_funs('varget', handles.fname, s1.Dataset(handles.netcdf_z_id).Name, [m-1 0 0], [1 rows1 cols1]);
			if (is_int8 || is_uint8),		Z = int16(Z);			end			% To avoid ovelflows
			cvlib_mex('add',grid1, Z);
		end
		figName1 = 'Stack Sum';

	elseif (strcmpi(opt, 'count'))			% Count the acumulated number of non NaNs along the 3rth dim of the 3D array
		
		% Check (GUESS!) if subsets were given as years instead of layer numbers
		if (is_PFV52 && subSet1(1) >= 1982 && subSet1(2) <= 2014)	% Though we only have them till 2011
			ind = find(tempos > subSet1(1));
			subSet1(1) = ind(1);
			ind = find(tempos > subSet1(2));
			subSet1(2) = nLayers - ind(1);
		end
	
		aguentabar(0,'title','Counting...','CreateCancelBtn')
		grid1 = alloc_mex(rows1, cols1, 'single');
		ini = 1+subSet1(1);		fim = nLayers-subSet1(2);	N = fim - ini + 1;	N10 = fix(N / 10);	c = 0;
		for (m = ini:fim)
			Z = nc_funs('varget', handles.fname, s1.Dataset(handles.netcdf_z_id).Name, [m-1 0 0], [1 rows1 cols1]);
			ind = isnan(Z);
			Z(~ind) = 1;	Z(ind) = 0;
			cvlib_mex('add',grid1, Z);
			p = fix(((m - ini) / N10));
			if (p)
				c = c + 1;				p = p + c - 1;		N10 = N10 + fix(N / 10);
				aguentabar(p/10)
			end
		end
		grid1(grid1 == 0) = NaN;
		aguentabar(1)
		figName1 = 'Stack Count';

	elseif (strcmpi(opt, 'diffstd'))
		if (nargin == 4),	subSet2 = [0 0];		end
		s2 = nc_funs('info', fname2);
		[X1,Y1,Z,head1] = nc_io(handles.fname, 'R');
		[X2,Y2,Z,head2,misc2] = nc_io(fname2,'R');
		rows2 = s2.Dataset(misc2.z_id).Size(end-1);
		cols2 = s2.Dataset(misc2.z_id).Size(end);
		nLayers2 = misc2.z_dim(1);
		n1 = max(1,subSet1(1)+1):(nLayers  - subSet1(2));	% Use these steps of file 1
		n2 = max(1,subSet2(1)+1):(nLayers2 - subSet2(2));	% Use these steps of file 2
		nSteps = min(numel(n1), numel(n2));					% Minimum number of steps that satisfies both files
		grid1 = alloc_mex(rows1, cols1, 'single');
		aguentabar(0,'title','Computing means and STDs')
		for (k = 0:nSteps-1)								% OK, now first compute meam of (A - B)
			tmp = nc_funs('varget', handles.fname, s2.Dataset(handles.netcdf_z_id).Name, [subSet1(1)+k 0 0], [1 rows1 cols1]);
			z2 = nc_funs('varget', fname2, s2.Dataset(misc2.z_id).Name, [subSet2(1)+k 0 0], [1 rows2 cols2]);
			z2 = make_compatible(z2, head2, head1);
			cvlib_mex('sub', tmp, z2);
			cvlib_mex('add', grid1, tmp);
			aguentabar((k+1)/(2*nSteps));
		end
		cvlib_mex('CvtScale', grid1, 1/nSteps, 0);			% Divide by N --- MEAN of (A- B)
		grid2 = alloc_mex(rows1, cols1, 'single');
		for (k = 0:nSteps-1)								% And now the STD
			tmp = nc_funs('varget', handles.fname, s1.Dataset(handles.netcdf_z_id).Name, [subSet1(1)+k 0 0], [1 rows1 cols1]);
			z2 = nc_funs('varget', fname2, s2.Dataset(misc2.z_id).Name, [subSet2(1)+k 0 0], [1 rows2 cols2]);
			z2 = make_compatible(z2, head2, head1);
			cvlib_mex('sub', tmp, z2);
			cvlib_mex('sub', tmp, grid1);					% Subtract mean
			cvlib_mex('mul', tmp, tmp);						% Square
			cvlib_mex('add', grid2, tmp);
			aguentabar(k/nSteps + 0.5);
		end
		cvlib_mex('CvtScale', grid2, 1/nSteps, 0);			% Divide by N
		cvlib_mex('pow', grid2, 0.5);						% sqrt --- STD of (A - B)
		aguentabar(1);
		clear tmp z2
		figName1 = 'Mean of differences';
		figName2 = 'STD of differences';
	end

	tmp.head = handles.head;
	if (isa(grid1,'single'))
		zz = grdutils(grid1,'-L');  tmp.head(5:6) = [zz(1) zz(2)];		% Singles & NaNs = BUGs in R13
	else
		tmp.head(5:6) = [double(min(grid1(:))) double(max(grid1(:)))];
	end
	tmp.X = linspace(tmp.head(1),tmp.head(2),cols1);
	tmp.Y = linspace(tmp.head(3),tmp.head(4),rows1);
	tmp.name = figName1;
	mirone(grid1, tmp)
	if (~isempty(grid2))
		tmp.head(5:6) = [double(min(grid2(:))) double(max(grid2(:)))];
		tmp.name = figName2;		mirone(grid2, tmp)
	end

% -----------------------------------------------------------------------------------------
function Z2 = make_compatible(Z2, head2, head1)
% Reinterpolate Z2 grid to be compatible in terms of -R and resolution with HEAD1 params
	opt_I = sprintf('-I%.8f/%.8f', head1(8:9));
	opt_R = sprintf('-R%.8f/%.8f/%.8f/%.8f', head1(1:4));
	Z2 = c_grdsample(Z2, head2, opt_R, opt_I, '-Q', '-Lg');
	
% ------------------------------------------------------------------------------------
function count_blooms(handles, Ncount, sub_set, grd_out)
% Given a 3D array with chlorophyl concentration, count events where there are at least
% NCOUNT consecutive growing elements along the 3rth dimension of the 3D array.
% For example, when cell (i,j,k) grows for 3 consecutive k, we store the value 3 (or higher
% the condition holds for k > 3), otherwise we write 0 at the cell position.
% This is a basic way of detecting blooms
%
% OPTIONS:
% NCOUNT  -> Minimum number of sucessive growing events so that number is stored in output file
%
% SUB_SET -> A two elements row vec with number of the offset of years where analysis start and stop.
%			For example [3 1] Starts analysis on forth year and stops on the before last year.
%			[0 0] Means using the all dataset.
%
% GRD_OUT	Name of the netCDF file where to store the result. Asked here if not provided.
%
% NOTE: Low memory footprint. Read laye, process it and write result

	if (nargin == 1)
		Ncount = 3;
	end
	if (nargin < 4 || (nargin == 4 && isempty(grd_out)))
		txt1 = 'netCDF grid format (*.nc,*.grd)';	txt2 = 'Select output netCDF grid';
		[FileName,PathName] = put_or_get_file(handles,{'*.nc;*.grd',txt1; '*.*', 'All Files (*.*)'},txt2,'put','.nc');
		if isequal(FileName,0),		return,		end
		grd_out = [PathName FileName];
	end

	if (nargin >= 3 && (numel(sub_set) == 2))
		jump_slices = sub_set(1);		stop_before_end_slices = sub_set(2);
	else
		jump_slices = 0;				stop_before_end_slices = 0;
	end

	n_slices = handles.number_of_timesteps;
	[z_id, s, rows, cols] = get_ncInfos(handles);
	n_slices = n_slices - (jump_slices + stop_before_end_slices);	% Number of layers to be used in this run

	aguentabar(0,'title','Counting blooms')
 	layer1 = nc_funs('varget', handles.fname, s.Dataset(z_id).Name, [jump_slices 0 0], [1 rows cols]);
 	layer2 = nc_funs('varget', handles.fname, s.Dataset(z_id).Name, [1 + jump_slices 0 0], [1 rows cols]);
	tmp = int16(layer2 > layer1);
	layer1 = layer2;

	for (m = 1:n_slices-2)		% Loop over remaining slices (first two were already consumed above)
 		layer2 = nc_funs('varget', handles.fname, s.Dataset(z_id).Name, [(m + 1 + jump_slices) 0 0], [1 rows cols]);
		ind = (layer2 > layer1);
		cvlib_mex('addS', tmp, 1)% Increase counter of consecutive growings
		tmp(~ind) = 0;			% Reset to zero those that interrupt the growing escalade
		t = tmp;
		ind = (tmp < Ncount);	% Find those that have a counting lower than the threshold.
		t(ind) = 0;				% Those that didn't (yet?) reach the threshould are written as zero

		thisLevel = handles.time(m+2);
		if (m == 1),	nc_io(grd_out, sprintf('w-%f/time',thisLevel), handles, reshape(t,[1 size(t)]))
		else,			nc_io(grd_out, sprintf('w%d\\%f', m-1, thisLevel), handles, t)
		end
		layer1 = layer2;

		h = aguentabar(m/(n_slices-2));	drawnow
		if (isnan(h)),	break,	end
	end

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
		
% ----------------------------------------------------------------------
function [s_flags, z_id_flags, msg] = checkFlags_compat(fname, number_of_timesteps, rows, cols)
% Check that the quality flags file FNAME is compatible with the other 3 input params
	msg = [];
	s_flags = nc_funs('info',fname);
	[X,Y,Z,head,misc] = nc_io(fname,'R');
	z_id_flags = misc.z_id;
	if ~(numel(head) == 9 && isfield(misc,'z_id'))
		msg = ['Blheak!! ' fname ' is is not a file with presumably with quality flags. By'];	return
	end
	if (numel(misc.z_dim) <= 2)
		msg = ['Ghrrr!! The ' fname ' is is not a 3D file. By'];		return
	end
	if (misc.z_dim(1) < number_of_timesteps)
		msg = 'Buhhuu!! The quality flags file has less "planes" than the-to-be-flagged-file. By';	return
	end
	if (~isequal([rows cols], [s_flags.Dataset(z_id_flags).Size(end-1) s_flags.Dataset(z_id_flags).Size(end)]))
		msg = 'Buhhuu!! quality flags and the-to-be-flagged-file have not the same size. By';
	end

% -------------------------------------------------------------------------------
function out = script_control(handles, auto)
% See if the OPTcontrol.txt file has an entry pointing to a file with parameters
% to run the aquaPlugin. If it has, read and parse that file.
%
% AUTO	If it's a char, than interpret it as the control script name directly
%		Any other type tells the program to search the name in OPTcontrol.txt

	if (~ischar(auto))
		opt_file = [handles.home_dir filesep 'data' filesep 'OPTcontrol.txt'];
		out = [];
		if (exist(opt_file, 'file') ~= 2),	return,		end

		fid = fopen(opt_file, 'r');
		c = (fread(fid,'*char'))';      fclose(fid);
		lines = strread(c,'%s','delimiter','\n');   clear c fid;
		fname = [];
		for (k = 1:numel(lines))
			if (~strncmp(lines{k},'MIR_AQUAPLUG',7)),	continue,	end
			fname = ddewhite(lines{k}(13:end));
			if (exist(fname,'file') ~= 2)
				errordlg(['Script file for aquaPlugin ' fname ' does not exist. Ignoring request'],'Error')
				return
			end
		end
	else
		fname = auto;
	end

	if (isempty(fname)),	return,		end		% OPTcontrol has not a MIR_AQUAPLUG key. Something will screw 

	try				% Wrap it in a try-catch so we have a chance to figure out the reason of eventual error
		fid = fopen(fname, 'r');
		if (fid < 0)
			errordlg(['File ' fname ' does not exist or can''t be open'])
			out = [];		return
		end
		c = (fread(fid,'*char'))';      fclose(fid);
		out = strread(c,'%s','delimiter','\n');   clear c fid;
		ind = true(1,numel(out));
		for (k = 1:numel(out))
			if (isempty(out{k}) || out{k}(1) == '#'),	continue,	end		% Jump comments
			ind(k) = false;
		end
		out(ind) = [];		% Remove the comment lines (if any)

		[t, r] = strtok(out{1});	% Fish the case selected in the control file.
		out{1} = str2double(r);
		for (k = 2:numel(out))
			if (strncmpi(out{k},'char',4))	% If we have a string argument, use it as it is
				[t, r] = strtok(out{k});
				out{k} = ddewhite(r);
			else
				out{k} = eval(out{k});		% Numeric arguments must be 'evaled'
			end
		end
	catch
		errordlg(lasterr, 'Errror')
		out = [];
	end

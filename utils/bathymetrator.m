function out = bathymetrator(opt_R, opt_I, fname, outfile)
% Create a binary file with the data of all files listed in the list file FNAME
% and falling inside the -R region defined by the GMT syntax OPT_R, with a bining
% size defined by the -I option stored in OPT_I.
% Files in FNAME file are assigned a hiearchy level in the second column.
% The result is written into OUTFILE.
% OUT = bathymetrator(...) returns a 3D array with the masks corresponding to each
% hierachy (actually the acumulated masks up to each level).

	if (nargin == 0)
		opt_R = '-R-9.2/-7.35/36.65/37.32';
		opt_I = '-I0.005';
		fname = 'C:\SVN\mironeWC\test_nomes_nc.txt';
		outfile = 'v:\lixo.b';
	end

	if (isempty(fname))
		hand.last_dir = cd;		% So that we can call put_or_get_file
		str1 = {'*.txt;*.TXT;*.dat;*.DAT', 'Text file (*.txt,*.TXT,*.dat,*.DAT)';'*.*', 'All Files (*.*)'};
		[FileName,PathName,hand] = put_or_get_file(hand,str1,'Select input text file name','get');
		if isequal(FileName,0),		return,		end
		fname = [PathName,FileName];
	end

	fid = fopen(fname,'r');
	todos = fread(fid,'*char');		fclose(fid);
	todos = strread(todos','%s','delimiter','\n');
	n_lines = numel(todos);
	files.level(n_lines) = 0;		files.fname = cell(n_lines,1);		c = false(n_lines,1);
	for (k = 1:n_lines)
		if (isempty(todos{k}) || todos{k}(1) == '#')
			c(k) = true;
			continue
		end
		% The string reversing is a trick to fish the last token first
		txt = todos{k}(end:-1:1);
		[t, r] = strtok(txt);
		files.level(k) = str2double(t(end:-1:1));
		files.fname{k} = ddewhite(r(end:-1:1));		% Supposedly these names can have blanks
	end
	files.level(c) = [];	files.fname(c) = [];	% Clean empty pre-allocated members

	[files, count] = find_level_groups(files);		% Count number of levels, number of files for each level and sort them
	n_levels = numel(count);

	mask  = mask_this_level(files.fname, count, 1, opt_R, opt_I);		% Mask of all files of first level
	masks = false(size(mask,1), size(mask,2), n_levels);
	old_m = mask;
	masks(:,:,1) = mask;		clear mask
	for (k = 2:n_levels)
		m = mask_this_level(files.fname, count, k, opt_R, opt_I);
		masks(:,:,k) = old_m | m;
		old_m = masks(:,:,k);		% Will be old at next iteration
	end

	%create_file(files, masks, [-9.06 36.85], [0.001 0.001], 'v:\lixo.b')
	create_file_(files, masks, opt_R, opt_I, outfile)

	if (nargout),	out = masks;	end

% -------------------------------------------------------------------------------------
function mask = mask_file(name, opt_R, opt_I)
% Make a mask out of all ensembles of the shapenc NAME file
	G = gmtmex(sprintf('grdmask "%s"?lonPointZ_1/latPointZ_1 %s %s -S0', name, opt_R, opt_I));
	mask = (G.z > 0);
	s = nc_funs('info',name);
	attribNames = {s.Attribute.Name};
	ind = strcmp(attribNames,'Number_of_main_ensembles');
	n_groups = s.Attribute(ind).Value;
	if (n_groups > 1)
		for (k = 2:n_groups)
			G = gmtmex(sprintf('grdmask "%s"?lonPointZ_%d/latPointZ_%d %s %s -S0', name, k, k, opt_R, opt_I));
			mask = mask | (G.z > 0);
		end
	end

% -------------------------------------------------------------------------------------
function mask = mask_this_level(fnames, count, lev, opt_R, opt_I)
% Make a mask out of all shapenc files having level LEV
	ind = sum(count(1:lev-1));		% Index of the LAST previously processed file
	mask = mask_file(fnames{ind+1}, opt_R, opt_I);
	for (n = 2:count(lev))			% Loop over number of the remaining files in this level
		name = fnames{ind+n};
		G = gmtmex(sprintf('grdmask "%s"?lonPointZ_1/latPointZ_1 %s %s -S0', name, opt_R, opt_I));
		mask = mask | (G.z > 0);
	end

% -------------------------------------------------------------------------------------
function [files, count] = find_level_groups(files)
% Return a COUNT N elements vector where N is the number of different levels and each element
% of the vector has the number of files with that level.
% The input FILES struct is also returned but with its members sorted by the contents of 'levels'
	n_cases = numel(files.level);			% Number of files
	[levels, IX] = sort(files.level);
	files.fname = files.fname(IX);
	[C,ia,ic] = unique(levels);				% Unique number of different levels
	files.level = ic;						% Levels are now from 1 to n-different-levels
	n_groups = numel(C);
	count = ones(n_groups, 1);
	j = 1;
	for (k = 1:n_cases-1)
		if (levels(k+1) == levels(k))
			count(j) = count(j) + 1;
		else
			j = j + 1;
		end
	end

% -------------------------------------------------------------------------------------
function create_file(files, masks, origin, inc, out_file)
% Scan all input files and copy only the points that pass through their mask level

	n_rows = size(masks, 1);	n_cols = size(masks, 2);
	fid = fopen(out_file,'wb');
	for (n = 1:numel(files.fname))
		disp(['Processing: ', files.fname{n}])
		s = nc_funs('info',files.fname{n});
		attribNames = {s.Attribute.Name};
		ind = strcmp(attribNames,'Number_of_main_ensembles');
		n_groups = s.Attribute(ind).Value;
		for (m = 1:n_groups)
			D = gmtmex(sprintf('gmtconvert "%s"?lonPointZ_%d/latPointZ_%d/z_%d', files.fname{n}, m,m,m));
			mask = masks(:,:,files.level(n));
			c = false(size(D.data, 1), 1);
			for (k = 1:size(D.data, 1))
				col = round((D.data(k,1) - origin(1)) / inc(1) -0.5) + 1;
				row = round((D.data(k,2) - origin(2)) / inc(2) -0.5) + 1;
				if (col < 1 || col > n_cols || row < 1 || row > n_rows || ~mask(row, col))
					c(k) = true;		% Out of the BB
				end
% 				if (col < 1 || col > n_cols || row < 1 || row > n_rows),	continue,	end		% Out of the BB
% 				if (mask(row, col))
% 					fwrite(fid, D.data(k,1:3), 'float32');
% 				end
			end
			D.data(c,:) = [];
			fwrite(fid, D.data', 'float32');
		end
	end
	fclose(fid);

% -------------------------------------------------------------------------------------
function create_file_(files, masks, opt_R, opt_I, out_file)
% Scan all input files and copy only the points that pass through their mask level

	[region, inc] = scan_opt_str(opt_R, opt_I);
	hdr = [region 0 1 0 inc];
	G = fill_grid_struct(masks(:,:,files.level(1)), hdr);
	append = '';		% To be used only once in gmtselect, than it will become '>'
	for (n = 1:numel(files.fname))
		disp(['Processing: ', files.fname{n}])
		s = nc_funs('info',files.fname{n});
		attribNames = {s.Attribute.Name};
		ind = strcmp(attribNames,'Number_of_main_ensembles');
		n_groups = s.Attribute(ind).Value;
		if (files.level(n) == 1)
			for (m = 1:n_groups)
				gmtmex(sprintf('gmtselect "%s"?lonPointZ_%d/latPointZ_%d/z_%d -bof %s >%s %s', files.fname{n}, m,m,m, opt_R, append, out_file))
				append = '>';
			end
		else
			for (m = 1:n_groups)
				G.z = single(masks(:,:,files.level(n)-1));
				gmtmex(sprintf('gmtselect "%s"?lonPointZ_%d/latPointZ_%d/z_%d -Ig -bof -G %s >%s %s', files.fname{n}, m,m,m, opt_R, append, out_file), G);
				append = '>';
			end
		end
	end

% -------------------------------------------------------------------------------------
function [region, inc] = scan_opt_str(opt_R, opt_I)
% Get the numeric values out of the -R & -I strings
 	ind = strfind(opt_R, '/');
 	region = [str2double(opt_R(3:ind(1)-1))  str2double(opt_R(ind(1)+1:ind(2)-1)) ...
		str2double(opt_R(ind(2)+1:ind(3)-1)) str2double(opt_R(ind(3)+1:end))];
	ind = strfind(opt_I, '/');
	if (isempty(ind))
		inc(2) = str2double(opt_I(3:end));
		inc(1) = inc(2);
	else
		inc = [str2double(opt_I(3:ind(1)-1)) str2double(opt_I(ind(1)+1:end))];
	end

% -------------------------------------------------------------------------------------------------
function G = fill_grid_struct(Z, head)
% Fill the Grid struct used in gmtmex. HEAD is the old 1x9 header vector.

	if (~isa(head, 'double')),	head = double(head);	end
	G.projection_ref_proj4 = '';
	G.projection_ref_wky = '';	
	G.range = head(1:6);
	G.inc = head(8:9);
	G.n_rows = size(Z,1);
	G.n_columns = size(Z,2);
	G.n_bands = size(Z,3);
	G.registration = head(7);
	G.no_data_value = NaN;
	G.title = '';
	G.remark = '';
	G.command = '';
	G.datatype = 'uint8';
	G.x = linspace(head(1), head(2), G.n_columns);
	G.y = linspace(head(3), head(4), G.n_rows);
	G.z = Z;
	G.x_units = '';
	G.y_units = '';
	G.z_units = '';	

% -------------------------------------------------------------------------------------------------
function list = bat_algarve()
% Example for the Algarve compilation
list = {
	'c:\j\bat\nc\algarve_LIDAR_50m.nc',	1 ...

	'c:\j\bat\nc\alvor_nearest.nc', 5 ...
	'c:\j\bat\nc\sagres_JGoncalves_5m.nc', 5 ...
	'c:\j\bat\nc\portimao_JGoncalves_5m.nc', 5 ...
	'c:\j\bat\nc\subnauta_SIMCO_5m.nc',	5 ...
	'c:\j\bat\nc\subnauta_vau_5m.nc',	5 ...
	'c:\j\bat\nc\marinertes_JG_40m.nc',	5 ...
	'c:\j\bat\nc\albufeira_JG_40m.nc',	5 ...
	'c:\j\bat\nc\quarteira_JG_40m.nc',	5 ...

	'c:\j\bat\nc\Linha_583A-591.nc',	5 ...
	'c:\j\bat\nc\Linha_592-600.nc',	5 ...
	'c:\j\bat\nc\Linha_601-608.nc',	5 ...
	'c:\j\bat\nc\Linha_609-612.nc',	5 ...

	'c:\j\bat\nc\SWIM_100m.nc',	10 ...

	'c:\j\bat\nc\CanalFaro_POLIS.nc',	15 ...

	'c:\j\bat\nc\CARTA603_IH.nc',	20 ...
	'c:\j\bat\nc\CARTA604_IH.nc',	20 ...
	'c:\j\bat\nc\CARTA605_IH.nc',	20 ...
	'c:\j\bat\nc\CARTA606_IH.nc',	20 ...

	'c:\j\bat\nc\guadiana_estuario_erwan.nc',	20 ...
	'c:\j\bat\nc\guadiana_rio_erwan.nc',	20 ...

	'c:\j\bat\nc\guadiana_EMERGE.nc',	25 ...

	'c:\j\bat\nc\portos2_IH.nc', 30 ...

	'c:\j\bat\nc\CarlosLoureiro.nc', 35 ...

	%'c:\j\bat\nc\CARTA7.nc',	40 ...
	%'c:\j\bat\nc\ATLGEO3.nc',	40 ...
	'c:\j\bat\nc\alg_bat_IH.nc', 40 ...
	};
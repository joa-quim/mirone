function chimoce(fname)
% Read a multi-segment file with two segments a write data in a NC formatted file
% First segment must contain 3 columns with station number, lon and lat. Example
% > BlaBla
% 10 -7.4 37.1167
% 11 -7.4 37.075
%
% The optional string in segment separator will be used as a Description attribute in nc file
%
% Second segment must have in segment separator the names of the variables separated by commas.
% Next comes the data in the orders specified above. Example:
% > dist, station, depth, temp, sal
% 0.0	10	-5	17.4	36.2
% 4.6	11	-5	19.7	36.3
% 9.3	12	-5	19.6	36.3
% 13.9	13	-5	19.4	36.3
% 18.5	14	-5	20.1	36.3
% 23.2	15	-5	21.1	36.4
% 27.8	16	-5	22.1	36.5
% 32.4	17	-5	22.6	36.6
% 37.1	18	-5	22.7	36.6
% 0.0	10	-10	17.3	36.2
% 4.6	11	-10	18.8	36.2
% 9.3	12	-10	19.1	36.2

	if (nargin == 0)
		fname  = 'v:\dist_m_z_T.dat';
	end

	[numeric_data, multi_segs_str] = text_read(fname,NaN,1,'>');

	% --------------- Decide the netCDF file name ---------------
	[PATO, FNAME] = fileparts(fname);
	fsep = filesep;
	if (isempty(PATO) || PATO(end) == '\'),		fsep = '';	end		% File is in the current directory
	fname = [PATO fsep FNAME '.nc'];
	% ------------------------------------------------------------

	if (numel(multi_segs_str) == 2)
		col_names = parse_col_names(multi_segs_str{2}(2:end));
		if (numel(multi_segs_str{1}) > 1)
			desc = ddewhite(multi_segs_str{1}(2:end));
		end
	else
		col_names = parse_col_names(multi_segs_str{1}(2:end));
	end

	% --- Compute accumulated distances along profile and either replace the existing 'dist' column or add a new one at the end
	lon_i = numeric_data{1}(1:end-1,2);		lat_i = numeric_data{1}(1:end-1,3);
	lon_f = numeric_data{1}(2:end,2);		lat_f = numeric_data{1}(2:end,3);
	dist  = [0; cumsum(vdist(lat_i,lon_i,lat_f,lon_f))] + eps;
	ind_dist = find(strcmp(col_names, 'dist'));
	ind_stat = strcmp(col_names, 'station');
	stations = numeric_data{2}(:,ind_stat);
	for (k = 1:size(numeric_data{1},1))
		stations(stations == stations(k)) = dist(k);
	end
	if (~isempty(ind_dist))
		numeric_data{2}(:,ind_dist) = stations;
	else
		numeric_data{2}(:,end+1) = stations;
		col_names{end+1} = 'dist';
	end
	% ----------------------------------------------------------------------------------

	nc_funs('create_empty', fname)

	% ---------------------------- Write the dimensions --------------------------------
	nc_funs('add_dimension', fname, 'depth',   size(numeric_data{2},1))
	nc_funs('add_dimension', fname, 'station', size(numeric_data{1},1))

	nc_global = -1;
	nc_funs('attput', fname, nc_global, 'Conventions', 'CF-1.0');
	nc_funs('attput', fname, nc_global, 'spatial_ref', '+proj=longlat');
	nc_funs('attput', fname, nc_global, 'ACprofiles', 'Water Chemistry');
	if (~isempty(desc)),	nc_funs('attput', fname, -1, 'Description', desc);		end

 	% ------------------------------ Write the variables ------------------------------------	
	c_names = {'st_num' 'lon' 'lat'};
	for (k = 1:3)						% First the stations coordinates
		write_var(fname, c_names{k}, 5, 'station', 'Station depth', 'meters', [], 'UTC time')
		nc_funs('varput', fname, c_names{k}, single(numeric_data{1}(:,k)), 0, size(numeric_data{1},1));
	end

	for (k = 1:numel(col_names))		% Next, the data itself
		write_var(fname, col_names{k}, 5, 'depth', 'Station depth', 'meters', [], 'UTC time', [], [])
		x = single(numeric_data{2}(:,k));
		nc_funs('varput', fname, col_names{k}, x, 0, numel(x));
	end

% ----------------------------------------------------------------------------------------------
function write_var(fname, name, tipo, dim, long_name, units, actual_range, comment, fillValue, missing_value)	
	varstruct.Name = name;
	varstruct.Nctype = tipo;
	if (~isempty(dim)),		varstruct.Dimension = {dim};	end
	nc_funs('addvar', fname, varstruct)

% ----------------------------------------------------------------------------------------------
function col_names = parse_col_names(str)
% ...
	ind = strfind(str, ',');
	col_names = cell(numel(ind) + 1, 1);

	col_names{1} = ddewhite(str(1:ind(1)-1));
	for (k = 2:numel(ind))
		col_names{k} = ddewhite(str(ind(k-1)+1:ind(k)-1));
	end
	col_names{k+1} = ddewhite(str(ind(k)+1:end));

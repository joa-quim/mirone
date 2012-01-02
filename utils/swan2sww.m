function swan2sww(fname, fname_sww)
% Create a .sww netCDF file from the files listed on the listfile FNAME
%
% FNAME_SWW, output name of the .sww netCDF file.
%
% FNAME, listfile with a list, one per record, of the grids created by Mirone-Swan
% which will be merged into the .sww file.
% An example listfile is like this
%
%sub_region_bat.grd
%tsu_time_00000.grd
%tsu_time_00003.grd
%tsu_time_00006.grd
%...
% Here the tsu_time_... grids were created by Mirone-Swan when the "momentum" option was selected
% Further from these grids this function expects also the the existence of tsu_time_?????_Uh and _Vh
% which were created by the same process.
% The file "sub_region_bat.grd" must be created by the user and its size must fit exactly
% that of the tsu_time_... grids.
%
% In the above listing filenames do not have path information. This works as long as the "listfile" file
% resides on the same directory as the grids that will be accessed. Otherwise, gridnames must
% be preceded by their full path.

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

	f_sep = filesep;
	if (nargin ~= 2),
		error('swans2sww: wrong number of arguments. I want 2')
	end
	
	nc_funs('create_empty', fname_sww)

    fid = fopen(fname);
	c = char(fread(fid))';      fclose(fid);
	names = strread(c,'%s','delimiter','\n');   clear c fid;

	% Deal with the case where names in the listfile don't have path
	pato = fileparts(names{1});
	if (isempty(pato))
		pato = fileparts(fname);
		if (~isempty(pato))
			for (k = 1:numel(names))
				names{k} = [pato filesep names{k}];
			end
		end
	end
	
	handles.grdname = names{1};				% just to have something in this "handles"
	[handles, X,Y,Z,head] = read_gmt_type_grids(handles, names{1});
	if (isempty(X)),	return,		end		% An error message has already been issued
	[ny, nx] = size(Z); 
	
	% -------------------- Dimensions
	nVolumes = (nx-1) * (ny-1) * 2;
	nPoints = numel(Z);
	nc_funs('add_dimension', fname_sww, 'number_of_volumes', nVolumes )
	nc_funs('add_dimension', fname_sww, 'number_of_vertices', 3 )
	nc_funs('add_dimension', fname_sww, 'numbers_in_range', 2 )
	nc_funs('add_dimension', fname_sww, 'number_of_points', nPoints )
	nc_funs('add_dimension', fname_sww, 'number_of_timesteps', length(names)-1 )
	
	% -------------------- Variables ------------------------------------------------
	x_varstruct.Name = 'x';		x_varstruct.Dimension = {'number_of_points'};
	y_varstruct.Name = 'y';		y_varstruct.Dimension = {'number_of_points'};
	nc_funs('addvar', fname_sww, x_varstruct)
	nc_funs('addvar', fname_sww, y_varstruct)

	varstruct.Dimension = {'number_of_points'};
	varstruct.Name = 'z';		varstruct.Nctype = 6;
	nc_funs('addvar', fname_sww, varstruct )
	varstruct.Name = 'elevation';		varstruct.Nctype = 6;
	nc_funs('addvar', fname_sww, varstruct )
	varstruct.Dimension = {'numbers_in_range'};
	varstruct.Name = 'elevation_range';
	nc_funs('addvar', fname_sww, varstruct )

	varstruct.Name = 'volumes';		varstruct.Nctype = 4;
	varstruct.Dimension = {'number_of_volumes' 'number_of_vertices'};
	nc_funs('addvar', fname_sww, varstruct )
	varstruct.Name = 'time';		varstruct.Nctype = 6;
	varstruct.Dimension = {'number_of_timesteps'};
	nc_funs('addvar', fname_sww, varstruct )
	varstruct.Name = 'stage';		varstruct.Dimension = {'number_of_timesteps' 'number_of_points'};
	nc_funs('addvar', fname_sww, varstruct )
	varstruct_range.Dimension = {'numbers_in_range'};
	varstruct_range.Nctype = 6;
	varstruct_range.Name = 'stage_range';			nc_funs('addvar', fname_sww, varstruct_range )
	varstruct.Name = 'xmomentum';					nc_funs('addvar', fname_sww, varstruct )
	varstruct_range.Name = 'xmomentum_range';		nc_funs('addvar', fname_sww, varstruct_range )
	varstruct.Name = 'ymomentum';					nc_funs('addvar', fname_sww, varstruct )
	varstruct_range.Name = 'ymomentum_range';		nc_funs('addvar', fname_sww, varstruct_range )
	% --------------------------------------------------------------------------------

	% -------------------------- Globals ---------------------------------------------
	nc_global = -1;
	nc_funs('attput', fname_sww, nc_global, 'institution', 'Mirone Tec' );
	nc_funs('attput', fname_sww, nc_global, 'description', 'Converted from Mirone-Swan' );
	nc_funs('attput', fname_sww, nc_global, 'xllcorner', head(1) );
	nc_funs('attput', fname_sww, nc_global, 'yllcorner', head(3) );
	nc_funs('attput', fname_sww, nc_global, 'zone', 29 );
	nc_funs('attput', fname_sww, nc_global, 'starttime', 0 );
	nc_funs('attput', fname_sww, nc_global, 'false_easting', 500000 );
	nc_funs('attput', fname_sww, nc_global, 'false_northing', 0 );
	nc_funs('attput', fname_sww, nc_global, 'datum', 'wgs84' );
	nc_funs('attput', fname_sww, nc_global, 'projection', 'UTM' );
	nc_funs('attput', fname_sww, nc_global, 'units', 'm' );
	% --------------------------------------------------------------------------------

	% Put the coords vectors
	X = X(:) - X(1);
	Y = flipud(Y(:) - Y(1));
	[X, Y] = meshgrid(X,Y);
	Z = flipud(Z);
	Z = Z(:)';
	nc_funs('varput', fname_sww, 'x', X(:) );
	nc_funs('varput', fname_sww, 'y', Y(:) );
	nc_funs('varput', fname_sww, 'z', Z );
	nc_funs('varput', fname_sww, 'elevation', Z );
	nc_funs('varput', fname_sww, 'elevation_range', head(5:6) );
	clear X Y Z;

	vertices = (reshape(0:(nx*ny-1),ny,nx));
	volumes = int32( zeros(nVolumes, 3) );
	i = 1;
	for (n = 1:nx-1)
		for (m = 1:ny-1)
			v1 = vertices(m,n);
			v2 = vertices(m+1,n);
			v3 = vertices(m+1,n+1);
			v4 = vertices(m,n+1);
			volumes(i,:) = [v1,v2,v3];		i = i + 1;
			volumes(i,:) = [v1,v3,v4];		i = i + 1;
		end
	end

	% Put the volumes variable
	nc_funs('varput', fname_sww, 'volumes', volumes, [0 0], [nVolumes 3] );
	clear volumes;

	% --------- Do a first round to get the time vector
	names(1) = [];				% It was the bathymetry file
	nTimes = length(names);
	tempo = zeros(1, nTimes);
	for (k = 1:nTimes)
		id = strfind(names{k}, '_time_');
		tempo(k) = str2double( names{k}(id+6:id+10) );
	end
	
	tempo = tempo - tempo(1);
	
	nc_funs('varput', fname_sww, 'time', tempo );

	% --------- Now do the hardcore writing
	% ----------------------------- STAGE -----------------------------------------
	minmax = [inf -inf];
	for (k = 1:nTimes)
		[pato, nome, ext] = fileparts(names{k});
		if (isempty(ext)),		thisName = [names{k} '.grd'];
		else					thisName = names{k};
		end
		if (exist(thisName,'file') ~= 2)
			disp(['swan2sww: Warning, could not find file ' thisName])
			continue
		end
		[handles,X,Y,Z,head] = read_gmt_type_grids(handles, thisName);		% Surface elevation
		if ( any(isnan(Z)) )
			disp(['swan2sww: Warning, file ' thisName ' has NaNs, which is deadly mortal to ANUGA'])
		end
		Z = flipud(Z);
		nc_funs('varput', fname_sww, 'stage', Z(:)', [k-1 0], [1 nPoints] );
		minmax(1) = min(minmax(1), head(5));		% Compute global min/max
		minmax(2) = max(minmax(2), head(6));
	end
	nc_funs('varput', fname_sww, 'stage_range', minmax(1:2) );

	% ----------------------------- XMOMENTUM --------------------------------------
	minmax = [inf -inf];
	for (k = 1:nTimes)
		[pato, nome, ext] = fileparts(names{k});
		if (isempty(ext)),		thisName = [names{k} '_Uh.grd'];
		else					thisName = [pato f_sep nome '_Uh' ext];
		end
		if (exist(thisName,'file') ~= 2)
			disp(['swan2sww: Warning, could not find file ' thisName])
			continue
		end
		[handles,X,Y,Z,head] = read_gmt_type_grids(handles, thisName);
		if ( any(isnan(Z)) )
			disp(['swan2sww: Warning, file ' thisName ' has NaNs, which is deadly mortal to ANUGA'])
		end
		Z = flipud(Z);
		nc_funs('varput', fname_sww, 'xmomentum', Z(:)', [k-1 0], [1 nPoints] );
		minmax(1) = min(minmax(1), head(5));		% Compute global min/max
		minmax(2) = max(minmax(2), head(6));
	end
	nc_funs('varput', fname_sww, 'xmomentum_range', minmax(1:2) );

	% ----------------------------- YMOMENTUM --------------------------------------
	minmax = [inf -inf];
	for (k = 1:nTimes)
		[pato, nome, ext] = fileparts(names{k});
		if (isempty(ext)),		thisName = [names{k} '_Vh.grd'];
		else					thisName = [pato f_sep nome '_Vh' ext];
		end
		if (exist(thisName,'file') ~= 2)
			disp(['swan2sww: Warning, could not find file ' thisName])
			continue
		end
		[handles,X,Y,Z,head] = read_gmt_type_grids(handles, thisName);
		if ( any(isnan(Z)) )
			disp(['swan2sww: Warning, file ' thisName ' has NaNs, which is deadly mortal to ANUGA'])
		end
		Z = flipud(Z);
		nc_funs('varput', fname_sww, 'ymomentum', Z(:)', [k-1 0], [1 nPoints] );
		minmax(1) = min(minmax(1), head(5));		% Compute global min/max
		minmax(2) = max(minmax(2), head(6));
	end
	nc_funs('varput', fname_sww, 'ymomentum_range', minmax(1:2) );

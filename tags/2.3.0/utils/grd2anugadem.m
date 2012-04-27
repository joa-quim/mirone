function grd2anugadem(grd_in, grd_out)
% Convert a grid (typically a COARDS or CF compliant, or any of the GMT formats) into the flavor that ANUGA apreciates
% 
% GRD_IN, input grid in a format recognized by the read_gmt_type_grids() function
% GRD_OUT, name stem for the to be created .dem netCDF file.
% Note: You must edit this file right below if the hardwired defaults info
%		is not at your own taste.

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
	
	% ------- Decide now what to put in some of the Global Attribs ----------
	% -- This doesn't have to be strickly true but the guy (anuga) must be 
	% -- fooled to beleave that it's dealing with UTM grids.
	institution = 'Mirone Tec';
	desc = 'netCDF non-compliant Converted from a GMT grid';
	NODATA_value = nan;			% It shouls always be true because nc_io or grdread_m take care but ...  
	UTMzone = 29;
	false_easting = 500000;
	false_northing = 0;
	projection = 'UTM';
	units = 'METERS';
	
	Nctype = 6;			% This is redundant for 6 (double) is the default, but we can change it
	% -----------------------------------------------------------------------

	if (~( nargin >= 2 && (ischar(grd_in) && ischar(grd_out)) ))
		error('grd2anugadem: wants two file names as input arguments')
	end
	
	[handles, X, Y, Z, head] = read_gmt_type_grids([],grd_in);

    % Make a good guess if coords are geographic
    geog = double( ( (head(1) >= -180 && head(2) <= 180) || (head(1) >= 0 && head(2) <= 360) )...
        && (head(3) >= -90 && head(4) <= 90) );
	if (geog)
		warndlg('grd2anugadem: At time of this writting (25-10-2007) ANUGA couldn''t deal with geogs grids.','WARNING')
	end
	if ( abs( diff(head(8:9)) ) > 1e-5 )
		warndlg('grd2anugadem: grid cells are not squared. They should and I don''t know the effect of it.','WARNING')
	end
	
	% Make sure grd_out has .dem extension
	[pato,name] = fileparts(grd_out);
	if (~isempty(pato)),	pato = [pato filesep];		end
	grd_out = [pato name '.dem'];
	
	nc_funs('create_empty', grd_out)
	
	% ----------------------- Write dimensions -------------------------------
	nx = round(diff(head(1:2)) / head(8) + ~head(7));
	ny = round(diff(head(3:4)) / head(9) + ~head(7));
	nc_funs('add_dimension', grd_out, 'number_of_rows', ny )
	nc_funs('add_dimension', grd_out, 'number_of_columns', nx )

	% ----------------------- Write variable ---------------------------------
	varstruct.Name = 'elevation';
	varstruct.Dimension = {'number_of_rows', 'number_of_columns'};
	varstruct.Nctype = Nctype;

	nc_funs('addvar', grd_out, varstruct)
	nc_funs('varput', grd_out, 'elevation', Z, [0 0], [ny nx] );

	% ----------------------- Write Globals ----------------------------------
	nc_global = -1;
	nc_funs('attput', grd_out, nc_global, 'institution', institution );
	nc_funs('attput', grd_out, nc_global, 'description', desc );
	nc_funs('attput', grd_out, nc_global, 'ncols', int32(nx) );
	nc_funs('attput', grd_out, nc_global, 'nrows', int32(ny) );
	if (head(7))		% Pixel reg grid
		xllcorner = head(1);		yllcorner = head(3);
	else				% Normal node reg grid
		xllcorner = head(1) - head(8) / 2;
		yllcorner = head(3) - head(9) / 2;
	end
	nc_funs('attput', grd_out, nc_global, 'xllcorner', xllcorner );
	nc_funs('attput', grd_out, nc_global, 'yllcorner', yllcorner );
	nc_funs('attput', grd_out, nc_global, 'cellsize', head(8) );
	nc_funs('attput', grd_out, nc_global, 'NODATA_value', NODATA_value );
	nc_funs('attput', grd_out, nc_global, 'zone', UTMzone );
	nc_funs('attput', grd_out, nc_global, 'false_easting', false_easting );
	nc_funs('attput', grd_out, nc_global, 'false_northing', false_northing );
	nc_funs('attput', grd_out, nc_global, 'projection', projection );
	nc_funs('attput', grd_out, nc_global, 'datum', 'WGS84' );
	nc_funs('attput', grd_out, nc_global, 'units', units );

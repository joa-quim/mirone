function conv2anugapts(fname_in, fname_out, opt)
% Convert input data into another flavor of netCDF file that pleases ANUGA: the PTS format
% PTS format files are just x, y, z points stored in a certain way
% 
% FNAME_IN, input grid in a format recognized by the read_gmt_type_grids() function
% FNAME_OUT, name stem for the to be created .pts netCDF file.
% OPT -> not used yet
%
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
	desc = 'NetCDF pts (x,y,z) format (ANUGA)';
	UTMzone = 29;
	false_easting = 500000;
	false_northing = 0;
	projection = 'UTM';
	units = 'METERS';

	Nctype = 6;			% This is redundant for 6 (double) is the default, but we can try to change it
	% -----------------------------------------------------------------------

	% ROW_MAJOR selects if data is stored with C ordering (row_major) or the
	% unfortunately still crawling column_major order (FORTRAN & MATLAB)
	row_major = true;

	% Testing argins
	if (nargin < 3),	opt = [];	end
	if (~( nargin >= 2 && (ischar(fname_in) && ischar(fname_out)) ))
		error('conv2anugapts: wants two file names as input arguments')
	end

	% ...
	if (nargin == 2 || opt(1) == 'g')
		try
			[handles, X, Y, Z, head] = read_gmt_type_grids([],fname_in);
		catch
			errordlg(['conv2anugapts: something screw up while reading grid -> ' lasterr],'ERROR');
			return
		end

        % Make a good guess if coords are geographic
        geog = double( ( (head(1) >= -180 && head(2) <= 180) || (head(1) >= 0 && head(2) <= 360) )...
            && (head(3) >= -90 && head(4) <= 90) );
		if (geog)
			warndlg('conv2anugapts: At time of this writting (25-10-2007) ANUGA couldn''t deal with geog data.','WARNING')
		end
		
		[nrows, ncols] = size(Z);					% I think those are irrelevant but ...

		% Shift origin to min(X), min(Y)
% 		X = X - min(X);		Y = Y - min(Y);
		
		if (row_major)
	 		Y = Y(end:-1:1);
			% minimalist meshgrid - singles frendly, which therefore alows saving LOTS of memory
			[Y, X] = meshgrid_free(Y, X);
			Z = (flipud(Z))';
		else
			% minimalist meshgrid - singles frendly, which therefore alows saving LOTS of memory
			[X, Y] = meshgrid_free(X, Y);
		end
			
		% This doesn't consume extra memory
		X = X(:);		Y = Y(:);		Z = Z(:);
		
		% If there are NaNs, drop them
		ind = isnan(Z);
		if (any(ind))
			X(ind) = [];		Y(ind) = [];	Z(ind) = [];
		end
		number_of_points = numel(Z);
		
		% Here I'm lost. _dem2pts uses xllcorner & yllcorner to compute X & Y as if they were grid registerd
		% e.g. x = j*cellsize + xllcorner
		% but since xllcorner & yllcorner come from an Arc grid they are pixel registerd so the _dem2pts
		% computation is wrong. Here I don't do this error, but will be the use of xll & yll?
		% Should I correct them to be grid reg?

		% Without pix/grid correction untill further knowledge
% 		if (head(7))		% Pixel reg grid
% 			xllcorner = head(1);		yllcorner = head(3);
% 		else				% Normal node reg grid
% 			xllcorner = head(1) - head(8) / 2;
% 			yllcorner = head(3) - head(9) / 2;
% 		end

		% Use absolute referencing (uff, this probably saves from the pixel reg bug)
		xllcorner = 0;
		yllcorner = 0;
		
	else
		error('conv2anugapts: other than input grid is not yet programed')
	end
	
	% Make sure fname_out has .pts extension
	[pato,name] = fileparts(fname_out);
	if (~isempty(pato)),	pato = [pato filesep];		end
	fname_out = [pato name '.pts'];
	
	nc_funs('create_empty', fname_out)
	
	% ----------------------- Write dimensions -------------------------------
	nc_funs('add_dimension', fname_out, 'number_of_points', number_of_points )
	nc_funs('add_dimension', fname_out, 'number_of_dimensions', 2 )

	% ----------------------- Create variables -------------------------------
	varstruct.Name = 'points';
	varstruct.Dimension = {'number_of_points', 'number_of_dimensions'};
	varstruct.Nctype = Nctype;
	nc_funs('addvar', fname_out, varstruct)
	varstruct.Name = 'elevation';
	varstruct.Dimension = {'number_of_points'};
	nc_funs('addvar', fname_out, varstruct)

	% ----------------------- Write Globals ----------------------------------
	nc_global = -1;
	nc_funs('attput', fname_out, nc_global, 'institution', institution );
	nc_funs('attput', fname_out, nc_global, 'description', desc );
	nc_funs('attput', fname_out, nc_global, 'ncols', int32(ncols) );		% Does this matter?
	nc_funs('attput', fname_out, nc_global, 'nrows', int32(nrows) );
	nc_funs('attput', fname_out, nc_global, 'xllcorner', xllcorner );
	nc_funs('attput', fname_out, nc_global, 'yllcorner', yllcorner );
	nc_funs('attput', fname_out, nc_global, 'zone', int32(UTMzone) );
	nc_funs('attput', fname_out, nc_global, 'false_easting', false_easting );
	nc_funs('attput', fname_out, nc_global, 'false_northing', false_northing );
	nc_funs('attput', fname_out, nc_global, 'projection', projection );
	nc_funs('attput', fname_out, nc_global, 'datum', 'WGS84' );
	nc_funs('attput', fname_out, nc_global, 'units', units );

	% ----------------------- Write variables ------------------------------------
	if (row_major),
		xy = [X Y];
	else
		xy = [Y X];
	end
	nc_funs('varput', fname_out, 'points', xy, [0 0], [number_of_points 2] );
	nc_funs('varput', fname_out, 'elevation', Z );

% --------------------------------------------------------------------------------	
function [X, Y] = meshgrid_free(X, Y)
	% A minimalist meshgrid that works with singles and therefore saves LOTS of memory
	
	X = X(:).'; 		% Make sure x is a full row vector.
	Y = Y(:);   		% Make sure y is a full column vector.
	nx = numel(X);
	X = X(ones(numel(Y), 1),:);
	Y = Y(:,ones(1, nx));
	

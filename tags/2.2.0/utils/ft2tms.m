function ft2tms(fname, fname_tms)
% Convert a MxN maregraphs file with [time wave_height] into the netCDF .tms format for ANUGA
% 
% FNAME, file name of the maregraphs file (normally one produced by Mirone-swan)
% FNAME_TMS, name stem for the to be created .tms netCDF files.
% Each column of wave heights will be written as an individual .tms file with the
% column number appended to the file name. Example "..../mareg_03.tms"

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

	if (~( nargin == 2 && (ischar(fname) && ischar(fname_tms)) ))
		error('ft2tms: wants two file names as input arguments')
	end
	
	data = text_read(fname);
	if (isempty(data))
		error('ft2tms: input data file is either empty or has nothing valuable.')
	end
	if (iscell(data))
		data = cat(1,data{:});
	end

	[pato, name] = fileparts(fname_tms);
	if (~isempty(pato)),	pato = [pato filesep];		end

	nCol = size(data,2) - 1;
	s_format = sprintf('_%%.%dd.tms', fix(log10(nCol))+1);		% Format string to have _1.tms, _2, or _01, ... _99, etc
	for (k = 1:nCol)
		thisName = [pato filesep name sprintf(s_format, k)];
		do_conversion(data(:,[1 k+1]), thisName)
	end
	
% ----------------------------------------------------------------------------------
function do_conversion(data, fname_tms)

	nc_funs('create_empty', fname_tms)
	
	% -------------------- Dimensions ----------------------------------------------
	number_of_timesteps = size(data,1);
	nc_funs('add_dimension', fname_tms, 'number_of_timesteps', number_of_timesteps )
	
	% -------------------- Variables ------------------------------------------------
	varstruct.Name = 'time';		varstruct.Dimension = {'number_of_timesteps'};
	nc_funs('addvar', fname_tms, varstruct)
	varstruct.Name = 'stage';
	nc_funs('addvar', fname_tms, varstruct)
	varstruct.Name = 'xmomentum';
	nc_funs('addvar', fname_tms, varstruct)
	varstruct.Name = 'ymomentum';
	nc_funs('addvar', fname_tms, varstruct )
	% --------------------------------------------------------------------------------

	% -------------------------- Globals ---------------------------------------------
	nc_global = -1;
	nc_funs('attput', fname_tms, nc_global, 'institution', 'Mirone Tec' );
	nc_funs('attput', fname_tms, nc_global, 'description', 'Converted from Mirone-Swan' );
	nc_funs('attput', fname_tms, nc_global, 'starttime', 0 );
	% --------------------------------------------------------------------------------

	% -------------------------- Put the Variables vectors ---------------------------
	nc_funs('varput', fname_tms, 'time', data(:,1)' );
	nc_funs('varput', fname_tms, 'stage', data(:,2)' );
	nc_funs('varput', fname_tms, 'xmomentum', zeros(1,number_of_timesteps) );
	nc_funs('varput', fname_tms, 'ymomentum', zeros(1,number_of_timesteps) );

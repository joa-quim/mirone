function gmt2mgd77_plus(fname, varargin)
% Convert old style .gmt files into the new MGD77+ format
%
% gmt2mgd77_plus(FNAME [,'name', NAME_NC, 'anom', 1|offset, 'meta', META_STRUCT])
%	FNAME is the name of the .gmt file (with or without the .gmt extension)
%		Between backets is a list of optional inputs described below.
%
%		P			V
%	'name'		Output .nc file name. If PV is omited the output name is obtained by replacing the .gmt with .nc
%
%	'anom'		1|F_offset		A Value of 1 means that the .gmt file has a magnetic anomaly.
%								A Value different from one is interpreted as meaning that the .gmt file has the
%								total field stored, to which a constant (F_offset) has been subtracted.
%								If this PV is not provided, an F_offset 40000 is assumed.
%
%	'swap'		whatever		Swap bytes of the .gmt file. Precious to recover old files
%
%	'meta'		A metadata structure.		Fileds of this structure are:
%								country			(Default 'Coxichina')
%								funding			(Default 'WWF')
%								chief			(Default 'Mestre Silva')
%								ship			(Default 'Nau Catarineta')
%								leg				(Default 'Chicken Leg')
%								port_departure	(Default 'Port Wine')
%								port_arrival	(Default 'Green Wine')
%								tow_dist		(Default '199')
%								survey_id		(Default '00000000')
%								DC_file_number	(Default '00000000')

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

	F_offset = 40000;		name_out = [];		meta = [];		is_anom = false;
	fsep = filesep;			opt_Y = '';

	for (k = 1:2:numel(varargin))
		if (strcmp(varargin{k}, 'name'))
			name_out = varargin{k+1};
		elseif (strcmp(varargin{k}, 'anom'))
			is_anom = varargin{k+1};
			if (is_anom == 1),		F_offset = 0;
			else					F_offset = is_anom;
			end
		elseif (strcmp(varargin{k}, 'meta'))
			meta = varargin{k+1};
		elseif (strcmp(varargin{k}, 'swap'))
			opt_Y = '-Y';
		end
	end

	% --------------- Default or input metadata --------------
	att.country = 'Coxichina';
	att.funding = 'WWF';
	att.chief = 'Mestre Silva';
	att.ship = 'Nau Catarineta';
	att.leg = 'Chicken Leg';
	att.port_departure = 'Port Wine';
	att.port_arrival = 'Green Wine';
	att.tow_dist = '199';
	att.survey_id = '00000000';
	att.DC_file_number = '00000000';
	if (~isempty(meta))		% get meta data
		if (isfield(meta, 'country')),			att.country = meta.country;		end
		if (isfield(meta, 'funding')),			att.funding = meta.funding;		end
		if (isfield(meta, 'chief')),			att.chief = meta.chief;			end
		if (isfield(meta, 'ship')),				att.ship = meta.ship;			end
		if (isfield(meta, 'leg')),				att.leg = meta.leg;				end
		if (isfield(meta, 'tow_dist')),			att.tow_dist = meta.tow_dist;	end
		if (isfield(meta, 'survey_id')),		att.survey_id = meta.survey_id;	end
		if (isfield(meta, 'port_departure')),	att.port_departure = meta.port_departure;	end
		if (isfield(meta, 'port_arrival')),		att.port_arrival = meta.port_arrival;		end
		if (isfield(meta, 'DC_file_number')),	att.DC_file_number = meta.DC_file_number;	end
	end
	if (~isempty(meta) && ~isa(meta.tow_dist,'char'))			% I made this mistake once already
		att.tow_dist = sprintf('%d', att.tow_dist);
	end
	% ------------------------------------------------------------
% 	att.country = 'Portugal';
% 	att.funding = 'EMEPC';
% 	att.chief = 'Marc Andre Gutcher';
% 	att.ship = 'D. Carlos I';
% 	att.leg = 'Chicken Leg';
% 	att.port_departure = 'Lisboa';
% 	att.port_arrival = 'Lisboa';
% 	att.tow_dist = '1';

	[PATO, FNAME] = fileparts(fname);
	if (isempty(PATO)),		fsep = [];	end		% File is in the current directory
	if (isempty(opt_Y))
		track = c_gmtlist([PATO fsep FNAME], '-Fsxygmt', '-G');
	else
		track = c_gmtlist([PATO fsep FNAME], '-Fsxygmt', '-G', opt_Y);
	end

	if (isempty(track.time))
		disp(['File ' FNAME ' is empty'])
		return
	end

	track.time = track.time + (date2jd(track.year) - date2jd(1970)) * 86400;	% Here we need time in seconds since 1970
	
	att.source_institution = track.agency;
	if (double(att.source_institution(end)) < 32),		att.source_institution(end) = ' ';		end		% Prevent trash
	tempo = clock;
	att.yearNow = tempo(1);		att.monthNow = tempo(2);	att.dayNow = tempo(3);
	
	% --------------- Decide the netCDF file name ---------------
	if (~isempty(name_out))
		fname = name_out;
	else
		fname = [PATO fsep FNAME '.nc'];
	end
	% ------------------------------------------------------------

	nc_funs('create_empty', fname)

	% ---------------------------- Write the dimensions --------------------------------
	nc_funs('add_dimension', fname, 'time', 0 )
	nc_funs('add_dimension', fname, 'id_dim', 8 )
	nc_funs('add_dimension', fname, 'sln_dim', 5 )
	nc_funs('add_dimension', fname, 'sspn_dim', 6 )
	% ----------------------------------------------------------------

	% -------------------------- Globals --------------------------	
	min_lon = track.info(5);		max_lon = track.info(6);		% Lon / Lat min/max
	min_lat = track.info(7);		max_lat = track.info(8);

	nc_global = -1;
	nc_funs('attput', fname, nc_global, 'Conventions', 'CF-1.0' );
	nc_funs('attput', fname, nc_global, 'Version', '2009.03.20' );
	nc_funs('attput', fname, nc_global, 'Author', 'Me' );
	nc_funs('attput', fname, nc_global, 'title', 'Cruise ' );
	nc_funs('attput', fname, nc_global, 'history', 'Conversion from .gmt to MGD77+ netCDF format' );
	nc_funs('attput', fname, nc_global, 'Survey_Identifier', att.survey_id);
	nc_funs('attput', fname, nc_global, 'Format_Acronym', 'MGD77' );
	nc_funs('attput', fname, nc_global, 'Data_Center_File_Number', att.DC_file_number );
	nc_funs('attput', fname, nc_global, 'Parameters_Surveyed_Code', '');
	nc_funs('attput', fname, nc_global, 'File_Creation_Year', num2str(att.yearNow) );
	nc_funs('attput', fname, nc_global, 'File_Creation_Month', num2str(att.monthNow) );
	nc_funs('attput', fname, nc_global, 'File_Creation_Day', num2str(att.dayNow) );
	nc_funs('attput', fname, nc_global, 'Source_Institution', att.source_institution);
	nc_funs('attput', fname, nc_global, 'Country', att.country);
	nc_funs('attput', fname, nc_global, 'Platform_Name', att.ship);
	nc_funs('attput', fname, nc_global, 'Platform_Type_Code', '1');
	nc_funs('attput', fname, nc_global, 'Platform_Type', 'SHIP');
	nc_funs('attput', fname, nc_global, 'Chief_Scientist', att.chief);
	nc_funs('attput', fname, nc_global, 'Project_Cruise_Leg', att.leg);
	nc_funs('attput', fname, nc_global, 'Funding', att.funding);
	nc_funs('attput', fname, nc_global, 'Survey_Departure_Year', num2str(track.info(11)) );
	nc_funs('attput', fname, nc_global, 'Survey_Departure_Month', num2str(track.info(10)) );
	nc_funs('attput', fname, nc_global, 'Survey_Departure_Day',  num2str(track.info(9)) );
	nc_funs('attput', fname, nc_global, 'Port_of_Departure', att.port_departure);
	nc_funs('attput', fname, nc_global, 'Survey_Arrival_Year', num2str(track.info(14)) );
	nc_funs('attput', fname, nc_global, 'Survey_Arrival_Month', num2str(track.info(13)) );
	nc_funs('attput', fname, nc_global, 'Survey_Arrival_Day', num2str(track.info(12)) );
	nc_funs('attput', fname, nc_global, 'Port_of_Arrival', att.port_arrival);
	nc_funs('attput', fname, nc_global, 'Navigation_Instrumentation', '');
	nc_funs('attput', fname, nc_global, 'Geodetic_Datum_Position_Determination_Method', '');
	nc_funs('attput', fname, nc_global, 'Bathymetry_Instrumentation', '');
	nc_funs('attput', fname, nc_global, 'Bathymetry_Add_Forms_of_Data', '');
	nc_funs('attput', fname, nc_global, 'Magnetics_Instrumentation', 'PROTON PRECESSION MAG');
	nc_funs('attput', fname, nc_global, 'Magnetics_Add_Forms_of_Data', '');
	nc_funs('attput', fname, nc_global, 'Gravity_Instrumentation', '');
	nc_funs('attput', fname, nc_global, 'Gravity_Add_Forms_of_Data', '');
	nc_funs('attput', fname, nc_global, 'Seismic_Instrumentation', '');
	nc_funs('attput', fname, nc_global, 'Seismic_Data_Formats', '');
	nc_funs('attput', fname, nc_global, 'Format_Type', '');
	nc_funs('attput', fname, nc_global, 'Format_Description', '(I1,A8,I3,I4,3I2,F5.3,F8.5,F9.5,I1,F6.4,F6.1,I2,i1,3F6.1,I1,F5.1,F6.0,F7.1,F6.1,F5.1,A5,A6,I1)');
	nc_funs('attput', fname, nc_global, 'Topmost_Latitude', num2str(ceil(max_lat)) );
	nc_funs('attput', fname, nc_global, 'Bottommost_Latitude', num2str(floor(max_lat)) );
	nc_funs('attput', fname, nc_global, 'Leftmost_Longitude', num2str(floor(min_lon)) );
	nc_funs('attput', fname, nc_global, 'Rightmost_Longitude', num2str(ceil(max_lon)) );
	nc_funs('attput', fname, nc_global, 'Bathymetry_Digitizing_Rate', '');
	nc_funs('attput', fname, nc_global, 'Bathymetry_Sampling_Rate', '');
	nc_funs('attput', fname, nc_global, 'Bathymetry_Assumed_Sound_Velocity', '15000');
	nc_funs('attput', fname, nc_global, 'Bathymetry_Datum_Code', '');
	nc_funs('attput', fname, nc_global, 'Bathymetry_Interpolation_Scheme', '');
	nc_funs('attput', fname, nc_global, 'Magnetics_Digitizing_Rate', '');
	nc_funs('attput', fname, nc_global, 'Magnetics_Sampling_Rate', '');
	nc_funs('attput', fname, nc_global, 'Magnetics_Sensor_Tow_Distance', att.tow_dist);
	nc_funs('attput', fname, nc_global, 'Magnetics_Sensor_Depth', '00010');
	nc_funs('attput', fname, nc_global, 'Magnetics_Sensor_Separation', '');
	nc_funs('attput', fname, nc_global, 'Magnetics_Ref_Field_Code', '');
	nc_funs('attput', fname, nc_global, 'Magnetics_Ref_Field', '');
	nc_funs('attput', fname, nc_global, 'Magnetics_Method_Applying_Res_Field', '');
	nc_funs('attput', fname, nc_global, 'Gravity_Digitizing_Rate', '');
	nc_funs('attput', fname, nc_global, 'Gravity_Sampling_Rate', '');
	nc_funs('attput', fname, nc_global, 'Gravity_Theoretical_Formula_Code', '');
	nc_funs('attput', fname, nc_global, 'Gravity_Theoretical_Formula', '');
	nc_funs('attput', fname, nc_global, 'Gravity_Reference_System_Code', '');
	nc_funs('attput', fname, nc_global, 'Gravity_Reference_System', '');
	nc_funs('attput', fname, nc_global, 'Gravity_Corrections_Applied', '');
	nc_funs('attput', fname, nc_global, 'Gravity_Departure_Base_Station', '');
	nc_funs('attput', fname, nc_global, 'Gravity_Departure_Base_Station_Name', '');
	nc_funs('attput', fname, nc_global, 'Gravity_Arrival_Base_Station', '');
	nc_funs('attput', fname, nc_global, 'Gravity_Arrival_Base_Station_Name', '');
	nc_funs('attput', fname, nc_global, 'Number_of_Ten_Degree_Identifiers', '5');
	nc_funs('attput', fname, nc_global, 'Ten_Degree_Identifier', '7302,7303,7400,7401,7402,9999,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,');
	nc_funs('attput', fname, nc_global, 'Additional_Documentation_1', '');
	nc_funs('attput', fname, nc_global, 'Additional_Documentation_2', '');
	nc_funs('attput', fname, nc_global, 'Additional_Documentation_3', '');
	nc_funs('attput', fname, nc_global, 'Additional_Documentation_4', '');
	nc_funs('attput', fname, nc_global, 'Additional_Documentation_5', '');
	nc_funs('attput', fname, nc_global, 'Additional_Documentation_6', '');
	nc_funs('attput', fname, nc_global, 'Additional_Documentation_7', '')
	nc_funs('attput', fname, nc_global, 'Source_Software', 'Mirone' )
	% -------------------------------------------------------------------------------------

	nons = int8(-128);		% A kind of nans

 	% ------------------------------ Write the variables ------------------------------------
	write_var(fname, 'time', 6, 'time', 'Time', 'seconds since 1970-01-01 00:00:00 0', [track.time(1) track.time(end)], 'UTC time, subtract TZ to get ship local time', [], [], [])
 	nc_funs('varput', fname, 'time', track.time, 0, numel(track.time));

	write_var(fname, 'drt', 1, [], 'Data Record Type', '', [], 'Normally 5', nons,  nons, [])
	nc_funs('varput', fname, 'drt', -128);
	write_var(fname, 'tz',  1, [], 'Time Zone Correction', 'hours', [], '-13 to +12 inclusive', nons,  nons, [])
	nc_funs('varput', fname, 'tz', 0);

	write_var(fname, 'lat', 4, 'time', 'Latitude', 'degrees_north', [min_lat max_lat], 'Negative south of Equator', int32(-2147483648), int32(-2147483648), 1e-7)
	nc_funs('varput', fname, 'lat', track.latitude);		% Scaling and type conversion is done inside nc_funs()

	write_var(fname, 'lon', 4, 'time', 'Longitude', 'degrees_east', [min_lon max_lon], 'Negative west of Greenwich', int32(-2147483648), int32(-2147483648), 2e-7)
	nc_funs('varput', fname, 'lon', track.longitude);		% Scaling and type conversion is done inside nc_funs()

	write_var(fname, 'ptc', 1, [], 'Position Type Code', 'hours', [], 'Observed (1), Interpolated (3), or Unspecified (9)', nons,  nons, [])
	nc_funs('varput', fname, 'ptc', -128);
	write_var(fname, 'twt', 4, [], 'Bathymetry Two-Way Travel-Time', 'seconds', [], 'Corrected for transducer depth, etc.', int32(-2147483648), int32(-2147483648), 1e-8)	
	nc_funs('varput', fname, 'twt', -2147483648);

	% ------------- Write TOPO ---------------------------------------------------------
	if (~all(isnan(track.topography)))
		write_var(fname, 'depth', 4, 'time', 'Bathymetry Corrected Depth', 'meter', [min(track.topography) max(track.topography)], 'Corrected for sound velocity variations (if known)', int32(-2147483648), int32(-2147483648), 1e-005)
		nc_funs('varput', fname, 'depth', -track.topography);
	else
		write_var(fname, 'depth', 4, [], 'Bathymetry Corrected Depth', 'meter', [], 'Corrected for sound velocity variations (if known)', int32(-2147483648), int32(-2147483648), 1e-005)
		nc_funs('varput', fname, 'depth', -2147483648);
	end
	nc_funs('attput', fname, 'depth', 'positive', 'down');
	% ----------------------------------------------------------------------------------

	write_var(fname, 'bcc', 1, [], 'Bathymetry Correction Code', '', [], '', nons,  nons, [])	
	nc_funs('varput', fname, 'bcc', -128);
	write_var(fname, 'btc', 1, [], 'Bathymetry Type Code', '', [], 'Observed (1), Interpolated (3), or Unspecified (9)', nons,  nons, [])	
	nc_funs('varput', fname, 'btc', -128);

	% ------------- Write TOTAL FIELD ----------------------------------------------------------
	if (~is_anom && ~all(isnan(track.magnetics)))
		write_var(fname, 'mtf1', 4, 'time', 'Magnetics First Sensor Total Field', 'gamma', [], 'Leading sensor', int32(-2147483648), int32(-2147483648), 0.0001)
		nc_funs('varput', fname, 'mtf1', track.magnetics + F_offset);
	else
		write_var(fname, 'mtf1', 4, [], 'Magnetics First Sensor Total Field', 'gamma', [], 'Leading sensor', int32(-2147483648), int32(-2147483648), 0.0001)
		nc_funs('varput', fname, 'mtf1', -2147483648);
	end
	% ----------------------------------------------------------------------------------

	write_var(fname, 'mtf2', 4, [], 'Magnetics Second Sensor Total Field', 'gamma', [], 'Trailing sensor', int32(-2147483648), int32(-2147483648), 0.0001)	
	nc_funs('varput', fname, 'mtf2', -2147483648);

	% ------------- Write Anomaly ----------------------------------------------------------
	%track.magnetics = track.gravity;
	%track.gravity = NaN;
	%is_anom = true;
	if (is_anom && ~all(isnan(track.magnetics)))
		write_var(fname, 'mag', 3, 'time', 'Magnetics Residual Field', 'gamma', [], 'Corrected for reference field (see header)', int16(-32768), int16(-32768), 0.1)	
		nc_funs('varput', fname, 'mag', track.magnetics);
	else
		write_var(fname, 'mag', 3, [], 'Magnetics Residual Field', 'gamma', [], 'Corrected for reference field (see header)', int16(-32768), int16(-32768), 0.1)
		nc_funs('varput', fname, 'mag', -32768);
	end
	% ----------------------------------------------------------------------------------

	write_var(fname, 'msens', 1, [], 'Magnetics Sensor For Residual Field', '', [], 'Magnetic sensor used: 1, 2, or Unspecified (9)', nons,  nons, [])	
	nc_funs('varput', fname, 'msens', -128);
	write_var(fname, 'diur', 3, [], 'Magnetics Diurnal Correction', 'gamma', [], 'Already applied to data', int16(-32768), int16(-32768), 0.1)	
	nc_funs('varput', fname, 'diur', -32768);
	write_var(fname, 'msd', 3, [], 'Magnetics Sensor Depth or Altitude', 'meter', [], 'Positive below sealevel', int16(-32768), int16(-32768), [])	
	nc_funs('varput', fname, 'msd', -32768);
	write_var(fname, 'gobs', 4, [], 'Gravity Observed', 'mGal', [], 'Corrected for Eotvos, drift, and tares', int32(-2147483648), int32(-2147483648), 1e-005)
	nc_funs('varput', fname, 'gobs', -2147483648);
	nc_funs('attput', fname, 'gobs', 'add_offset', 980000);
	write_var(fname, 'eot', 3, [], 'Gravity Eotvos Correction', 'mGal', [], '7.5 V cos (lat) sin (azim) + 0.0042 V*V', int16(-32768), int16(-32768), 0.1)
	nc_funs('varput', fname, 'eot', -32768);
	
	% ------------- Write FAA ----------------------------------------------------------
	if (~all(isnan(track.gravity)))
		write_var(fname, 'faa', 3, 'time', 'Gravity Free-Air Anomaly', 'mGal', [min(track.gravity) max(track.gravity)], 'Observed - theoretical', int16(-32768), int16(-32768), 0.1)
		nc_funs('varput', fname, 'faa', track.gravity);			% Scaling and type conversion is done inside nc_funs()
	else
		write_var(fname, 'faa', 3, [], 'Gravity Free-Air Anomaly', 'mGal', [], 'Observed - theoretical', int16(-32768), int16(-32768), 0.1)
		nc_funs('varput', fname, 'faa', -32768);
	end
	% ----------------------------------------------------------------------------------
	
	write_var(fname, 'nqc', 1, [], 'Navigation Quality Code', '', [], 'Suspected by (5) source agency, (6) NGDC, or no problems found (9)',  nons,  nons, [])
	nc_funs('varput', fname, 'nqc', -128);
	write_var(fname, 'id', 1, 'id_dim', 'Survey ID', '', [], 'Identical to ID in header', nons,  nons, [])
	nc_funs('varput', fname, 'id', int8(att.survey_id));
	write_var(fname, 'sln', 1, 'sln_dim', 'Seismic Line Number', '', [], 'For cross-referencing with seismic data', nons,  nons, [])
	nc_funs('varput', fname, 'sln', int8('99999') );
	write_var(fname, 'sspn', 1, 'sspn_dim', 'Seismic Shot-Point Number', '', [], 'For cross-referencing with seismic data', nons,  nons, [])
	nc_funs('varput', fname, 'sspn', int8('999999') );

% ----------------------------------------------------------------------------------------------
function write_var(fname, name, tipo, dim, long_name, units, actual_range, comment, fillValue, missing_value, scale_factor)	
	varstruct.Name = name;
	varstruct.Nctype = tipo;
	if (~isempty(dim)),		varstruct.Dimension = {dim};	end
	nc_funs('addvar', fname, varstruct)
	nc_funs('attput', fname, varstruct.Name, 'long_name', long_name );
	nc_funs('attput', fname, varstruct.Name, 'units', units);
	nc_funs('attput', fname, varstruct.Name, 'comment', comment);
	if (~isempty(actual_range)),	nc_funs('attput', fname, varstruct.Name, 'actual_range', actual_range);	end
	if (~isempty(fillValue)),		nc_funs('attput', fname, varstruct.Name, '_FillValue', fillValue);		end
	if (~isempty(missing_value)),	nc_funs('attput', fname, varstruct.Name, 'missing_value', missing_value);	end
	if (~isempty(scale_factor)),	nc_funs('attput', fname, varstruct.Name, 'scale_factor', scale_factor);	end

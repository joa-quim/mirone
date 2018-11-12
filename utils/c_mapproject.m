function out = c_mapproject(data, varargin)
% Temporary function to easy up transition from GMT4 to GMT5.2

% $Id: c_mapproject.m 11303 2018-05-28 21:39:31Z Joaquim Luis $

	cmd = 'mapproject';
	for (k = 1:numel(varargin))
		cmd = sprintf('%s %s', cmd, varargin{k});
	end
	out = gmtmex(cmd, data);
	gmtmex('destroy')


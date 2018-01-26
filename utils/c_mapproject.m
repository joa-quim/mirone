function out = c_mapproject(data, varargin)
% Temporary function to easy up transition from GMT4 to GMT5.2

% $Id: c_mapproject.m 10230 2018-01-26 01:34:32Z j $

	cmd = 'mapproject';
	for (k = 1:numel(varargin))
		cmd = sprintf('%s %s', cmd, varargin{k});
	end
	out = gmtmex(cmd, data);
	gmtmex('destroy')


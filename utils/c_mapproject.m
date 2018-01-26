function out = c_mapproject(data, varargin)
% Temporary function to easy up transition from GMT4 to GMT5.2

% $Id$

	cmd = 'mapproject';
	for (k = 1:numel(varargin))
		cmd = sprintf('%s %s', cmd, varargin{k});
	end
	out = gmtmex(cmd, data);
	gmtmex('destroy')


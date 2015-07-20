function out = _mapproject(data, varargin)
% Temporary function to easy up transition from GMT4 to GMT5.2

% $Id$

	global gmt_ver
	
	if (gmt_ver == 4)
		out = mapproject_m(data, varargin{:});
	else
		cmd = 'mapproject';
		for (k = 1:numel(varargin))
			cmd = sprintf('%s %s', cmd, varargin{k});
		end
		gmtmex('create')
		out = gmtmex(cmd, data);
		gmtmex('destroy')
	end

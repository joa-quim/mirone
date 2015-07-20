function out = c_shoredump(varargin)
% Temporary function to easy up transition from GMT4 to GMT5.2

% $Id$

	global gmt_ver
	if (isempty(gmt_ver)),		gmt_ver = 4;	end		% For example, if calls do not come via mirone.m
	
	if (gmt_ver == 4)
		out = shoredump(varargin{:});
	else
		cmd = 'pscoast';
		no_W = false;
		for (k = 1:numel(varargin))
			cmd = sprintf('%s %s', cmd, varargin{k});
			if (varargin{k}(2) == 'N' || varargin{k}(2) == 'I')
				no_W = true;
			end
		end
		if (~no_W),		cmd = [cmd ' -W'];	end		% If no Rivers or Borders, than Coastlines (-W) 
		gmtmex('create')
		out = gmtmex(cmd, data);
		gmtmex('destroy')
	end

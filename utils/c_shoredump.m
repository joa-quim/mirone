function out = c_shoredump(varargin)
% Temporary function to easy up transition from GMT4 to GMT5.2

% $Id: c_shoredump.m 7925 2016-06-23 00:22:14Z j $

	global gmt_ver
	if (isempty(gmt_ver)),		gmt_ver = 4;	end		% For example, if calls do not come via mirone.m
	
	if (gmt_ver == 4)
		out = shoredump(varargin{:})';
	else
		cmd = 'pscoast -M';
		no_W = false;
		maxLat = 0;
		for (k = 1:numel(varargin))
			if (varargin{k}(2) == 'R')		% We must fck saddly check if -R spans on [0 360] and if yes acting accordingly.
				ind = strfind(varargin{k}, '/');
				maxLat = str2double(varargin{k}(ind(1)+1:ind(2)-1));
			end
			if (maxLat > 180)
				cmd = sprintf('%s %s --FORMAT_GEO_OUT=+D', cmd, varargin{k});
			else
				cmd = sprintf('%s %s', cmd, varargin{k});
			end
			if (varargin{k}(2) == 'N' || varargin{k}(2) == 'I')
				no_W = true;
			end
		end
		if (~no_W),		cmd = [cmd ' -W'];	end		% If no Rivers or Borders, than Coastlines (-W) 
		out = gmtmex(cmd);
	end

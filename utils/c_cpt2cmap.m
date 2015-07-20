function [cmap, range] = c_cpt2cmap(fname, varargin)
% Temporary function to easy up transition from GMT4 to GMT5.2

% $Id$

	global gmt_ver
	
	if (gmt_ver == 4)
		if (nargout == 1)
			cmap = cpt2cmap(['-C' fname], varargin{:});
		else
			[cmap, range] = cpt2cmap(['-C' fname], varargin{:});
		end
	else
		cmd = ['read -Tc ' fname ' '];
		for (k = 1:numel(varargin))
			cmd = sprintf('%s %s', cmd, varargin{k});
		end
		gmtmex('create')
		C = gmtmex(cmd);
		cmap = C.colormap;
		gmtmex('destroy')
		if (nargout == 2)
			range = C.range;
		end
	end

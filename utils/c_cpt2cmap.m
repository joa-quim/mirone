function [cmap, range] = c_cpt2cmap(fname, varargin)
% Temporary function to easy up transition from GMT4 to GMT5

% $Id: c_cpt2cmap.m 11303 2018-05-28 21:39:31Z Joaquim Luis $
	
	if (strcmp(fname(1:2), '-C')),	fname = fname(3:end);	end		% Strip the -C from name. It should only be added here

	cmd = ['read -Tc ' fname ' '];
	for (k = 1:numel(varargin))
		cmd = sprintf('%s %s', cmd, varargin{k});
	end
	C = gmtmex(cmd);
	cmap = C.colormap;
	gmtmex('destroy')
	if (nargout == 2)
		range = C.range;
	end


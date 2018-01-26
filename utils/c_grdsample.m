function [Zout, hdr] = c_grdsample(Zin, head, varargin)
% Temporary function to easy up transition from GMT4 to GMT5

% The tsu_funs still calls grdsample directly as a system call

% $Id$

	G = fill_grid_struct(Zin, head);
	cmd = 'grdsample -n+c';
	for (k = 1:numel(varargin))
		cmd = sprintf('%s %s', cmd, varargin{k});
	end
	Zout = gmtmex(cmd, G);
	gmtmex('destroy')
	if (nargout == 1)
		Zout = Zout.z;
	else
		hdr = [Zout.range Zout.registration Zout.inc];
		Zout = Zout.z;
	end

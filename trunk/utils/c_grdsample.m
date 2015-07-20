function Zout = c_grdsample(Zin, head, varargin)
% Temporary function to easy up transition from GMT4 to GMT5.2

% O tsu_funs still calls grdsample directly as a system call

% $Id$

	global gmt_ver
	
	if (gmt_ver == 4)
		Zout = grdsample_m(Zin, head, varargin{:});
	else
		G = fill_grid_struct(Zin, head);
		cmd = 'grdsample';
		for (k = 1:numel(varargin))
			cmd = sprintf('%s %s', cmd, varargin{k});
		end
		gmtmex('create')
		Zout = gmtmex(cmd, G);
		Zout = Zout.z;
		gmtmex('destroy')
	end

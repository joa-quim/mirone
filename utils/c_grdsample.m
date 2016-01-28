function Zout = c_grdsample(Zin, head, varargin)
% Temporary function to easy up transition from GMT4 to GMT5.2

% O tsu_funs still calls grdsample directly as a system call

% $Id: c_grdsample.m 4726 2015-07-20 23:11:57Z j $

	global gmt_ver
	if (isempty(gmt_ver)),		gmt_ver = 4;	end		% For example, if calls do not come via mirone.m
	
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

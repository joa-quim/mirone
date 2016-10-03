function Zout = c_grdtrend(Zin, head, varargin)
% Temporary function to easy up transition from GMT4 to GMT5.2

% $Id: c_grdtrend.m 9840 2016-10-03 01:18:59Z j $

	global gmt_ver
	if (isempty(gmt_ver)),		gmt_ver = 4;	end		% For example, if calls do not come via mirone.m
	
	if (gmt_ver == 4)
		Zout = grdtrend_m(Zin, head, varargin{:});
	else
		G = fill_grid_struct(Zin, head);
		cmd = 'grdtrend';
		for (k = 1:numel(varargin))
			cmd = sprintf('%s %s', cmd, varargin{k});
		end
		Zout = gmtmex(cmd, G);
		gmtmex('destroy')
		Zout = Zout.z;
	end

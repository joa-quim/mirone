function Zout = c_grdtrend(Zin, head, varargin)
% Temporary function to easy up transition from GMT4 to GMT5.2

% $Id: c_grdtrend.m 11303 2018-05-28 21:39:31Z Joaquim Luis $

	G = fill_grid_struct(Zin, head);
	cmd = 'grdtrend';
	for (k = 1:numel(varargin))
		cmd = sprintf('%s %s', cmd, varargin{k});
	end
	Zout = gmtmex(cmd, G);
	gmtmex('destroy')
	Zout = Zout.z;


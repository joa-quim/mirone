function [Zout, head] = _grdproject(Zin, head, varargin)
% Temporary function to easy up transition from GMT4 to GMT5.2

% $Id$

	global gmt_ver
	
	if (gmt_ver == 4)
		if (nargout == 1)
			Zout = grdproject_m(Zin, head, varargin{:});
		else
			[Zout, head] = grdproject_m(Zin, head, varargin{:});
		end
	else
		G = fill_grid_struct(Zin, head);
		cmd = 'grdproject';
		for (k = 1:numel(varargin))
			cmd = sprintf('%s %s', cmd, varargin{k});
		end
		gmtmex('create')
		Z = gmtmex(cmd, G);
		Zout = Z.z;
		gmtmex('destroy')
		if (nargout == 2)
			head = Z.hdr;
		end
	end

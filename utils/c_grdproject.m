function [Zout, head] = c_grdproject(Zin, head, varargin)
% Temporary function to easy up transition from GMT4 to GMT5.2

% $Id: c_grdproject.m 7928 2016-06-23 00:27:25Z j $

	global gmt_ver
	if (isempty(gmt_ver)),		gmt_ver = 4;	end		% For example, if calls do not come via mirone.m
	
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
		Z = gmtmex(cmd, G);
		Zout = Z.z;
		gmtmex('destroy')
		if (nargout == 2)
			head = Z.hdr;
		end
	end

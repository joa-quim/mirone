function [Zout, hdr] = c_grdsample(Zin, head, varargin)
% Temporary function to easy up transition from GMT4 to GMT5

% The tsu_funs still calls grdsample directly as a system call

% $Id: c_grdsample.m 7897 2016-05-17 20:46:24Z j $

	global gmt_ver
	if (isempty(gmt_ver)),		gmt_ver = 4;	end		% For example, if calls do not come via mirone.m
	
	if (gmt_ver == 4)
		Zout = grdsample_m(Zin, head, varargin{:});
		if (nargout == 2)
			warndlg('Requesting two outputs from GMT4 grdsample_m MEX is not supported. Expect ...','WarnError')
			hdr = [];
		end
	else
		G = fill_grid_struct(Zin, head);
		cmd = 'grdsample';
		for (k = 1:numel(varargin))
			cmd = sprintf('%s %s', cmd, varargin{k});
		end
		Zout = gmtmex(cmd, G);
		if (nargout == 1)
			Zout = Zout.z;
		else
			hdr = [Zout.range Zout.registration Zout.inc];
			Zout = Zout.z;
		end
	end

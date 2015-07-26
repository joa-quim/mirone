function out = c_gmtlist(fname, varargin)
% Temporary function to easy up transition from GMT4 to GMT5.2

% $Id$

	global gmt_ver
	if (isempty(gmt_ver)),		gmt_ver = 4;	end		% For example, if calls do not come via mirone.m
	
	if (gmt_ver == 4)
		out = gmtlist_m(fname, varargin{:});
	else
		cmd = ['x2sys_datalist -Tgmt ' fname];
		for (k = 1:numel(varargin))
			cmd = sprintf('%s %s', cmd, varargin{k});
		end
		gmtmex('create')
		out = gmtmex(cmd);
		gmtmex('destroy')
	end

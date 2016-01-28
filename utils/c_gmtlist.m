function out = c_gmtlist(fname, varargin)
% Temporary function to easy up transition from GMT4 to GMT5.2

% $Id: c_gmtlist.m 4763 2015-07-28 23:58:47Z j $

	global gmt_ver
	if (isempty(gmt_ver)),		gmt_ver = 4;	end		% For example, if calls do not come via mirone.m
	
	if (gmt_ver == 4)
		out = gmtlist_m(fname, varargin{:});
	else
		ind = strfind(fname, '.gmt');
		if (isempty(ind)),		fname = [fname '.gmt'];		end		% Here we need the extension
		cmd = ['x2sys_datalist -Tgmt ' fname];
		for (k = 1:numel(varargin))
			if (strcmp(varargin{k}, '-G')),		continue,	end		% no -G here
			if (strcmp(varargin{k}(1:2), '-F'))
				varargin{k} = '-Ftime,lon,lat,faa,mag,top';
			end
			cmd = sprintf('%s %s', cmd, varargin{k});
		end
		gmtmex('create')
		out = gmtmex(cmd);
		gmtmex('destroy')
	end

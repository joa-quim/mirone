function out = c_gmtlist(fname, varargin)
% Temporary function to easy up transition from GMT4 to GMT5.2

% $Id$

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
	out = gmtmex(cmd);
	gmtmex('destroy')


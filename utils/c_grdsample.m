function [Zout, hdr] = c_grdsample(Zin, head, varargin)
% Temporary function to easy up transition from GMT4 to GMT5

% The tsu_funs still calls grdsample directly as a system call

% $Id: c_grdsample.m 10327 2018-03-22 18:14:50Z j $

	if (head(5) == 0 && head(6) == 0)
		zMinMax = grdutils(Zin,'-L');
		head(5) = zMinMax(1);		head(6) = zMinMax(2);
	end

	G = fill_grid_struct(Zin, head);
	cmd = 'grdsample';
	for (k = 1:numel(varargin))
		cmd = sprintf('%s %s', cmd, varargin{k});
	end

	ind = strfind(cmd, '-n');
	if (isempty(ind))
		cmd = [cmd ' -n+c'];
	else
		[t, r] = strtok(cmd(ind(1):end));
		cmd = [cmd(1:ind(1)-1) t '+c' r];
	end

	Zout = gmtmex(cmd, G);
	gmtmex('destroy')
	if (nargout == 1)
		Zout = Zout.z;
	else
		hdr = [Zout.range Zout.registration Zout.inc];
		Zout = Zout.z;
	end

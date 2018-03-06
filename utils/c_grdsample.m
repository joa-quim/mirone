function [Zout, hdr] = c_grdsample(Zin, head, varargin)
% Temporary function to easy up transition from GMT4 to GMT5

% The tsu_funs still calls grdsample directly as a system call

% $Id$

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

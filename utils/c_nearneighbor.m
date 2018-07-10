function [out, hdr] = c_nearneighbor(data, varargin)
% Temporary function to easy up transition from GMT4 to GMT5.2

% A merda é que 'data' pode ser um Mx3 ou entao um Mx1 e os X,Y[,W] devem tar no varargin

% $Id: c_nearneighbor.m 11368 2018-07-10 11:20:42Z j $

	% nearneighbor in GMT5 is less elastic in terms of how to swallow the data array(s)
	if (size(data,2) == 3 || size(data,2) == 4 || isa(data,'char'))
		k0 = 1;
	elseif (numel(varargin) >= 2 && isnumeric(varargin{1}) && size(varargin{1},2) == 1 && ...
			isnumeric(varargin{2}) && size(varargin{2},2) == 1)
		if (numel(varargin) >= 3 && isnumeric(varargin{3}) && size(varargin{3},2) == 1)
			data = [data varargin{1} varargin{2} varargin{3}];		% Input is now Mx4
			k0 = 4;
		else
			data = [data varargin{1} varargin{2}];					% Input is now Mx3
			k0 = 3;
		end
	end
	cmd = 'nearneighbor';
	for (k = k0:numel(varargin))
		if (strcmp(varargin{k}, '-e')),		continue,	end			% -e option exists only in GMT4 mex
		cmd = sprintf('%s %s', cmd, varargin{k});
	end
	G = gmtmex(cmd, data);
	gmtmex('destroy')
	if (nargout == 1)
		out = G.z;
	else
		hdr = [G.range G.registration G.inc];
		out = G.z;
	end

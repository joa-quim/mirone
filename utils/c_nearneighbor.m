function [out, hdr] = c_nearneighbor(data, varargin)
% Temporary function to easy up transition from GMT4 to GMT5.2

% A merda é que 'data' pode ser um Mx3 ou entao um Mx1 e os X,Y[,W] devem tar no varargin

% $Id: c_nearneighbor.m 4733 2015-07-22 17:03:43Z j $

	global gmt_ver
	if (isempty(gmt_ver)),		gmt_ver = 4;	end		% For example, if calls do not come via mirone.m
	
	if (gmt_ver == 4)
		[out, hdr] = nearneighbor_m(data, varargin{:});
	else
		% nearneighbor in GMT5 is less elastic in terms of how to swallow the data array(s)
		if (size(data,2) == 3 || size(data,2) == 4)
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
			cmd = sprintf('%s %s', cmd, varargin{k});
		end
		gmtmex('create')
		G = gmtmex(cmd, data);
		out = G.z;
		hdr = G.hdr;
		gmtmex('destroy')
	end

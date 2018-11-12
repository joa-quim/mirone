function [mask, head, X, Y] = c_grdlandmask(varargin)
% Temporary function to easy up transition from GMT4 to GMT5.2

% $Id: c_grdlandmask.m 11303 2018-05-28 21:39:31Z Joaquim Luis $

	cmd = 'grdlandmask';
	for (k = 1:numel(varargin))
		cmd = sprintf('%s %s', cmd, varargin{k});
	end
	cmd = strrep(cmd, '-e','');		% This -e option is only for the GMT4 mex
	mask = gmtmex(cmd);
	gmtmex('destroy')
	if (nargout > 1)
		head = [mask.range mask.registration mask.inc];	X = mask.x;		Y = mask.y;
	end
	mask = mask.z;

function [X, Y, Z, head] = _grdread(fname, varargin)
% Temporary function to easy up transition from GMT4 to GMT5.2

% $Id$

	global gmt_ver
	
	if (gmt_ver == 4)
		[X, Y, Z, head] = grdread_m(fname, 'single', varargin{:});
	else
		gmtmex('create')
		Zout = gmtmex(['read -Tg' fname]);
		X = Zout.x;
		Y = Zout.y;
		Z = Zout.z;
		head = Zout.hdr;
		gmtmex('destroy')
		head(10) = 0;		% Old grdread_m kept trace of data type on tenth element to preserve int16 types on output
	end

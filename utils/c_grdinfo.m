function info = c_grdinfo(fname, opt)
% Temporary function to easy up transition from GMT4 to GMT5.2

% $Id: c_grdinfo.m 7953 2016-09-10 01:28:20Z j $

	global gmt_ver
	if (isempty(gmt_ver)),		gmt_ver = 4;	end		% For example, if calls do not come via mirone.m
	
	if (gmt_ver == 4)
		info = grdinfo_m(fname, opt);
	else
		flen = numel(fname);
		if (~strcmp(opt, 'hdr_struct'))
			info = gmtmex(['grdinfo -C ' fname]);
		else
			s = gmtmex(['grdinfo ' fname]);
			info.Title = s.text{1}(flen+9:end);
			info.Command = s.text{2}(flen+11:end);
			info.Remark = s.text{3}(flen+10:end);
			info.Registration = s.text{4}(flen+3:end);
			info.Format = s.text{5}(flen+3:end);
			info.X_info = s.text{6}(flen+3:end);
			info.Y_info = s.text{7}(flen+3:end);
			info.Z_info = s.text{8}(flen+3:end);
			info.Scale  = s.text{9}(flen+3:end);
			info.hdr    = gmtmex(['grdinfo -C ' fname]);	% This is now a struct with 'data','text','header','comment','proj4','wkt'
		end
	end

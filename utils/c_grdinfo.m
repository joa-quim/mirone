function info = c_grdinfo(fname, opt)
% Temporary function to easy up transition from GMT4 to GMT5.2

% $Id: c_grdinfo.m 7737 2016-01-05 14:29:51Z j $

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
			info.Title = s{1}(flen+9:end);
			info.Command = s{2}(flen+11:end);
			info.Remark = s{3}(flen+10:end);
			info.Registration = s{4}(flen+3:end);
			info.Format = s{5}(flen+3:end);
			info.X_info = s{6}(flen+3:end);
			info.Y_info = s{7}(flen+3:end);
			info.Z_info = s{8}(flen+3:end);
			info.Scale  = s{9}(flen+3:end);
			info.hdr    = gmtmex(['grdinfo -C ' fname]);
		end
	end

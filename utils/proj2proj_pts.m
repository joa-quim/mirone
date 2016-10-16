function [out, msg, opt_R] = proj2proj_pts(handles, data, varargin)
% project DATA array to either the projection currently stored in appdata or selected in varargin
%
% This function is intended to replace GEOG2PROJECTED_PTS (and drop GMT proj) but we are only starting.
% 
% ...
% DATA can have > 3 columns but that's a bad idea

%	Copyright (c) 2004-2014 by J. Luis
%
% 	This program is part of Mirone and is free software; you can redistribute
% 	it and/or modify it under the terms of the GNU Lesser General Public
% 	License as published by the Free Software Foundation; either
% 	version 2.1 of the License, or any later version.
% 
% 	This program is distributed in the hope that it will be useful,
% 	but WITHOUT ANY WARRANTY; without even the implied warranty of
% 	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% 	Lesser General Public License for more details.
%
%	Contact info: w3.ualg.pt/~jluis/mirone
% --------------------------------------------------------------------

% $Id: proj2proj_pts.m 9866 2016-10-16 18:37:47Z j $

	opt_R = [];
	[srcWKT, srcProj4, srcGMT, dstProj4, dstWKT, limits, pad, msg] = parseIn(handles, varargin{:});
	%if (~isempty(msg)),		return,		end

	if (isempty(dstProj4) && isempty(dstWKT))		% Seek proj info in figure's appdata
		dstWKT = getappdata(handles.figure1,'ProjWKT');
		if (isempty(dstWKT))
			proj4 = getappdata(handles.figure1,'Proj4');
			if (~isempty(proj4))
				dstWKT = ogrproj(proj4);
			end
		end
		% If 'dstWKT' is still empty it means no proj info stored in image. Try 'geog'
		if (isempty(dstWKT) && handles.geog)
			if (~isequal(srcProj4, '+proj=longlat'))	% At least do not try to proj geog -> geog
				dstWKT = ogrproj('+proj=longlat');
			end
		end
	elseif ( isempty(dstWKT) && ~isempty(dstProj4) )
		dstWKT = ogrproj(dstProj4);
	end

	if (isempty(dstWKT))		% NOTHING to do here. Return whatever we got in except if we have a limits request
		if (~isempty(limits))
			data = trim_lims(data, limits, pad);
		end
		out = data;
		if (~isequal(srcProj4, '+proj=longlat'))	% Above we checked that if not the case, a dstWKT was assigned
			msg = 'unknown conversion';
		end
		return
	end

	if (isempty(srcWKT) && ~isempty(srcProj4))
		if (~strcmp(srcProj4, '+proj4=latlong'))	% Avoid a BUG somewhere as this proj4 string causes the error:
			projStruc.SrcProjSRS = srcProj4;		% OGRPROJ: Translating source SRS failed.
		end											% And by not using it, internally Geog WGS84 is assumed.
	elseif (~isempty(srcWKT))
		projStruc.SrcProjWKT = srcWKT;
	else
		warning('proj2proj_pts:zz','No source referencing system in input data. Returning as is.')		% The GMT case
		out = data;
		return
	end
	projStruc.DstProjWKT = dstWKT;

	if (~isa(data,'double'))		% Ghrrrr! I must make ogrproj work also with singles
		data = double(data);
	end

	if (size(data,2) > 3)			% If we have more than 3 columns send only the first 3
		out = ogrproj(data(:,1:3), projStruc);
		out = [out data(:,4:end)];		% And rebuild
	else
		out = ogrproj(data, projStruc);
	end
	
	if (~isempty(limits))
		out = trim_lims(out, limits, pad);
	end
	
	if (nargout == 3)
		x_min = min(out(:,1));		x_max = max(out(:,1));
		y_min = min(out(:,2));		y_max = max(out(:,2));
		opt_R = sprintf('-R%f/%f/%f/%f',x_min, x_max, y_min, y_max);
	end


% ------------------------------------------------------------------------------------------
function [srcWKT, srcProj4, srcGMT, dstProj4, dstWKT, limits, pad, msg] = parseIn(handles, varargin)
% ...
	msg = [];		limits = [];	pad = 0;
	srcWKT = [];	dstWKT = [];	srcProj4 = [];	srcGMT = [];	dstProj4 = [];

	for (k = 1:2:numel(varargin))
		if (strcmp(varargin{k}, 'srcWKT'))
			srcWKT = varargin{k+1};
		elseif (strcmp(varargin{k}, 'srcProj4'))
			srcProj4 = varargin{k+1};
		elseif (strcmp(varargin{k}, 'srcGMT'))
			srcGMT = varargin{k+1};
		elseif (strcmp(varargin{k}, 'dstProj4'))
			dstProj4 = varargin{k+1};
		elseif (strcmp(varargin{k}, 'dstWKT'))
			dstWKT = varargin{k+1};
		elseif (strncmp(varargin{k}, 'lim', 3))
			limits = varargin{k+1};
			if (numel(limits) < 4)
				msg = 'geog2projected_pts: ERROR, "lims" vector must have 4 elements';
				return
			end
		elseif (strcmp(varargin{k}, 'pad'))
			pad = varargin{k+1};
		end
	end

% ------------------------------------------------------------------------------------------
function data = trim_lims(data, lims, pad)
% Get rid of data that are outside the map limits
	if (pad ~= 0)
		lims = lims + [-pad pad -pad pad];      % Extend limits by PAD value
	end
	ind = (data(:,1) < lims(1) | data(:,1) > lims(2));
	data(ind,:) = [];
	ind = (data(:,2) < lims(3) | data(:,2) > lims(4));
	data(ind,:) = [];


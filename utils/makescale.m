function h = makescale(hAx, hObject)
%MAKESCALE plot a scale line.
%
%   MAKESCALE(HAX) creates a scale on the axis specificed by the handle of the
%       Mirone HAX axes and is based on its limits. The units will be either km,
%       m, or Nm in case that unit is slected in 'properties' and map is in geogs
%
%   H = MAKESCALE(...) returns a 2x1 vector with the handles of the line and text.
%
%   Note: This function should work for other terrestial bodies as well
%       A spherical body is assumed but it's the major axis that is used
%       in computations. The scale is correct at the point of insertion.

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

% $Id$

	handles = guidata(hAx);
	if (handles.no_file),		return,		end

	if strcmp(get(hObject,'checked'),'on')
		oldScale = findobj(hAx,'Tag','MapScale');
		if (~isempty(oldScale)),	delete(oldScale),	end
		set(hObject,'checked','off')
		return
	else
		set(hObject,'checked','on')
	end

	img_lims = getappdata(handles.axes1,'ThisImageLims');
	x_lim = img_lims(1:2);
	y_lim = img_lims(3:4);

	DY = diff(y_lim);
	x0 = x_lim(1) + 0.02 * diff(x_lim);
	y0 = y_lim(1) + 0.02 * DY;
	[slon,slat,str] = scale_line(handles, x0, y0, x_lim);
	
	% Make the scale
	hline = line('xdata',slon, 'ydata',slat, 'Parent',hAx, 'color','k', 'LineWidth',3,'Tag','MapScale');
	htext = text(slon(end),slat(1)+.01*DY,str,'HorizontalAlignment','center', ...
		'VerticalAlignment','bottom','FontSize',12,'FontWeight','bold','Tag','MapScale');

	if (nargout > 0),	h = [hline htext];	end

% -----------------------------------------------------------------------
function [x,y,str] = scale_line(handles, x0, y0, x_lim)
% Compute a scale ruler

	if (handles.geog)
		km_per_deg = pi * handles.DefineEllipsoide(1) / 180 * cos(y0*pi/180) * 1e-3;	% km / deg
		f = 1 / km_per_deg;
		ruler_len = dtick(diff(x_lim) * km_per_deg, true);
		str = sprintf('%g km', ruler_len);
		if (handles.DefineMeasureUnit(1) == 'n')	% If default unit is Nautical miles
			ruler_len = ruler_len / 2;
			str = sprintf('%g Nm', ruler_len);
		else
			if (ruler_len < 1)	% We are in m
				f = f * 1000;
				str = sprintf('%g m', ruler_len * 1000);
			end
		end
	else
		ruler_len = dtick(diff(x_lim) * 1e-3, false);
		str = sprintf('%g km', ruler_len);
		if (ruler_len < 1)		% We are in m
			f = 1000;
			str = sprintf('%g m', ruler_len * 1000);
		end
	end

	x = x0 + ruler_len * f;
	x = [x0 x];		y = [y0 y0];

% -----------------------------------------------------------------------
function dd = dtick(dlim,geog)
% From François Beauducel's function of the same name in his DEM program (BSD Licensed)

	if (geog && dlim <= 2/60)		% less than 2 minutes: base 36
		m = 10^floor(log10(dlim*36))/36;
	elseif (geog && dlim <= 2)		% less than 2 degrees: base 6
		m = 10^floor(log10(dlim*6))/6;
	else		% more than few degrees or not degrees: decimal rules
		m = 10^floor(log10(dlim));
	end
	p = ceil(dlim/m);
	if     (p <= 1),	dd = 0.1*m;
	elseif (p == 2)		dd = 0.2*m;
	elseif (p <= 5)		dd = 0.5*m;
	else				dd = m;
	end

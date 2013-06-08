function save_seismicity(h_mir_fig, h_events, opt)
% This function saves a seismicity object present in a Mirone figure,
% H_MIR_FIG -> handle to the Mirone figure containing the:
% H_EVENTS  -> handle to the seismicity object (The data is fished from this object handles)
%   OR
% H_EVENTS  -> [5 x n] double array with rows holding [x; y; events_time; events_mag; events_dep]
% OPT ... handle to a polygon (for the moment). In this case H_EVENTS must be []
% WARNING: The memory monster is at the wild

%	Copyright (c) 2004-2012 by J. Luis
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

% ordem = save_seismic_fields('write');   % Ask the order of fields for saving
% if (isempty(ordem))     return;     end
% ind = find(ordem == 0);                 % See if there are un-wanted fields
% if (isempty(ordem))    ordem(ind) = []; end
% fld = {'lon' 'lat' 'dep' 'mag', 'year' 'month' 'day' 'hour' 'minut' 'second'};
% for (k=1:length(ordem))
%     for (m=1:10)
%        if (ordem(k) == m)
%            tmp{k} = fld{m};
%            continue
%        end
%     end
% end

	handles_mir = guidata(h_mir_fig);       % Get the Mirone handles structure

	[FileName,PathName] = put_or_get_file(handles_mir, ...
		{'*.dat;*.DAT', 'Seimicity file (*.dat,*.DAT)'; '*.*', 'All Files (*.*)'},'Select File name','put','.dat');
	if isequal(FileName,0),		return,		end
	fname = [PathName FileName];

	if (~ishandle(h_events))				% Data was transmited in input
		x = h_events(1,:);    y = h_events(2,:);	% I hope that what guys from TMW say is true
		events_time = h_events(3,:);				% and this doesn't use any extra memory
		events_mag = h_events(4,:);
		events_dep = h_events(5,:);
	else									% We have to fish the data from appdata
        % Note: We have to tranpose those for fwrite (unbelivable stupid memory consuming fwrite)
        if (isempty(h_events))				% For polygons we don't know yet the seismicity handle
			h_events = findobj(h_mir_fig,'Tag','Earthquakes');
        end
        events_time = (getappdata(h_events,'SeismicityTime'))';
        events_mag  = (double(getappdata(h_events,'SeismicityMag')) / 10)';
        events_dep  = (double(getappdata(h_events,'SeismicityDepth')) / 10)';
        x = (get(h_events,'XData'));       y = (get(h_events,'YData'));
        if (nargin == 3 && ishandle(opt))		% OPT must be a handle to a polygon
			IN = inpolygon(x,y,get(opt,'XData'),get(opt,'YData'));
			if (any(IN))
				x(~IN) = [];            y(~IN) = [];    events_time(~IN) = [];
				events_mag(~IN) = [];   events_dep(~IN) = [];
			else
				warndlg('I don''t know about your eyes, but I don''t find any events inside the polygon.','Chico Clever')
				return
			end
		end
	end

	fid = fopen(fname, 'w');
	if (fid < 0),    errordlg(['Can''t open file:  ' fname],'Error');    return;     end

	year = fix(events_time);            % integer part
	frac = events_time - year;          % fraction part
	jdYear = date2jd(year);             % True Julian day of the 1st January of YEAR
	frac = frac .* (365 + isleapyear(year));
	[dumb,month, day, hour, minute] = jd2date((jdYear+frac),[]);
	fprintf(fid,'%.3f\t%.3f\t%d\t%02d\t%02d\t%.1f\t%.1f\t%02d\t%02d\n',...
		[x; y; year; month; day; events_mag; events_dep; hour; minute]);
	fclose(fid);

%--------------------------------------------------------------------------
function t = isleapyear(year)
	t = ( ~rem(year, 4) & rem(year, 100) ) | ~rem(year, 400);

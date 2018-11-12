function external_drive(handles, prog_name, varargin)
% Programatic interface to call sub-GUIs and run their callbacks

%	Copyright (c) 2004-2018 by J. Luis
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

% $Id: external_drive.m 11303 2018-05-28 21:39:31Z Joaquim Luis $

% ex: mirone('v:\pes.grd','-Cgrdsample_mir,guidata(gcf)','-Xcheckbox_Option_Q,1','-Xedit_x_inc,+0.04','-Xedit_y_inc,+0.04','-Xpush_OK')

try
	varargin = varargin{1};			% Because varargin was actually a cell of cells
	show_error = false;				% Because some errors are expected (when callbacks do not exist)
	n_arg = numel(varargin);
	for (k = 1:n_arg)
		if (isempty(varargin{k})),	continue,	end		% Shouldn't but hard to avoid
		if (strcmp(varargin{k}, 'debug'))
			show_error = true;
			continue
		end

		if (n_arg > 1 && ~isa(varargin{2}, 'char'))		% WTF case ???
			try
				feval(prog_name, [varargin{1} '_CB'], handles.(varargin{1}), handles, varargin{2:end});
			catch
				if (show_error)
					disp(['Errored at argument: ' varargin{k}])
					disp(lasterr)
				end
			end
			break
		end

		ind = strfind(varargin{k}, ',');
		% If it has a ",n" suffix take it to mean a set(...,'Val',n)
		if (isempty(ind))					% For example when running pushbutton (e.g. the "OK" button)
			try
				feval(prog_name, [varargin{k} '_CB'], handles.(varargin{k}), handles);
				handles = guidata(handles.figure1);
			catch
				if (show_error)
					disp(['Errored at argument: ' varargin{k}])
					disp(lasterr)
				end
			end
		else
			% If it has a ",n" suffix take it to mean a set(...,'Val',n), but if it's a
			% ",+<...> assume that we are changing an edit box
			butt = varargin{k}(1:ind-1);
			if (varargin{k}(ind+1) == '+')
				set(handles.(butt), 'Str', varargin{k}(ind+2:end))
			else
				set(handles.(butt), 'Val', str2double(varargin{k}(ind+1:end)))
			end
			try					% Because there is no guaranty that the callback even exists
				feval(prog_name, [butt '_CB'], handles.(butt), handles);
				handles = guidata(handles.figure1);
			catch
				if (show_error)
					disp(['Errored at argument: ' varargin{k}])
					disp(lasterr)
				end
			end
		end
	end
catch
	disp(lasterr)
end
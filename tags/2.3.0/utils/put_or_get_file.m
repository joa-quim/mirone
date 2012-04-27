function [FileName,PathName,handles] = put_or_get_file(handles,str1,str2,type, ext)
% Use this function to select input or output filename
% EXT if provided forces 'FileName' to output with that wished extension

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

	%return_to = handles.home_dir;		% Old behavior
	return_to = cd;						% New behavior. Return to where it was.
	if (strcmp(type,'get'))
		cd(handles.last_dir)
		[FileName,PathName] = uigetfile(str1,str2);
	elseif (strcmp(type,'put'))
		cd(handles.work_dir)
		[FileName,PathName] = uiputfile(str1,str2);
		if (isequal(FileName,0))
			cd(return_to)
			return
		end
		if (nargin == 5)		% Check that 'FileName' goes with the desired extension
			[PATH,FNAME,EXT] = fileparts(FileName);
			if (ext(1) ~= '.'),		ext = ['.' ext];		end
			if (isempty(EXT)),		FileName = [FNAME ext];	end
		end
	end
	pause(0.01);
	
	if (PathName ~= 0)
		handles.last_dir = PathName;
		if (nargout < 3),	guidata(handles.figure1,handles);	end
	end
	cd(return_to)

function [FileName,PathName,handles] = put_or_get_file(handles,str1,str2,type)
	% Use this function to select input or output filename
	if (strcmp(type,'get'))
		cd(handles.last_dir)
		[FileName,PathName] = uigetfile(str1,str2);
	elseif (strcmp(type,'put'))
		cd(handles.work_dir)
		[FileName,PathName] = uiputfile(str1,str2);
	end
	pause(0.01);
	
	if (PathName ~= 0)
		handles.last_dir = PathName;
		if (nargout < 3),	guidata(handles.figure1,handles);	end
	end
	cd(handles.home_dir);       % allways go home

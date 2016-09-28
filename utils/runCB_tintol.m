function runCB_tintol(handles)
% CallBack function that is executed at the end of mirone:show_image and hides
% all uicontrols not used in TINTOL.
% HANDLES is either the Mirone handles, or the Mirone Fig handle

	if (~isa(handles, 'struct'))		% Assume it's Mirone Fig handle and get handles
		handles = guidata(handles);
	end
	set([handles.Image handles.Tools handles.Plates handles.MagGrav ...
		handles.Seismology handles.GMT handles.GridTools], 'Vis', 'off')
% 	set([handles.DrawGeogCirc handles.Anaglyph handles.MBplaning], 'Vis', 'off')	% Stupid compiler screws here
 	set([handles.Anaglyph handles.MBplaning], 'Vis', 'off')
	
	set(findobj(handles.Tsunamis, 'Label', 'TINTOL'), 'Vis','off')			% When in TINTOL we don't want to see this anymore
	set(findobj(handles.Tsunamis, 'Label', 'Plot Stations'), 'Vis','on')	% And these were hiden waiting for showing up here
	set(findobj(handles.Tsunamis, 'Label', 'Import Stations'), 'Vis','on')
	set(findobj(handles.Tsunamis, 'Label', 'Swan'), 'Vis','off')
	
function varargout = cine_animations(varargin)
% CINE_ANIMATIONS

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
	
% $Id$

	if (isempty(varargin)),		return,		end
 
	handles = varargin{1};

	hUIs = cine_animations_LayoutFcn(handles.figure1);
	s = load([handles.path_data 'mirone_icons.mat'],'play_ico','pause_ico');
	img = cat(3, s.play_ico, s.pause_ico);
	set(hUIs(2),'CData',s.play_ico)
	set(hUIs(2),'UserData',img)					% We will store it here to toggle icons too
	setappdata(hUIs(1),'AnimaFamily', hUIs)		% 'resizetrue' will use them to set their units to
												% 'normalized'. We can't do it here because its too soon
	if (nargout),	varargout{1} = hUIs;	end

% ---------------------------------------------------------------------------------------
function slider_CB(hObject, handles, opt)
% 
	ud = get(handles.hUIcinanim(1),'UserData');
	if (get(hObject,'Min') == 0)		% First call. Time to set the limits that only now are known
		set(hObject,'Min',1, 'Max',ud(3), 'Val',1, 'SliderStep',[1 5]*(1 / (ud(3) - 1)))
	end

	if (nargin == 3)						% Called by toggle_playPause_CB
		set( hObject,'Value',ud(1), 'Tooltip',sprintf('%d',ud(1)) )
	else									% Manual frame setting
		val = round(get(hObject,'Value'));
		ind = val : val + ud(2) - 1;
		img = handles.cinemaImgs(:,:,ind);
		set(handles.hImg, 'CData',img)
		ud(1) = val + ud(2);
		if (ud(1) > ud(3)),		ud(1) = 1;		end
		set(handles.hUIcinanim(1),'UserData', ud)
		set(hObject,'Tooltip',sprintf('%d',val))
	end

% ---------------------------------------------------------------------------------------
function toggle_playPause_CB(hObject, handles)
% Toggle between pause and play at frame rates per second (so hardware allows)

	icons = get(handles.hUIcinanim(2),'UserData');
	non_stop = get(hObject,'Val');
	if (non_stop)
		set(hObject,'Tooltip','Pause')			% 
		set(handles.hUIcinanim(2), 'CData', icons(:,:,4:6))
	else
		set(handles.hUIcinanim(2), 'CData', icons(:,:,1:3))		% The 'Play' icon
	end
	fps = edit_frameRate_CB(handles.hUIcinanim(4), handles);	% Get frame rate
	fps_cum = 0;
	while (non_stop)
		ud = get(handles.hUIcinanim(1),'UserData');
		ind = ud(1) : ud(1) + ud(2) - 1;
		img = handles.cinemaImgs(:,:,ind);
		set(handles.hImg, 'CData',img)
		ud(1) = ud(1) + ud(2);
		fps_cum = fps_cum + fps;
		if (fps_cum > 0.5)						% Update slider only at >= 0.5 sec intervals
			slider_CB(handles.hUIcinanim(3), handles, ud(1))
			fps_cum = 0;
		end
		if (ud(1) > ud(3)),		ud(1) = 1;		end
		set(handles.hUIcinanim(1),'UserData', ud)
		pause(1 / fps)
		non_stop = get(hObject,'Val');			% This allows us to stop the weel on a second button press
	end

	% If it comes here it's because the turning round was interrupted
	set(hObject,'Tooltip','Play')

% ---------------------------------------------------------------------------------------
function fps = edit_frameRate_CB(hObject, handles)
% Just make sure that user does not do very silly things.
% There is a small complication that the number has also the '(fps)' text string
% that needs to be stripped of before converting to numeric

	xx = double(get(hObject,'Str'));
	ind = (xx >= 48 & xx <= 57);		% Numbers are ascii [40:57] ('012...9')
	xx = str2double( char(xx(ind)) );
	if (isnan(xx) || xx <= 0)
		xx = 5;		set(hObject,'Str','5 (fps)')
	elseif (xx > 24)					% Should never occur
		xx = 24;	set(hObject,'Str','24 (fps)')		% mpeg norm
	end
	if (nargout),	fps = xx;	end

% ---------------------------------------------------------------------
function h = cine_animations_LayoutFcn(hMirFig)

	frameWidth = 290;
	posFig = get(hMirFig,'Pos');
	
	x0 = posFig(3) - frameWidth + 1;

	h(1) = uicontrol('Parent',hMirFig, 'Position',[x0 0 frameWidth 22], 'Style','frame','Tag','cinanima');

	h(2) = uicontrol('Parent',hMirFig, 'Position',[x0+5 0 21 21],...
	'Call',{@cinanim_CB, hMirFig},...
	'Style','togglebutton',...
	'Tooltip','Play',...
	'SelectionHighlight','off',...
	'Tag','toggle_playPause');

	h(3) = uicontrol('Parent',hMirFig, 'Position',[x0+25 2 211 16],...
	'Call',{@cinanim_CB, hMirFig},...
	'BackgroundColor',[0.9 0.9 0.9],...
	'Style','slider',...
	'Min',0,...
	'Max',1,...			% This will set to the case dependent correct value latter
	'Value',0,...
	'Tag','slider');

	h(4) = uicontrol('Parent',hMirFig, 'Position',[x0+237 1 52 21],...
	'Call',{@cinanim_CB, hMirFig},...
	'BackgroundColor',[1 1 1],...
	'String','5 (fps)',...
	'Style','edit',...
	'Tooltip','Frames per second (decimals allowed)',...
	'Tag','edit_frameRate');


function cinanim_CB(hObject, evt, h1)
% This function is executed by the Call and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(h1));

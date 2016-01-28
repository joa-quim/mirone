function show_manguito(varargin)
% Show a user stupid error inspiring image for a brieve moment

%	Copyright (c) 2004-2013 by J. Luis
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

% $Id: show_manguito.m 3885 2013-02-26 00:34:27Z j $

	hObject = figure('Vis','off');
	show_manguito_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(hObject,'center')

	if (nargin == 0)
		mir_dirs = getappdata(0,'MIRONE_DIRS');
		if (~isempty(mir_dirs))
			home_dir = mir_dirs.home_dir;		% Start in values
		else
			home_dir = cd;
		end
		manguito = [home_dir '/data/manguito_galactico.jpg'];
	else
		manguito = [varargin{1}.path_data 'manguito_galactico.jpg'];
	end
	
	% Load image
	manguito = imread(manguito);
	image(manguito,'Parent',handles.axes1)
	set(handles.axes1,'Visible', 'off');
	set(hObject,'Vis','on');

	if (strncmp(computer,'PC',2))
		WindowAPI(hObject, 'Clip', true)	% Clip borders. splash screen
	end

	if (nargin == 2)
		pause(varargin{2})		% User transmitted pause time
	else
		pause(1)
	end

	delete(hObject)

% ---------------------------------- 
function show_manguito_LayoutFcn(h1)

set(h1, 'Position',[256 432 250 250],...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'MenuBar','none',...
'Name','Manguito Galactico',...
'NumberTitle','off',...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1');

axes('Parent',h1,...
'Units','pixels',...
'Position',[1 1 250 250],...
'CameraPosition',[0.5 0.5 9.16025403784439],...
'Tag','axes1');

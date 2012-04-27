function [X,Y,Z,head,m,n] = load_grd(handles, opt)
% Load a grid either from memory, or re-read it (GMT only) if it's too big to fit in it (biger than handles.grdMaxSize)
%
% OPT --> 'silent'		Shut up, even in case of an error message
% OPT --> 'double'		Resturn Z as a double (default is Z original type)

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

	% Fake image with bg_color. Nothing loadable
	if (handles.image_type == 20),  X=[];	Y=[];	Z=[];	head=[];   return,	end     
	
	if (nargin == 1),	opt = '';	end
	X = getappdata(handles.figure1,'dem_x');	Y = getappdata(handles.figure1,'dem_y');
	Z = getappdata(handles.figure1,'dem_z');	head = handles.head;
	if (nargout > 4),		[m,n] = size(Z);	end
	err_msg = '';

	if (isempty(X) && handles.image_type == 1 && ~handles.computed_grid)
		[handles, X, Y, Z, head] = read_gmt_type_grids(handles,handles.grdname);
		if (nargout == 6),	[m,n] = size(Z);	end
	elseif (isempty(X) && handles.image_type == 4 && ~handles.computed_grid)
		Z = [];			% Check for this to detect a (re)-loading error
		err_msg = 'Grid was not on memory. Increase "Grid max size" and start over again.';
	elseif (isempty(X) && isempty(handles.grdname))
		Z = [];			% Check for this to detect a (re)-loading error
		err_msg = 'Grid could not be reloaded. You probably need to increase "Grid max size"';
	end
	if (strncmp(opt,'double',1) && ~isa(Z,'double')),	Z = double(Z);		end
	if (~isempty(err_msg) && ~strcmp(opt,'silent')),	errordlg(err_msg,'ERROR'),		end

function InOut2WS(handles, opt)
% OPT = 'direct' => copy X,Y,Z,I,head to the base workspace
% OPT = 'GRID_inverse' => import X,Y,Z,head from the base workspace
% OPT = 'IMG_inverse' => import X,Y,I from the base workspace
% OPT = 'loadmat' => import a .mat file with a 2D array named Z or z
% OPT = 'clear' => clear the base workspace

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

if ~strcmp(opt,'loadmat')
	if (handles.no_file),      return;      end
	is_grid = handles.validGrid;
	vars = evalin('base','who');
end

try     % So many thing may go wrong
	if (strcmp(opt,'direct') && is_grid)
		[X,Y,Z,head] = load_grd(handles);
		if isempty(Z),		return,		end			% An error message was already issued
		assignin('base','X',X);			assignin('base','Y',Y);
		assignin('base','Z',double(Z)); assignin('base','head',head);
		assignin('base','img',get(handles.hImg,'CData'));
		assignin('base','FigHand',handles.figure1);
		assignin('base','is_grid',is_grid);			% Tell base WS if we have a grid or a image
	elseif (strcmp(opt,'direct') && ~is_grid)
		assignin('base','img',flipdim(get(handles.hImg,'CData'),1));
		[m,n,k] = size(get(handles.hImg,'CData'));
		assignin('base','is_grid',is_grid);         % Tell base WS if we have a grid or a image
		assignin('base','head',handles.head);
	elseif strcmp(opt,'GRID_inverse')
		if (~any(strcmp('is_grid',vars)))
			errordlg('You screwed up my memory if we are dealing with a grid or a image.','Error'),		return
		end
		is_grid = evalin('base','is_grid');
		if (is_grid == 0)
			errordlg('You are confused. There is no grid available here.','Error');     return
		end
		if (~any(strcmp('X',vars)) || ~any(strcmp('Y',vars)))
			errordlg('You screw up one off X or Y coordinate vectors. Restart or give up.','Error'),	return
		end
		Xb = evalin('base','X');            Yb = evalin('base','Y');
		if (~any(strcmp('Z',vars)) || ~any(strcmp('head',vars)))
			errordlg('You screw up either the Z or head arrays. Restart or give up.','Error');  return
		end
		Zb = single(evalin('base','Z'));    head_b = evalin('base','head');
		if (size(Zb,1) ~= length(Yb) || size(Zb,2) ~= length(Xb))
			errordlg('Z array size is not compatible with the X & Y vectors. Restart or give up.','Error'),	return
		elseif (length(head_b) < 9)
			errordlg('You screw up the header variable. Rebuild it or give up.','Error'),	return
		end
		[X,Y,Z,head_b] = load_grd(handles);
		if (isequal(Z,Zb)),     return;     end     % Z did not change. Return
		% OK, if we get here, its time to use the new grid's values
		zzz = grdutils(Zb,'-L');  z_min = zzz(1);     z_max = zzz(2);     clear zzz;
		head_b(5:6) = [z_min z_max];
		tmp.X = Xb;    tmp.Y = Yb;    tmp.head = head_b;    tmp.name = 'ML computed grid';
		mirone(Zb,tmp);
	elseif strcmp(opt,'IMG_inverse')
		if (~any(strcmp('img',vars)))
			errordlg('You screw up the IMG array var. Restart or give up.','Error'),	return
		end
		I = evalin('base','img');
		if ( handles.image_type == 2)
			mirone(I)
		else
			if (~any(strcmp('head',vars)))
				errordlg('You screw up the necessary header vector. Restart or give up.','Error'),	return
			end
			head = evalin('base','head');
			tmp.X = head(1:2);		tmp.Y = head(3:4);		tmp.head = head;
			if (~isa(I,'logical') && ndims(I) == 2),	tmp.cmap = get(handles.figure1,'Colormap'),		end
			tmp.name = 'ML imported img';
			mirone(I,tmp);
		end

	elseif strcmp(opt,'loadmat')            % Load a .mat file containing a 2D array
		[FileName,PathName] = put_or_get_file(handles, ...
			{'*.mat;','mat file format (*.mat)'; '*.*', 'All Files (*.*)'},'Select .mat file','get');
		if isequal(FileName,0),		return,		end

		try                                             % We do a try-catch because evalin does not work on the compiled version
			load([PathName FileName])
			if (~isa(Z,'single')),   Z = single(Z);      end
		catch
			try
				if (~isa(z,'single')),   Z = single(z);  clear z;    end
			catch
				errordlg('The "Open .mat 2D array" requires that the array variable is named "Z" or "z"','ERROR')
				return
			end
		end
		[m,n,k] = size(Z);
		if (m < 2 || n < 2)
			errordlg('You must be joking. This is a vector, not a true 2D array.','ERROR')
			return
		end
		if (ndims(Z) > 2),  Z = Z(:,:,1);    end
		[zz] = grdutils(Z,'-L');  z_min = zz(1);     z_max = zz(2);     clear zz;
		head = [1 n 1 m z_min z_max 0 1 1];     % Minimalist header
		tmp.X = 1:n;    tmp.Y = 1:m;    tmp.head = head;    tmp.name = 'ML imput grid';
		mirone(Z,tmp);
	else
		evalin('base','clear')
	end
catch
    errordlg(lasterr,'Error')
end


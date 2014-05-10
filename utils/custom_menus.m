function custom_menus(hFig, path_data)
% Creates Menus from data read from the 'Custom_menu_def.txt' file
%
%	This functions is called by mirone at each Mirone figure creation time. Than:
%		hFig		handle of the mirone fig
%		path_data	Mirone data dir
%
%	The Custom_menu_def.txt has a form of this type (example)
%		Geography,Basins[/sub],dir

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

	fname = [path_data 'Custom_menu_def.txt'];
	if (exist(fname,'file') ~= 2),	return,		end
	fid = fopen(fname);
	if (fid < 0),	return,		end

	c = (fread(fid,'*char'))';      fclose(fid);
	menus = strread(c,'%s','delimiter','\n');   clear c fid;
	% Remove eventual empty lines
	m = numel(menus);    c = false(m,1);
	for (k = 1:m)
		if (isempty(menus{k})),     c(k) = true;    end
	end
	menus(c) = [];

	%Geography,Basins[/sub],dir

	% Parse the file contents
	m = numel(menus);    c = false(m,1);
	mainMenu = cell(m,1);
	subMenu1  = cell(m,1);
	subMenu2 = cell(m,1);
	dataDir  = cell(m,1);
	for (k = 1:m)
		if (menus{k}(1) == '#'),	c(k) = true;	continue,	end
		ind = strfind(menus{k}, ',');
		mainMenu{k} = menus{k}(1:ind(1)-1);
		subMenu1{k}  = menus{k}(ind(1)+1:ind(2)-1);	% This one must further be parsed for Basins[/sub]
		inds = strfind(subMenu1{k}, '/');
		if (~isempty(inds))
			subMenu2{k} = subMenu1{k}(inds(1)+1:end);
			subMenu1{k}(inds(1):end) = [];
		end
		dataDir{k} = menus{k}(ind(2)+1:end);
	end
	mainMenu(c) = [];	subMenu1(c) = [];	subMenu2(c) = [];	dataDir(c) = [];   % Remove comment lines
	m = numel(mainMenu);			% Update counting

	hSbub1 = zeros(1,m);    hSbub2  = zeros(1,m);
	n = 1;
	for (k = 1:m)
		hThis = findobj(hFig,'Tag',mainMenu{k});
		if (isempty(hThis))
			% Make a new one on the main toolbar. Allow it?
		end
		
		if (~hSbub1(n)),	hSbub1(n) = uimenu('Parent',hThis,'Label',['Custom -> ' subMenu1{n}]);	end
		if (k == 1),		set(hSbub1(1), 'Sep','on'),		end		% Only first has Separator
		if (k > 1 && ~strcmp(subMenu1{k}, subMenu1{k-1}))
			n = n + 1;
		end
		if (isempty(subMenu2{k}))
			set(hSbub1(k), 'Call',['hand=guidata(gcf);hand.last_dir=''',dataDir{k},''';mirone(''TransferB_CB'',hand,''guessType'')'])
		else
			hSbub2(k) = uimenu('Parent',hSbub1(n),'Label',subMenu2{k}, 'Call', ...
				['hand=guidata(gcf);hand.last_dir=''',dataDir{k},''';mirone(''TransferB_CB'',hand,''guessType'')']);
		end
	end

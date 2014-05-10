% --------------------------- PROJECTIONS MENU ------------------------------------
function projection_menu(hFig, hProj, home_dir)
% Creates the Projection Menu from data read from the 'SRSproj_def.txt' file
%
%	This functions is called by mirone_uis at each Mirone figure creation time. Than:
%		hFig		handle of the mirone fig
%		hProj		handle of "Projections" uimenu
%		home_dir	Mirone home dir
%
%	Alternatively it may be called at a later time but with 2 arguments ONLY. Than
%		hFig		The Mirone handles
%		hProj		A GMT -J type projection info (used by grdroject)

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

	if (nargin == 2)
		setPRJ([],[], hFig, 1, hProj)
		return
	end

	fid = fopen([home_dir filesep 'data' filesep 'SRSproj_def.txt']);
	if (fid < 0),	return,		end

	c = (fread(fid,'*char'))';      fclose(fid);
	menus = strread(c,'%s','delimiter','\n');   clear c fid;
	% Remove eventual empty lines
	m = numel(menus);    c = false(m,1);
	for (k = 1:m)
		if (isempty(menus{k})),     c(k) = true;    end
	end
	menus(c) = [];

	% Parse the file contents
	m = numel(menus);    c = false(m,1);
	mainMenu = cell(m,1);
	subMenu = cell(m,1);
	projStr = cell(m,1);
	for (k = 1:m)
		if (menus{k}(1) == '#'),	c(k) = true;	continue,	end
		[t,r] = getMenuLabel(menus{k});
		mainMenu{k} = t;
		[t,r] = getMenuLabel(r);
		if (numel(t) == 1),     subMenu{k} = [];    % We don't have a submenu
		else                    subMenu{k} = t;
		end
		if (~strcmp(mainMenu{k},'None')),   projStr{k} = ripblanks(r);      end
	end
	mainMenu(c) = [];    subMenu(c) = [];    projStr(c) = [];   % Remove comments lines
	m = numel(mainMenu);       % Update counting

	% Detect repeated main menus which occur when we have subMenus
	primos = true(m,1);
	for (k = 2:m),		primos(k) = ~strcmp(mainMenu{k}, mainMenu{k-1});    end

	% Parse the projStr string and assign it either to a GMT type -J proj or proj4 string
	projCOMM = cell(m,3);       % 3 for -J -C -T in case of a GMT proj 
	for (k = 1:m)
		if (isempty(projStr{k})),	continue,	end      % The 'None' line

		if (projStr{k}(1) == '+')   % We have a Proj4 string
			projCOMM{k,1} = projStr{k};
			continue
		end

		[t,r] = strtok(projStr{k});
		if ( strcmp(t(2:3),'Ju') && ~strcmp(t(end-1:end),'/1') )	% Risky but the Ju<zone>/1:1 crashes mapproject_m
			t = sprintf('%s/1', t);
		else
			if ( strcmp(t(end-1:end),'/1') )
				t = sprintf('%s:1', t);		% Append scale
			elseif ( ~strcmp(t(end-2:end), '1:1') )
				t = sprintf('%s/1:1', t);	% Append scale
			end
		end
		projCOMM{k,1} = t;
		if (~isempty(r))            % Either a -C<...> or -T
			j = 0;
			while (r)
				[t,r] = strtok(r);
				projCOMM{k,2+j} = t;
				j = j + 1;
			end
		else
			projCOMM{k,2} = '-C';    % A must have
		end
	end

	hMain = zeros(1,m);    hSec  = zeros(1,m);
	for (k = 1:m)
		if (isempty(subMenu{k}))    % No subMenu, we have than a direct projection setting
			hMain(k) = uimenu('Parent',hProj,'Label',mainMenu{k},'Call',{@setPRJ,hFig,k,projCOMM});
		else
			if (primos(k))          % Parent of a Submenu. No proj string settings
				hMain(k) = uimenu('Parent',hProj,'Label',mainMenu{k});
				hUI = hMain(k);     % Make a copy to use on SubMenus
				hSec(k) = uimenu('Parent',hUI,'Label',subMenu{k},'Call',{@setPRJ,hFig,k,projCOMM});
			else                    % Child of a Parent with a submenu. We have a proj string to set
				hSec(k) = uimenu('Parent',hUI,'Label',subMenu{k},'Call',{@setPRJ,hFig,k,projCOMM});
			end
		end
	end

	hMain(hMain == 0) = [];
	hSec(hSec == 0) = [];

	set(hMain(1), 'checked', 'on')
	setappdata(hFig, 'ProjList', [hMain hSec])

% -------------------------------------------------------------------------------
function [t,r] = getMenuLabel(s)
    % S contains the Menu label and proj string. Return the Label in T and the remaining in R 
    [t,r] = strtok(s);
    if (t(1) == '"')    % A composite name (for example "Portuguese Mess")
        ind = strfind(r,'"');
        t = [t(2:end) r(1:ind(1)-1)];
        r = r(ind(1)+1:end);
    end
        
% -------------------------------------------------------------------------------
function s = ripblanks(s)
% removes leading and trailing white space from S
	[r, c] = find(~isspace(s));
	s = s(min(c):max(c));

% -------------------------------------------------------------------------------
function setPRJ(obj,nikles, hFig, k, projCOMM)
% Set the projection string in Figure's appdata and a checkmark on the selected projection

	if (isempty(obj))		% 'External' call. Probably after a grid reprojected with grdproject
		if (~isa(projCOMM, 'cell'))		projCOMM = {projCOMM};		end
		if (numel(projCOMM) == 1)		projCOMM{end+1} = '-C';		end
		handles = hFig;
		hFig = handles.figure1;
	else
		handles = guidata(hFig);
		projList = getappdata(hFig,'ProjList');
	    set(projList,'checked','off');
		set(obj,'checked','on')
	end

	if (strcmp(get(obj,'Label'),'None'))
		setappdata(hFig,'ProjGMT','')
		setappdata(hFig,'Proj4','')
		handles.is_projected = 0;
	elseif (projCOMM{k}(1) == '+')
		setappdata(hFig,'Proj4',projCOMM{k,1})  % Here the columns 2 and 3 are empty
		setappdata(hFig,'ProjGMT','')
		handles.is_projected = 1;
	elseif (projCOMM{k}(1) == '-')
		prj = projCOMM(k,:);
		% If doesn't have the third element, remove it because no empties in mapproject
		if (isempty(prj{end})),     prj(end) = [];  end
		setappdata(hFig,'ProjGMT',prj)
		setappdata(hFig,'Proj4','')
		handles.is_projected = 1;
	else
		set(obj,'checked','off')
		setappdata(hFig,'ProjGMT','')
		setappdata(hFig,'Proj4','')
		handles.is_projected = 0;
		errordlg(['Bad projection string: ' projCOMM{k,:}],'Error')
	end
	guidata(hFig, handles)
	setAxesDefCoordIn(handles);     % Sets the value of the axes uicontextmenu that selects whether project or not
	aux_funs('toProjPT',handles)    % Set vars to be used when PT->geog conversion, used in PIXVAL_STSBAR
	if (isempty(getappdata(hFig,'DispInGeogs')))
		setappdata(hFig,'DispInGeogs',0)     % We need to set it now (this is the default)
	end
